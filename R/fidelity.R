#' Affine-Invariant Riemannian Distance Between SPD Matrices
#'
#' @description
#' Computes
#' \deqn{d_{\mathrm{AI}}(\mathbf{C}_1, \mathbf{C}_2) =
#'   \big\lVert \log\!\big( \mathbf{C}_1^{-1/2} \mathbf{C}_2 \mathbf{C}_1^{-1/2}
#'   \big) \big\rVert_F.}
#' Unlike \code{\link{cov_logeuclidean_distance}}, the affine-invariant
#' distance is the geodesic distance on the SPD manifold under the
#' canonical affine-invariant metric of Pennec et al.\ (2006). The two
#' distances agree at the identity but can diverge for matrices far from
#' the identity, and reporting both gives a multi-metric fidelity audit.
#'
#' @param C1,C2 Numeric SPD matrices of identical shape.
#' @param eigenvalue_floor Positive scalar. Floor used when computing
#'   the matrix inverse square root.
#' @return A non-negative scalar.
#' @references
#' Pennec, X., Fillard, P., Ayache, N. (2006). A Riemannian framework
#' for tensor computing. Int. J. Comput. Vis., 66(1), 41--66.
#' @export
cov_affine_invariant_distance <- function(C1, C2, eigenvalue_floor = 1e-12) {
  if (!is.matrix(C1) || !is.matrix(C2) || !all(dim(C1) == dim(C2))) {
    stop("C1 and C2 must be matrices of the same shape.")
  }
  S <- 0.5 * (C1 + t(C1))
  eig <- eigen(S, symmetric = TRUE)
  vals <- pmax(eig$values, eigenvalue_floor)
  inv_sqrt <- eig$vectors %*% diag(vals^(-0.5), nrow = length(vals)) %*%
    t(eig$vectors)
  middle <- inv_sqrt %*% (0.5 * (C2 + t(C2))) %*% inv_sqrt
  sqrt(sum(.spd_logm(middle)^2))
}


#' Pearson Correlation of Eigenvalue Spectra Between Two SPD Matrices
#'
#' @description
#' Returns the Pearson correlation of the eigenvalue spectra of \code{C1}
#' and \code{C2} after sorting in descending order. A value near 1
#' indicates that the augmented covariance preserves the spectral shape
#' of the real reference; values below 1 indicate spectral distortion.
#'
#' @param C1,C2 Numeric SPD matrices of identical shape.
#' @return A correlation in [-1, 1] (or \code{NA} if either spectrum
#'   has zero variance).
#' @importFrom stats cor
#' @export
cov_eigenvalue_correlation <- function(C1, C2) {
  if (!is.matrix(C1) || !is.matrix(C2) || !all(dim(C1) == dim(C2))) {
    stop("C1 and C2 must be matrices of the same shape.")
  }
  v1 <- sort(eigen(0.5 * (C1 + t(C1)), symmetric = TRUE,
                   only.values = TRUE)$values, decreasing = TRUE)
  v2 <- sort(eigen(0.5 * (C2 + t(C2)), symmetric = TRUE,
                   only.values = TRUE)$values, decreasing = TRUE)
  if (stats::sd(v1) == 0 || stats::sd(v2) == 0) return(NA_real_)
  stats::cor(v1, v2)
}


#' Trace Ratio Between Two Covariance Matrices
#'
#' @description
#' Computes \eqn{\mathrm{tr}(\mathbf{C}_{\text{aug}}) /
#'   \mathrm{tr}(\mathbf{C}_{\text{real}})}. The ideal value is 1
#' (total covariance power preserved); deviations measure power
#' inflation (>1) or attenuation (<1).
#'
#' @param C_aug Numeric SPD matrix.
#' @param C_real Numeric SPD matrix of identical shape.
#' @return A positive scalar.
#' @export
cov_trace_ratio <- function(C_aug, C_real) {
  if (!is.matrix(C_aug) || !is.matrix(C_real) ||
      !all(dim(C_aug) == dim(C_real))) {
    stop("C_aug and C_real must be matrices of the same shape.")
  }
  sum(diag(C_aug)) / max(sum(diag(C_real)), 1e-12)
}


#' Condition-Number Ratio Between Two Covariance Matrices
#'
#' @description
#' Computes the ratio of condition numbers \eqn{\kappa(\mathbf{C}_{\text{aug}})
#'   / \kappa(\mathbf{C}_{\text{real}})}. Off-manifold projection or
#' isotropic perturbation can compress or stretch the spectrum
#' asymmetrically; this ratio detects numerical/geometric instability
#' that the trace ratio alone may miss.
#'
#' @param C_aug Numeric SPD matrix.
#' @param C_real Numeric SPD matrix of identical shape.
#' @return A positive scalar.
#' @export
cov_condition_ratio <- function(C_aug, C_real) {
  if (!is.matrix(C_aug) || !is.matrix(C_real) ||
      !all(dim(C_aug) == dim(C_real))) {
    stop("C_aug and C_real must be matrices of the same shape.")
  }
  k_aug <- {
    e <- eigen(0.5 * (C_aug + t(C_aug)), symmetric = TRUE,
               only.values = TRUE)$values
    max(e) / max(min(e), 1e-12)
  }
  k_real <- {
    e <- eigen(0.5 * (C_real + t(C_real)), symmetric = TRUE,
               only.values = TRUE)$values
    max(e) / max(min(e), 1e-12)
  }
  k_aug / max(k_real, 1e-12)
}


#' Anchor-to-Synthetic Perturbation Distance
#'
#' @description
#' Mean log-Euclidean distance from each synthetic covariance to its
#' originating anchor. Used by the E0 amplitude-matching diagnostic to
#' verify that the off-manifold and on-manifold methods produce
#' comparable per-anchor perturbation magnitudes; a ratio
#' \eqn{\rho = d_{E0}/d_{G0}} close to 1 indicates the matching
#' contract is satisfied (Shen \& Degras, 2026).
#'
#' @param anchors List of SPD anchor covariance matrices.
#' @param synthetic List of SPD synthetic covariance matrices.
#' @param anchor_idx Integer vector of length \code{length(synthetic)}
#'   indexing into \code{anchors}.
#' @return A non-negative scalar.
#' @export
cov_anchor_perturbation_distance <- function(anchors, synthetic,
                                             anchor_idx) {
  if (!is.list(anchors) || !is.list(synthetic)) {
    stop("anchors and synthetic must be lists of SPD matrices.")
  }
  if (length(synthetic) != length(anchor_idx)) {
    stop("anchor_idx must have length equal to length(synthetic).")
  }
  if (length(synthetic) == 0L) return(NA_real_)
  anchor_log <- lapply(anchors, .spd_logm)
  d <- numeric(length(synthetic))
  for (k in seq_along(synthetic)) {
    diff <- .spd_logm(synthetic[[k]]) - anchor_log[[anchor_idx[k]]]
    d[k] <- sqrt(sum(diff^2))
  }
  mean(d)
}


#' Six-Metric Covariance Fidelity Audit
#'
#' @description
#' One-call wrapper that returns the six fidelity metrics used in the
#' multi-dataset audit of Shen \& Degras (2026): log-Euclidean distance
#' to a reference, affine-invariant Riemannian distance to a reference,
#' eigenvalue-spectrum correlation, trace ratio, condition-number ratio,
#' and anchor-to-synthetic perturbation distance. Returns the row both
#' in raw form and dimension-normalised by \eqn{\sqrt{p(p+1)/2}} so that
#' comparisons across datasets with different channel counts are not
#' dominated by the trivial dimension scaling.
#'
#' @param anchors List of SPD anchor covariance matrices.
#' @param synthetic List of SPD synthetic covariance matrices.
#' @param anchor_idx Integer vector mapping each synthetic to its anchor.
#' @param reference Optional SPD matrix used as the per-class reference for
#'   the LE / AI / spectrum / trace / condition metrics. Defaults to the
#'   element-wise mean of \code{anchors}.
#' @return A named list of metric values (raw and normalised).
#' @references
#' Shen, Y., Degras, D. (2026). Covariance Geometry as a Safety
#' Constraint for Cross-Session BCI Augmentation: A Multi-Dataset
#' Non-Inferiority and Fidelity Audit. \emph{(this manuscript).}
#' @seealso \code{\link{cov_logeuclidean_distance}},
#'   \code{\link{cov_affine_invariant_distance}},
#'   \code{\link{cov_eigenvalue_correlation}},
#'   \code{\link{cov_trace_ratio}},
#'   \code{\link{cov_condition_ratio}},
#'   \code{\link{cov_anchor_perturbation_distance}}.
#' @export
audit_covariance_fidelity <- function(anchors, synthetic, anchor_idx,
                                      reference = NULL) {
  if (length(synthetic) == 0L) {
    return(list(n_synthetic = 0L,
                log_euclidean_to_reference = NA_real_,
                affine_invariant_to_reference = NA_real_,
                eigenvalue_correlation = NA_real_,
                trace_ratio = NA_real_,
                condition_number_ratio = NA_real_,
                anchor_perturbation_distance = NA_real_,
                normalized = list()))
  }
  if (is.null(reference)) {
    reference <- 0.5 * (Reduce("+", anchors) / length(anchors) +
                          t(Reduce("+", anchors) / length(anchors)))
  }
  syn_mean <- Reduce("+", synthetic) / length(synthetic)
  syn_mean <- 0.5 * (syn_mean + t(syn_mean))

  le_to_ref <- mean(vapply(synthetic,
                           function(C) cov_logeuclidean_distance(C, reference),
                           numeric(1)))
  ai_to_ref <- mean(vapply(synthetic,
                           function(C) cov_affine_invariant_distance(C, reference),
                           numeric(1)))
  eig_corr  <- cov_eigenvalue_correlation(syn_mean, reference)
  tr_ratio  <- cov_trace_ratio(syn_mean, reference)
  cond_rat  <- cov_condition_ratio(syn_mean, reference)
  anc_dist  <- cov_anchor_perturbation_distance(anchors, synthetic,
                                                anchor_idx)

  p <- nrow(reference)
  d_norm <- sqrt(p * (p + 1) / 2)

  list(
    n_synthetic = length(synthetic),
    log_euclidean_to_reference = le_to_ref,
    affine_invariant_to_reference = ai_to_ref,
    eigenvalue_correlation = eig_corr,
    trace_ratio = tr_ratio,
    condition_number_ratio = cond_rat,
    anchor_perturbation_distance = anc_dist,
    normalized = list(
      log_euclidean = le_to_ref / d_norm,
      affine_invariant = ai_to_ref / d_norm,
      anchor_perturbation = anc_dist / d_norm,
      feature_dim = p * (p + 1) / 2
    )
  )
}
