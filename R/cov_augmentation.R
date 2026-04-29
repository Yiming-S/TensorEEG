#' Convert an EEG Tensor to a List of Trial Covariance Matrices
#'
#' @description
#' Computes one regularised symmetric positive-definite (SPD) covariance
#' matrix per trial from a 3rd-order EEG tensor. The output list is the
#' canonical input to the manifold-aware augmentation routines in this
#' package and to downstream Riemannian classifiers
#' (e.g.\ \code{pyriemann}'s \code{MDM}).
#'
#' @details
#' Given a tensor
#' \eqn{\mathcal{X} \in \mathbb{R}^{T \times C \times K}}, the trial-wise
#' covariance for trial \eqn{k} is
#' \deqn{\mathbf{C}_k = \frac{1}{T-1} \mathbf{X}_k^{\top} \mathbf{X}_k
#'                       + \lambda \mathbf{I}_C, \quad
#'        \mathbf{X}_k \in \mathbb{R}^{T \times C}, \quad
#'        \mathbf{C}_k \in \mathcal{S}_{++}^{C}.}
#' The ridge term \eqn{\lambda\mathbf{I}_C} guarantees strict positive
#' definiteness for short or low-rank trials. With
#' \code{centre = TRUE}, each trial is mean-centred per channel before
#' the cross-product, matching the convention of \code{stats::cov} and
#' \code{pyriemann.utils.covariance}.
#'
#' @param X_tensor Numeric array with dimensions
#'   \code{[n_time, n_channels, n_trials]} (the format returned by
#'   \code{\link{sim_eeg_master}}).
#' @param ridge Numeric. Diagonal regularisation \eqn{\lambda} added to
#'   each covariance to keep it strictly SPD.
#' @param centre Logical. If \code{TRUE}, mean-centre each trial per
#'   channel before computing the cross-product.
#'
#' @return A list of length \code{n_trials}. Each element is a
#'   \code{n_channels x n_channels} SPD matrix.
#'
#' @seealso \code{\link{augment_cov_riemannian}},
#'   \code{\link{cov_logeuclidean_distance}}.
#' @export
tensor_to_cov <- function(X_tensor, ridge = 1e-6, centre = TRUE) {
  if (!is.array(X_tensor) || length(dim(X_tensor)) != 3L) {
    stop("X_tensor must be a 3D array [time x channel x trial].")
  }
  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) ||
      ridge < 0) {
    stop("ridge must be a single non-negative finite number.")
  }
  if (!is.logical(centre) || length(centre) != 1L || is.na(centre)) {
    stop("centre must be TRUE or FALSE.")
  }

  dims <- dim(X_tensor)
  n_time <- dims[1]
  n_channels <- dims[2]
  n_trials <- dims[3]
  if (n_time < 2L) {
    stop("Each trial must have at least 2 time samples to estimate covariance.")
  }

  cov_list <- vector("list", n_trials)
  I_C <- diag(n_channels)
  for (k in seq_len(n_trials)) {
    Xk <- X_tensor[, , k, drop = FALSE]
    dim(Xk) <- c(n_time, n_channels)
    if (centre) {
      Xk <- sweep(Xk, 2, colMeans(Xk), "-")
    }
    C <- crossprod(Xk) / (n_time - 1L)
    C <- 0.5 * (C + t(C))
    cov_list[[k]] <- C + ridge * I_C
  }
  cov_list
}


# Internal: symmetric matrix logarithm via eigendecomposition. Faster and
# more numerically stable on SPD inputs than expm::logm, which uses Schur
# and may complain about complex eigenvalues from rounding.
.spd_logm <- function(C) {
  C <- 0.5 * (C + t(C))
  eig <- eigen(C, symmetric = TRUE)
  vals <- eig$values
  if (any(vals <= 0)) {
    stop("Covariance is not strictly positive definite (min eigenvalue ",
         signif(min(vals), 3), "). Increase 'ridge' in tensor_to_cov().")
  }
  eig$vectors %*% diag(log(vals), nrow = length(vals)) %*% t(eig$vectors)
}


# Internal: symmetric matrix exponential via eigendecomposition. Returns
# an SPD matrix when the input is symmetric.
.spd_expm <- function(Z) {
  Z <- 0.5 * (Z + t(Z))
  eig <- eigen(Z, symmetric = TRUE)
  eig$vectors %*% diag(exp(eig$values), nrow = length(eig$values)) %*% t(eig$vectors)
}


#' Log-Euclidean Distance Between SPD Matrices
#'
#' @description
#' Computes the log-Euclidean distance
#' \deqn{d_{\mathrm{LE}}(\mathbf{C}_1, \mathbf{C}_2)
#'   = \lVert \log(\mathbf{C}_1) - \log(\mathbf{C}_2) \rVert_F,}
#' a metric on the SPD manifold equivalent (up to constants) to the
#' geodesic distance under the affine-invariant metric in a neighbourhood
#' of the identity. It is the standard fidelity score for evaluating
#' covariance-aware augmentation against a real reference.
#'
#' @param C1,C2 Numeric SPD matrices of the same shape.
#'
#' @return A non-negative scalar.
#'
#' @seealso \code{\link{tensor_to_cov}},
#'   \code{\link{augment_cov_riemannian}}.
#' @references
#' Arsigny, V., Fillard, P., Pennec, X., Ayache, N. (2007). Geometric
#' means in a novel vector space structure on symmetric positive-definite
#' matrices. SIAM J. Matrix Anal. Appl., 29(1), 328--347.
#' @export
cov_logeuclidean_distance <- function(C1, C2) {
  if (!is.matrix(C1) || !is.matrix(C2) || !all(dim(C1) == dim(C2))) {
    stop("C1 and C2 must be matrices of the same shape.")
  }
  Z1 <- .spd_logm(C1)
  Z2 <- .spd_logm(C2)
  sqrt(sum((Z1 - Z2)^2))
}


#' Riemannian (Log-Euclidean) Augmentation of Trial Covariance Matrices
#'
#' @description
#' Implements the covariance-aware augmentation procedure used to
#' validate the EEG simulation framework on public BCI benchmarks. For
#' each anchor covariance matrix \eqn{\mathbf{C}_i} the function maps to
#' the log-Euclidean tangent space, adds class-preserving symmetric
#' Gaussian perturbation and an optional session-drift component, and
#' maps back to \eqn{\mathcal{S}_{++}^{C}}. Output matrices are SPD by
#' construction.
#'
#' @details
#' For each anchor \eqn{\mathbf{C}_i \in \mathcal{S}_{++}^{C}} and for
#' each of \code{n_aug} replicates the procedure is:
#' \enumerate{
#'   \item \eqn{\mathbf{Z}_i \gets \log(\mathbf{C}_i)} (symmetric matrix
#'         logarithm).
#'   \item Sample a symmetric perturbation
#'         \eqn{\mathbf{E}_i = (\mathbf{G}_i + \mathbf{G}_i^{\top})/2}
#'         with \eqn{\mathbf{G}_i \sim \mathcal{N}(0, \sigma^2 \mathbf{I})}.
#'   \item Optionally add a session-drift term
#'         \eqn{\mathbf{D}_{\Delta}(\theta)} (also symmetric) when
#'         \code{drift} is provided.
#'   \item \eqn{\tilde{\mathbf{Z}}_i \gets \mathbf{Z}_i + \mathbf{E}_i
#'                                + \mathbf{D}_{\Delta}(\theta).}
#'   \item \eqn{\tilde{\mathbf{C}}_i \gets \exp(\tilde{\mathbf{Z}}_i).}
#' }
#' The returned object is a flat list of \code{length(cov_list) * n_aug}
#' SPD matrices in anchor-major order; companion vectors carry the source
#' anchor index, the augmentation replicate index, and the propagated
#' class label.
#'
#' @param cov_list List of SPD matrices (e.g.\ output of
#'   \code{\link{tensor_to_cov}}). All matrices must share the same
#'   dimension.
#' @param n_aug Integer. Number of synthetic covariance replicates per
#'   anchor (default 5).
#' @param sigma Numeric. Standard deviation of the symmetric Gaussian
#'   tangent-space perturbation (default 0.10).
#' @param drift Optional numeric matrix of the same dimension as the
#'   covariance matrices. If supplied, added to every replicate (after
#'   symmetrisation) to simulate cross-session drift in tangent space.
#' @param labels Optional vector of length \code{length(cov_list)}. Each
#'   replicate inherits its anchor's label.
#' @param seed Optional integer. Sets the RNG seed for reproducibility.
#'
#' @return A list with elements
#' \describe{
#'   \item{\code{cov}}{List of length \code{length(cov_list) * n_aug} of
#'     SPD covariance matrices.}
#'   \item{\code{anchor}}{Integer vector of anchor indices.}
#'   \item{\code{replicate}}{Integer vector of within-anchor replicate
#'     indices.}
#'   \item{\code{labels}}{Vector of inherited labels (or \code{NULL} if
#'     \code{labels} was not supplied).}
#'   \item{\code{params}}{List of \code{sigma}, \code{n_aug}, and
#'     \code{seed} used.}
#' }
#'
#' @seealso \code{\link{tensor_to_cov}},
#'   \code{\link{cov_logeuclidean_distance}}.
#' @references
#' Pennec, X., Fillard, P., Ayache, N. (2006). A Riemannian framework
#' for tensor computing. Int. J. Comput. Vis., 66(1), 41--66.
#' Arsigny, V., Fillard, P., Pennec, X., Ayache, N. (2007). Geometric
#' means in a novel vector space structure on symmetric positive-definite
#' matrices. SIAM J. Matrix Anal. Appl., 29(1), 328--347.
#' @examples
#' set.seed(1)
#' sim <- sim_eeg_master(n_trials = 8, n_time = 200, n_channels = 6,
#'                       n_sources = 4, seed = 1, verbose = FALSE)
#' cov_list <- tensor_to_cov(sim$data)
#' aug <- augment_cov_riemannian(cov_list, n_aug = 3, sigma = 0.05,
#'                               labels = sim$labels, seed = 42)
#' length(aug$cov)
#' table(aug$anchor)
#' @importFrom stats rnorm
#' @export
augment_cov_riemannian <- function(cov_list, n_aug = 5L, sigma = 0.10,
                                   drift = NULL, labels = NULL,
                                   seed = NULL) {
  if (!is.list(cov_list) || length(cov_list) < 1L) {
    stop("cov_list must be a non-empty list of SPD matrices.")
  }
  if (!.is_whole_number(n_aug) || n_aug < 1L) {
    stop("n_aug must be a positive integer.")
  }
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) ||
      sigma < 0) {
    stop("sigma must be a single non-negative finite number.")
  }
  n_aug <- as.integer(round(n_aug))
  .warn_synthetic_real_ratio(n_aug, length(cov_list))

  p <- nrow(cov_list[[1L]])
  if (is.null(p) || p < 1L) {
    stop("First element of cov_list is not a square matrix.")
  }
  for (k in seq_along(cov_list)) {
    Ck <- cov_list[[k]]
    if (!is.matrix(Ck) || nrow(Ck) != p || ncol(Ck) != p) {
      stop(sprintf("cov_list[[%d]] must be a %d x %d matrix.", k, p, p))
    }
  }

  drift_sym <- NULL
  if (!is.null(drift)) {
    if (!is.matrix(drift) || nrow(drift) != p || ncol(drift) != p) {
      stop(sprintf("drift must be a %d x %d matrix.", p, p))
    }
    drift_sym <- 0.5 * (drift + t(drift))
  }

  if (!is.null(labels)) {
    if (length(labels) != length(cov_list)) {
      stop("labels must have length equal to length(cov_list).")
    }
  }

  if (!is.null(seed)) {
    if (!.is_whole_number(seed)) {
      stop("seed must be a single finite integer.")
    }
    set.seed(as.integer(round(seed)))
  }

  total <- length(cov_list) * n_aug
  out_cov <- vector("list", total)
  anchor_idx <- integer(total)
  rep_idx <- integer(total)
  out_labels <- if (is.null(labels)) NULL else rep(labels[1], total)

  pos <- 1L
  for (i in seq_along(cov_list)) {
    Z_i <- .spd_logm(cov_list[[i]])
    for (r in seq_len(n_aug)) {
      G <- matrix(stats::rnorm(p * p, sd = sigma), p, p)
      E <- 0.5 * (G + t(G))
      Z_aug <- Z_i + E
      if (!is.null(drift_sym)) {
        Z_aug <- Z_aug + drift_sym
      }
      out_cov[[pos]] <- .spd_expm(Z_aug)
      anchor_idx[pos] <- i
      rep_idx[pos] <- r
      if (!is.null(out_labels)) {
        out_labels[pos] <- labels[i]
      }
      pos <- pos + 1L
    }
  }

  list(
    cov = out_cov,
    anchor = anchor_idx,
    replicate = rep_idx,
    labels = out_labels,
    params = list(sigma = sigma, n_aug = n_aug, seed = seed)
  )
}


# Internal: project a symmetric matrix onto the SPD cone by clipping
# eigenvalues at a small floor.
.spd_project <- function(M, floor = 1e-6) {
  S <- 0.5 * (M + t(M))
  eig <- eigen(S, symmetric = TRUE)
  vals <- pmax(eig$values, floor)
  eig$vectors %*% diag(vals, nrow = length(vals)) %*% t(eig$vectors)
}


# Internal: anchor-to-augmented log-Euclidean distance vector for a paired
# (anchors, augmented, anchor_idx) tuple.
.anchor_log_distances <- function(anchors_log, augmented_cov, anchor_idx) {
  d <- numeric(length(augmented_cov))
  for (k in seq_along(augmented_cov)) {
    Z_aug <- .spd_logm(augmented_cov[[k]])
    diff <- Z_aug - anchors_log[[anchor_idx[k]]]
    d[k] <- sqrt(sum(diff^2))
  }
  d
}


# Internal: budget-ratio warning. When the synthetic-to-real ratio
# (n_aug per anchor) is large, isotropic tangent jitter can dominate
# the LDA pooled-covariance estimate at moderate calibration budgets;
# reproducible failure mode documented for G0 in Shen & Degras (2026).
# The warning is restricted to the regime where the failure was actually
# observed (10 <= n_anchor < 40 with r >= 3); toy-sized anchor lists
# (n_anchor < 10) are typically unit-test fixtures and are skipped.
.warn_synthetic_real_ratio <- function(n_aug, n_anchor,
                                       ratio_threshold = 3L,
                                       lower_anchor = 10L,
                                       upper_anchor = 40L) {
  if (n_aug >= ratio_threshold &&
      n_anchor >= lower_anchor &&
      n_anchor < upper_anchor) {
    warning(sprintf(
      paste0("Synthetic-to-real ratio %d:1 with only %d real anchors. ",
             "Isotropic tangent jitter (G0) showed a reproducible ",
             "budget-dependent failure on BNCI2014_001 at ncal=30 with ",
             "r=3 in Shen & Degras (2026); consider G2 (geodesic mixup) ",
             "or a smaller n_aug if downstream accuracy is the goal."),
      as.integer(n_aug), as.integer(n_anchor)),
      call. = FALSE)
  }
  invisible(NULL)
}


#' Amplitude-Matched Euclidean Covariance Perturbation (E0)
#'
#' @description
#' Off-manifold negative-control augmentation that perturbs covariance
#' matrices in the ambient Euclidean space and projects back to the SPD
#' cone. The Gaussian noise scale \eqn{\sigma_E} is calibrated so that
#' the median anchor-to-augmented log-Euclidean distance matches that
#' of \code{\link{augment_cov_riemannian}} (G0) at a reference \code{sigma}.
#' This removes the ``unfair magnitude'' critique that an uncalibrated
#' Euclidean baseline would otherwise invite. Synthetic covariances
#' inherit anchor labels with no synthetic label perturbation.
#'
#' @details
#' For each anchor \eqn{\mathbf{C}_i \in \mathcal{S}_{++}^{C}}:
#' \enumerate{
#'   \item Compute the target median log-Euclidean distance from the
#'         G0 generator at \code{g0_sigma}.
#'   \item Search a single scalar \eqn{\sigma_E} on a log-spaced grid
#'         (default 41 points spanning \eqn{[10^{-3}, 10^{2}]}) so
#'         that \eqn{\mathrm{median}(d_{\mathrm{LE}}(\mathbf{C}_i,
#'         \tilde{\mathbf{C}}_{ij}))} matches the G0 target.
#'   \item Final E0 stack: sample symmetric Gaussian noise
#'         \eqn{\mathbf{H}_{ik} = (\mathbf{N}_{ik} + \mathbf{N}_{ik}^{\top})/2},
#'         \eqn{\mathbf{N}_{ik}\sim\mathcal{N}(0,\sigma_E^2 \mathbf{I})},
#'         add to \eqn{\mathbf{C}_i}, and project back to the SPD cone
#'         by clipping eigenvalues at \code{eigenvalue_floor}.
#' }
#'
#' @param cov_list List of SPD anchor covariance matrices.
#' @param n_aug Integer. Number of synthetic replicates per anchor.
#' @param g0_sigma Numeric. Tangent perturbation scale that defines the
#'   matching target for E0 (default 0.10, paired with the default of
#'   \code{\link{augment_cov_riemannian}}).
#' @param labels Optional vector of anchor labels (inherited verbatim).
#' @param seed Optional integer RNG seed for reproducibility.
#' @param eigenvalue_floor Positive scalar floor used by the SPD
#'   projection.
#' @param tolerance_relative Relative tolerance for the amplitude match
#'   (target \eqn{|\rho - 1| \le} this; default 0.10).
#'
#' @return A list with the same shape as
#'   \code{\link{augment_cov_riemannian}} plus a \code{diagnostic}
#'   element giving the achieved amplitude-match \eqn{\rho =
#'   \mathrm{median}(d_{E0}) / \mathrm{median}(d_{G0})}.
#' @seealso \code{\link{augment_cov_riemannian}}.
#' @importFrom stats rnorm median
#' @export
augment_cov_amplitude_matched_euclidean <- function(cov_list, n_aug = 5L,
                                                     g0_sigma = 0.10,
                                                     labels = NULL,
                                                     seed = NULL,
                                                     eigenvalue_floor = 1e-6,
                                                     tolerance_relative = 0.10) {
  if (!is.list(cov_list) || length(cov_list) < 1L) {
    stop("cov_list must be a non-empty list of SPD matrices.")
  }
  if (!.is_whole_number(n_aug) || n_aug < 1L) {
    stop("n_aug must be a positive integer.")
  }
  if (!is.numeric(g0_sigma) || length(g0_sigma) != 1L || g0_sigma <= 0) {
    stop("g0_sigma must be a single positive number.")
  }
  if (!is.numeric(eigenvalue_floor) || eigenvalue_floor <= 0) {
    stop("eigenvalue_floor must be positive.")
  }
  n_aug <- as.integer(round(n_aug))

  p <- nrow(cov_list[[1L]])
  for (k in seq_along(cov_list)) {
    Ck <- cov_list[[k]]
    if (!is.matrix(Ck) || nrow(Ck) != p || ncol(Ck) != p) {
      stop(sprintf("cov_list[[%d]] must be a %d x %d matrix.", k, p, p))
    }
  }

  if (!is.null(seed)) {
    if (!.is_whole_number(seed)) stop("seed must be a single integer.")
    set.seed(as.integer(round(seed)))
  }

  anchors_log <- lapply(cov_list, .spd_logm)

  # Target distance from G0 (deterministic given the same seed call).
  total <- length(cov_list) * n_aug
  G0_log_distances <- numeric(total)
  pos <- 1L
  for (i in seq_along(cov_list)) {
    Z_i <- anchors_log[[i]]
    for (r in seq_len(n_aug)) {
      G <- matrix(stats::rnorm(p * p, sd = g0_sigma), p, p)
      E <- 0.5 * (G + t(G))
      G0_log_distances[pos] <- sqrt(sum(E^2))
      pos <- pos + 1L
    }
  }
  target_distance <- stats::median(G0_log_distances)

  # E0 sigma search: log-spaced grid + nearest-distance pick.
  sigma_grid <- 10^seq(-3, 2, length.out = 41)
  search_anchors_log <- anchors_log

  distance_at_sigma <- function(sigma_e) {
    out <- vector("list", total)
    idx <- integer(total)
    pos <- 1L
    for (i in seq_along(cov_list)) {
      Ci <- cov_list[[i]]
      for (r in seq_len(n_aug)) {
        N <- matrix(stats::rnorm(p * p, sd = sigma_e), p, p)
        H <- 0.5 * (N + t(N))
        out[[pos]] <- .spd_project(Ci + H, floor = eigenvalue_floor)
        idx[pos] <- i
        pos <- pos + 1L
      }
    }
    list(cov = out, idx = idx,
         med = stats::median(.anchor_log_distances(search_anchors_log,
                                                    out, idx)))
  }

  achieved <- vapply(sigma_grid, function(s) distance_at_sigma(s)$med,
                     numeric(1))
  best_idx <- which.min(abs(achieved - target_distance))
  sigma_e <- sigma_grid[best_idx]

  result <- distance_at_sigma(sigma_e)
  rho <- result$med / target_distance

  out_labels <- if (is.null(labels)) NULL else labels[result$idx]
  list(
    cov = result$cov,
    anchor = result$idx,
    replicate = rep(seq_len(n_aug), times = length(cov_list)),
    labels = out_labels,
    params = list(g0_sigma = g0_sigma, sigma_e = sigma_e, n_aug = n_aug,
                  seed = seed, eigenvalue_floor = eigenvalue_floor),
    diagnostic = list(target_distance_g0 = target_distance,
                      achieved_distance_e0 = result$med,
                      rho = rho,
                      tolerance = tolerance_relative,
                      success = abs(rho - 1) <= tolerance_relative)
  )
}


#' Empirical Tangent Gaussian Augmentation (G1)
#'
#' @description
#' Geometry-preserving augmentation that adds class-aware Gaussian
#' perturbation in the log-Euclidean tangent space. For each class
#' \eqn{c}, the empirical tangent covariance \eqn{\hat\Sigma_c} is
#' estimated from the half-vectorised log-covariances of the same-class
#' anchors via Ledoit--Wolf-style shrinkage; samples are drawn from
#' \eqn{\mathcal{N}(\mathbf{0}, \hat\Sigma_c)} and added at scale
#' \code{sigma}.
#'
#' @param cov_list List of SPD anchor covariance matrices.
#' @param labels Vector of anchor labels (required for class-aware
#'   dispersion estimation).
#' @param n_aug Integer. Number of synthetic replicates per anchor.
#' @param sigma Numeric tangent perturbation scale (default 0.10).
#' @param shrinkage Numeric in [0, 1]. Linear shrinkage intensity toward
#'   the diagonal (default 0.10). Set to 0 to use the raw sample
#'   covariance when class size permits.
#' @param seed Optional integer RNG seed.
#'
#' @return A list with elements \code{cov}, \code{anchor}, \code{replicate},
#'   \code{labels}, \code{params}.
#' @importFrom stats rnorm cov
#' @export
augment_cov_empirical_tangent <- function(cov_list, labels, n_aug = 5L,
                                          sigma = 0.10, shrinkage = 0.10,
                                          seed = NULL) {
  if (!is.list(cov_list) || length(cov_list) < 2L) {
    stop("cov_list must contain at least two anchors.")
  }
  if (length(labels) != length(cov_list)) {
    stop("labels must have length equal to length(cov_list).")
  }
  if (!.is_whole_number(n_aug) || n_aug < 1L) stop("n_aug must be positive integer.")
  if (sigma < 0) stop("sigma must be non-negative.")
  if (shrinkage < 0 || shrinkage > 1) stop("shrinkage must lie in [0, 1].")
  n_aug <- as.integer(round(n_aug))

  p <- nrow(cov_list[[1L]])
  d_dim <- p * (p + 1L) / 2L
  triu <- which(upper.tri(matrix(0, p, p), diag = TRUE), arr.ind = TRUE)
  weight <- ifelse(triu[, 1] == triu[, 2], 1, sqrt(2))

  if (!is.null(seed)) {
    if (!.is_whole_number(seed)) stop("seed must be a single integer.")
    set.seed(as.integer(round(seed)))
  }

  anchors_log <- lapply(cov_list, .spd_logm)
  vech_log <- function(Z) Z[triu] * weight
  unvech <- function(z) {
    M <- matrix(0, p, p)
    M[triu] <- z / weight
    M[lower.tri(M)] <- t(M)[lower.tri(M)]
    0.5 * (M + t(M))
  }

  z_anchors <- t(vapply(anchors_log, vech_log, numeric(d_dim)))
  classes <- sort(unique(labels))
  cov_per_class <- vector("list", length(classes))
  names(cov_per_class) <- as.character(classes)
  for (cls in classes) {
    mask <- labels == cls
    z_cls <- z_anchors[mask, , drop = FALSE]
    if (nrow(z_cls) < 2L) {
      cov_per_class[[as.character(cls)]] <- diag(d_dim)
      next
    }
    sample_cov <- stats::cov(z_cls)
    target <- diag(diag(sample_cov))   # diagonal target for shrinkage
    cov_per_class[[as.character(cls)]] <-
      (1 - shrinkage) * sample_cov + shrinkage * target +
      1e-10 * diag(d_dim)
  }

  total <- length(cov_list) * n_aug
  out_cov <- vector("list", total)
  anchor_idx <- integer(total)
  rep_idx <- integer(total)
  out_labels <- labels[1]
  out_labels <- rep(out_labels, total)

  pos <- 1L
  for (i in seq_along(cov_list)) {
    cls <- labels[i]
    Sigma_c <- cov_per_class[[as.character(cls)]]
    chol_c <- tryCatch(chol(Sigma_c), error = function(e) NULL)
    z_i <- z_anchors[i, ]
    Z_i_log <- anchors_log[[i]]
    for (r in seq_len(n_aug)) {
      eps <- if (!is.null(chol_c)) {
        as.numeric(stats::rnorm(d_dim) %*% chol_c)
      } else {
        ev <- eigen(Sigma_c, symmetric = TRUE)
        as.numeric(stats::rnorm(d_dim) * sqrt(pmax(ev$values, 0)) %*% t(ev$vectors))
      }
      z_aug <- z_i + sigma * eps
      Z_aug <- unvech(z_aug)
      out_cov[[pos]] <- .spd_expm(Z_aug)
      anchor_idx[pos] <- i
      rep_idx[pos] <- r
      out_labels[pos] <- cls
      pos <- pos + 1L
    }
  }

  list(cov = out_cov, anchor = anchor_idx, replicate = rep_idx,
       labels = out_labels,
       params = list(sigma = sigma, shrinkage = shrinkage,
                     n_aug = n_aug, seed = seed))
}


#' Log-Euclidean Geodesic Mixup Augmentation (G2)
#'
#' @description
#' Class-aware log-Euclidean geodesic interpolation between same-class
#' anchor covariances. For each anchor and replicate, a same-class
#' partner \eqn{\mathbf{C}_j} and a Beta-distributed weight
#' \eqn{\alpha \sim \mathrm{Beta}(a, a)} are sampled, and the synthetic
#' covariance is
#' \deqn{\tilde{\mathbf{C}} = \exp\!\big( (1-\alpha) \log \mathbf{C}_i
#'   + \alpha \log \mathbf{C}_j \big).}
#' G2 was the most stable geometry-preserving variant in the audit of
#' Shen \& Degras (2026), with pooled paired-difference standard deviation
#' \eqn{\hat\tau = 0.08} pp across four datasets at \eqn{n_\mathrm{cal}=20}.
#'
#' @param cov_list List of SPD anchor covariance matrices.
#' @param labels Vector of anchor labels.
#' @param n_aug Integer. Number of synthetic replicates per anchor.
#' @param beta_alpha Positive numeric. Symmetric Beta concentration
#'   (default 1.0 = uniform on \eqn{[0, 1]}).
#' @param seed Optional integer RNG seed.
#'
#' @return A list with \code{cov}, \code{anchor}, \code{replicate},
#'   \code{labels}, \code{params}, plus a \code{trace} sublist recording
#'   the partner index and \eqn{\alpha} for each synthetic trial.
#' @importFrom stats rbeta
#' @export
augment_cov_geodesic_mixup <- function(cov_list, labels, n_aug = 5L,
                                       beta_alpha = 1.0, seed = NULL) {
  if (!is.list(cov_list) || length(cov_list) < 2L) {
    stop("cov_list must contain at least two anchors.")
  }
  if (length(labels) != length(cov_list)) {
    stop("labels must have length equal to length(cov_list).")
  }
  if (!.is_whole_number(n_aug) || n_aug < 1L) stop("n_aug must be positive integer.")
  if (beta_alpha <= 0) stop("beta_alpha must be positive.")
  n_aug <- as.integer(round(n_aug))

  p <- nrow(cov_list[[1L]])
  if (!is.null(seed)) {
    if (!.is_whole_number(seed)) stop("seed must be a single integer.")
    set.seed(as.integer(round(seed)))
  }
  anchors_log <- lapply(cov_list, .spd_logm)

  classes <- sort(unique(labels))
  same_class <- list()
  for (cls in classes) same_class[[as.character(cls)]] <- which(labels == cls)

  total <- length(cov_list) * n_aug
  out_cov <- vector("list", total)
  anchor_idx <- integer(total)
  rep_idx <- integer(total)
  partner_idx <- integer(total)
  alpha_vals <- numeric(total)
  out_labels <- labels[1]
  out_labels <- rep(out_labels, total)

  pos <- 1L
  for (i in seq_along(cov_list)) {
    cls <- labels[i]
    candidates <- setdiff(same_class[[as.character(cls)]], i)
    if (length(candidates) == 0L) candidates <- same_class[[as.character(cls)]]
    Z_i <- anchors_log[[i]]
    for (r in seq_len(n_aug)) {
      partner <- if (length(candidates) == 1L) candidates else
        sample(candidates, size = 1L)
      alpha <- stats::rbeta(1, beta_alpha, beta_alpha)
      Z_j <- anchors_log[[partner]]
      Z_aug <- (1 - alpha) * Z_i + alpha * Z_j
      out_cov[[pos]] <- .spd_expm(Z_aug)
      anchor_idx[pos] <- i
      rep_idx[pos] <- r
      partner_idx[pos] <- partner
      alpha_vals[pos] <- alpha
      out_labels[pos] <- cls
      pos <- pos + 1L
    }
  }

  list(cov = out_cov, anchor = anchor_idx, replicate = rep_idx,
       labels = out_labels,
       trace = list(partner = partner_idx, alpha = alpha_vals),
       params = list(beta_alpha = beta_alpha, n_aug = n_aug, seed = seed))
}


# Internal: affine-invariant Riemannian (geometric) mean of an SPD stack.
.riemann_mean <- function(cov_list, max_iter = 30L, tol = 1e-6) {
  p <- nrow(cov_list[[1L]])
  arith <- Reduce("+", cov_list) / length(cov_list)
  M <- 0.5 * (arith + t(arith))
  for (it in seq_len(max_iter)) {
    eig <- eigen(M, symmetric = TRUE)
    vals <- pmax(eig$values, 1e-12)
    inv_sqrt <- eig$vectors %*% diag(vals^(-0.5), p) %*% t(eig$vectors)
    sqrt_M  <- eig$vectors %*% diag(vals^( 0.5), p) %*% t(eig$vectors)
    tangent_sum <- matrix(0, p, p)
    for (k in seq_along(cov_list)) {
      W <- inv_sqrt %*% cov_list[[k]] %*% inv_sqrt
      tangent_sum <- tangent_sum + .spd_logm(W)
    }
    tangent_mean <- tangent_sum / length(cov_list)
    if (sqrt(sum(tangent_mean^2)) < tol) break
    M <- 0.5 * (sqrt_M %*% .spd_expm(tangent_mean) %*% sqrt_M +
                  t(sqrt_M %*% .spd_expm(tangent_mean) %*% sqrt_M))
  }
  M
}


#' Riemannian Alignment of Source Covariances Toward a Target Manifold (A0)
#'
#' @description
#' Transductive alignment baseline. Whitens the source-session anchor
#' covariances toward the affine-invariant Riemannian mean of an
#' \emph{unlabeled} target-session covariance stack via
#' \deqn{\mathbf{W} = \mathbf{M}_T^{1/2} \mathbf{M}_S^{-1/2},
#'   \quad \tilde{\mathbf{C}}_i = \mathbf{W} \mathbf{C}_i \mathbf{W}^{\top}.}
#' A0 consults unlabeled target covariances and is therefore a
#' \emph{transductive} method; it does not enter the non-transductive
#' primary non-inferiority family of Shen \& Degras (2026) and is
#' reported there as a transfer-learning upper bound.
#'
#' @param source_cov_list List of SPD source-session covariances (anchors).
#' @param target_cov_list List of SPD target-session covariances
#'   (\strong{unlabeled} target data; only their geometry is used).
#' @param labels Optional vector of anchor labels (passed through).
#'
#' @return A list with elements \code{cov} (aligned anchors),
#'   \code{labels}, and \code{params}.
#' @export
augment_cov_alignment_riemannian <- function(source_cov_list,
                                             target_cov_list,
                                             labels = NULL) {
  if (!is.list(source_cov_list) || length(source_cov_list) < 1L) {
    stop("source_cov_list must be a non-empty list of SPD matrices.")
  }
  if (!is.list(target_cov_list) || length(target_cov_list) < 1L) {
    stop("target_cov_list must be a non-empty list of SPD matrices.")
  }
  p <- nrow(source_cov_list[[1L]])
  if (nrow(target_cov_list[[1L]]) != p) {
    stop("source and target covariances must have the same dimension.")
  }

  M_S <- .riemann_mean(source_cov_list)
  M_T <- .riemann_mean(target_cov_list)

  eig_s <- eigen(M_S, symmetric = TRUE)
  vals_s <- pmax(eig_s$values, 1e-12)
  src_inv_sqrt <- eig_s$vectors %*% diag(vals_s^(-0.5), p) %*% t(eig_s$vectors)

  eig_t <- eigen(M_T, symmetric = TRUE)
  vals_t <- pmax(eig_t$values, 1e-12)
  tgt_sqrt <- eig_t$vectors %*% diag(vals_t^(0.5), p) %*% t(eig_t$vectors)

  W <- tgt_sqrt %*% src_inv_sqrt

  aligned <- lapply(source_cov_list, function(Ci) {
    A <- W %*% Ci %*% t(W)
    0.5 * (A + t(A))
  })

  list(cov = aligned, labels = labels,
       params = list(uses_target = "unlabeled_covariance_only"))
}
