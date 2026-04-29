#' Read a Calibration Manifest CSV
#'
#' @description
#' Reads the per-row calibration manifest produced by the Python protocol
#' driver
#' (\code{scripts/protocol/manifest.py::write_manifest}; one row per
#' \code{(dataset, subject, budget, resample, method, class)}). The
#' manifest is the audit trail that lets a third party replay synthetic
#' covariances from the recorded seeds; this function loads it into R for
#' use with \code{\link{replay_from_manifest}}.
#'
#' @param path File path to a \code{calibration_manifest.csv}.
#' @return A \code{data.frame} with columns matching
#'   \code{scripts/protocol/manifest.py::CSV_COLUMNS}, including
#'   \code{real_trial_ids} parsed as a list-column of integer vectors.
#' @importFrom utils read.csv
#' @export
read_calibration_manifest <- function(path) {
  if (!file.exists(path)) stop("manifest file not found: ", path)
  df <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("dataset", "subject", "budget_per_class", "resample_id",
                "method_code", "class_label", "n_real_trials",
                "real_trial_ids", "real_seed", "method_seed",
                "n_synthetic_per_anchor", "target_session_ids_used",
                "status", "notes")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0L) {
    stop("manifest is missing required columns: ",
         paste(missing, collapse = ", "))
  }
  df$real_trial_ids <- lapply(df$real_trial_ids, .parse_int_list_json)
  df
}


# Internal: parse a JSON-like integer list "[1, 2, 3]" into a numeric
# vector. We avoid jsonlite to keep TensorEEG's Imports minimal.
.parse_int_list_json <- function(s) {
  s <- gsub("[\\[\\]\\s]", "", s, perl = TRUE)
  if (!nzchar(s)) return(integer(0))
  as.integer(strsplit(s, ",", fixed = TRUE)[[1L]])
}


#' Replay Synthetic Covariance Stack From a Manifest Row Group
#'
#' @description
#' Reconstructs the synthetic covariance stack for one
#' \code{(dataset, subject, budget, resample, method)} cell of a
#' calibration manifest, using the per-method seed recorded in the
#' manifest and the source-session covariance list provided by the
#' caller. This is the R-side analogue of
#' \code{scripts/protocol/run_fidelity.py::_regenerate_synthetic} and
#' makes the protocol's reproducibility infrastructure available to R
#' users without leaving the language.
#'
#' @details
#' The manifest is grouped by \code{(dataset, subject, budget,
#' resample, method)} into class-specific rows; this function takes one
#' such group, reads the recorded \code{real_trial_ids} per class to
#' assemble anchors from \code{source_cov_list}, and dispatches to the
#' appropriate augmentation routine
#' (\code{\link{augment_cov_riemannian}} for G0,
#' \code{\link{augment_cov_amplitude_matched_euclidean}} for E0,
#' \code{\link{augment_cov_empirical_tangent}} for G1, or
#' \code{\link{augment_cov_geodesic_mixup}} for G2). R0 returns an
#' empty stack and the anchors verbatim. A0 is intentionally not
#' supported here because it is transductive and requires unlabeled
#' target covariances rather than seeds.
#'
#' @param manifest_group A \code{data.frame} subset corresponding to a
#'   single \code{(dataset, subject, budget, resample, method)} cell.
#'   Typically obtained as a row group from
#'   \code{\link{read_calibration_manifest}}.
#' @param source_cov_list List of SPD source-session covariance
#'   matrices, indexed by trial id. The function selects elements by
#'   \code{real_trial_ids}.
#' @param source_labels Vector of source-session labels of length
#'   \code{length(source_cov_list)}.
#' @param ratio_per_anchor Augmentation ratio (\code{n_aug}) used in
#'   the original run. Default 3 to match the primary protocol.
#' @param g0_sigma Tangent perturbation scale used by G0 (and as the
#'   amplitude-matching reference for E0). Default 0.15.
#' @param g1_sigma Empirical-tangent scale for G1. Default 0.15.
#' @param g2_beta_alpha Symmetric Beta concentration for G2. Default 1.
#'
#' @return A list with elements \code{cov} (synthetic stack),
#'   \code{labels}, \code{anchor}, \code{anchors_used} (the source
#'   covariances actually consumed by this cell), and \code{params}.
#' @seealso \code{\link{read_calibration_manifest}},
#'   \code{\link{audit_covariance_fidelity}}.
#' @export
replay_from_manifest <- function(manifest_group, source_cov_list,
                                 source_labels,
                                 ratio_per_anchor = 3L,
                                 g0_sigma = 0.15,
                                 g1_sigma = 0.15,
                                 g2_beta_alpha = 1.0) {
  if (!is.data.frame(manifest_group) || nrow(manifest_group) < 1L) {
    stop("manifest_group must be a non-empty data.frame.")
  }
  cells <- unique(manifest_group[, c("dataset", "subject",
                                     "budget_per_class", "resample_id",
                                     "method_code")])
  if (nrow(cells) != 1L) {
    stop("manifest_group must contain exactly one ",
         "(dataset, subject, budget, resample, method) cell.")
  }
  method <- cells$method_code[1]
  method_seed <- as.integer(manifest_group$method_seed[1])
  if (length(source_cov_list) != length(source_labels)) {
    stop("source_cov_list and source_labels must have equal length.")
  }

  # Assemble anchors from real_trial_ids in sorted-class order.
  classes <- sort(unique(manifest_group$class_label))
  anchors <- list()
  anchor_labels <- integer(0)
  for (cls in classes) {
    rows <- manifest_group[manifest_group$class_label == cls, , drop = FALSE]
    if (nrow(rows) != 1L) {
      stop(sprintf("manifest cell has %d rows for class %d; expected 1.",
                   nrow(rows), cls))
    }
    ids <- rows$real_trial_ids[[1]]
    # Manifest indices are 0-based (Python). R is 1-based.
    r_ids <- as.integer(ids) + 1L
    if (any(r_ids < 1L) || any(r_ids > length(source_cov_list))) {
      stop("manifest real_trial_ids reference indices outside source range.")
    }
    anchors <- c(anchors, source_cov_list[r_ids])
    anchor_labels <- c(anchor_labels, rep(as.integer(cls), length(r_ids)))
  }

  if (method == "R0") {
    return(list(cov = list(), labels = integer(0), anchor = integer(0),
                anchors_used = anchors,
                params = list(method = "R0", method_seed = method_seed)))
  }
  if (method == "G0") {
    res <- augment_cov_riemannian(anchors, n_aug = ratio_per_anchor,
                                   sigma = g0_sigma,
                                   labels = anchor_labels,
                                   seed = method_seed)
    res$anchors_used <- anchors
    return(res)
  }
  if (method == "G1") {
    res <- augment_cov_empirical_tangent(anchors, anchor_labels,
                                          n_aug = ratio_per_anchor,
                                          sigma = g1_sigma,
                                          seed = method_seed)
    res$anchors_used <- anchors
    return(res)
  }
  if (method == "G2") {
    res <- augment_cov_geodesic_mixup(anchors, anchor_labels,
                                       n_aug = ratio_per_anchor,
                                       beta_alpha = g2_beta_alpha,
                                       seed = method_seed)
    res$anchors_used <- anchors
    return(res)
  }
  if (method == "E0") {
    res <- augment_cov_amplitude_matched_euclidean(
      anchors, n_aug = ratio_per_anchor,
      g0_sigma = g0_sigma,
      labels = anchor_labels,
      seed = method_seed)
    res$anchors_used <- anchors
    return(res)
  }
  if (method == "A0") {
    stop("A0 (Riemannian alignment) is transductive and cannot be ",
         "replayed from manifest seeds alone. Use ",
         "augment_cov_alignment_riemannian() with the unlabeled ",
         "target covariance stack instead.")
  }
  stop("unknown method_code in manifest_group: ", method)
}
