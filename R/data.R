#' Example SPD Anchor Covariances for Demonstrations
#'
#' @description
#' A small two-class stack of symmetric positive-definite (SPD) covariance
#' matrices used by the package's vignettes, demos, and integration tests.
#' Each anchor is a 6 x 6 covariance computed from a deterministic
#' synthetic 6-channel trial: class 0 carries a 10 Hz alpha-like
#' oscillation and class 1 carries a 6 Hz lower-frequency oscillation,
#' with a small per-channel frequency jitter and white-noise floor so
#' that within-class variability is non-trivial. A 1e-3 ridge is added on
#' the diagonal to guarantee strict positive definiteness across all
#' anchors.
#'
#' This data set is intentionally small so that
#' \code{\link{augment_cov_riemannian}},
#' \code{\link{augment_cov_amplitude_matched_euclidean}},
#' \code{\link{augment_cov_empirical_tangent}},
#' \code{\link{augment_cov_geodesic_mixup}}, and
#' \code{\link{audit_covariance_fidelity}} can all run in seconds.
#'
#' @format A list of length 8. Each element is a 6 x 6 numeric SPD
#'   matrix. Anchors 1--4 belong to class 0 (10 Hz); anchors 5--8 belong
#'   to class 1 (6 Hz). The companion label vector is shipped at
#'   \code{system.file("extdata", "example_labels.rds", package = "TensorEEG")}.
#'
#' @source Generated deterministically by
#'   \code{data-raw/generate_examples.R} with seed 2026 and the
#'   per-anchor seeds 101--104 (class 0) and 201--204 (class 1).
#'
#' @examples
#' data(example_anchors)
#' length(example_anchors)
#' dim(example_anchors[[1]])
#' labels <- readRDS(system.file("extdata", "example_labels.rds",
#'                               package = "TensorEEG"))
#' table(labels)
"example_anchors"
