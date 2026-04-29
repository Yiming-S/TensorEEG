# ----------------------------------------------------------------------------
# data-raw/generate_examples.R
#
# Regenerates the example assets shipped with TensorEEG:
#   data/example_anchors.rda       -- lazy-loaded list of SPD anchor matrices
#   inst/extdata/example_anchors.rds   -- same anchors, file-based access
#   inst/extdata/example_labels.rds    -- class labels for the anchors
#   inst/extdata/example_manifest.csv  -- minimal calibration manifest CSV
#
# Run from the package root with:
#   Rscript data-raw/generate_examples.R
#
# The output is deterministic given the seeds used here. The anchors are
# small enough (p = 6, N = 8) that vignettes and demos can run in a few
# seconds, but large enough to exercise the augmentation and fidelity
# routines on a realistic SPD-manifold geometry.
# ----------------------------------------------------------------------------

set.seed(2026L)

# Load the package source files directly so this script does not depend
# on a previously installed TensorEEG version.
pkg_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."),
                          mustWork = FALSE)
if (!dir.exists(file.path(pkg_root, "R"))) {
  pkg_root <- normalizePath(".")
}
src_files <- list.files(file.path(pkg_root, "R"), pattern = "[.]R$",
                        full.names = TRUE)
for (f in src_files) source(f)

# ---------------------------------------------------------------------------
# 1. SPD anchor stack: 8 anchors, p = 6 channels, two classes.
#    Class A is built from a 10 Hz alpha-like correlation pattern;
#    class B from a slower 6 Hz pattern. Each class has 4 noisy realisations
#    so that augmentation routines (G1, G2) have within-class variability
#    to learn from.
# ---------------------------------------------------------------------------
p <- 6L
n_per_class <- 4L

build_anchor <- function(base_freq, jitter_seed) {
  set.seed(jitter_seed)
  # 6-channel synthetic trial: cos waves at slightly perturbed frequencies
  # plus a small amount of white noise per channel.
  fs <- 250
  n_time <- 500
  t <- seq_len(n_time) / fs
  X <- matrix(0, n_time, p)
  for (k in seq_len(p)) {
    f_k <- base_freq + stats::runif(1, -0.5, 0.5)
    phase_k <- stats::runif(1, 0, 2 * pi)
    X[, k] <- cos(2 * pi * f_k * t + phase_k) +
              0.3 * stats::rnorm(n_time)
  }
  C <- stats::cov(X)
  # Mild ridge to guarantee strict positive definiteness across all anchors.
  C + 1e-3 * diag(p)
}

anchor_seeds_A <- 101:104
anchor_seeds_B <- 201:204

anchors_A <- lapply(anchor_seeds_A, function(s) build_anchor(10, s))
anchors_B <- lapply(anchor_seeds_B, function(s) build_anchor(6,  s))

example_anchors <- c(anchors_A, anchors_B)
example_labels <- c(rep(0L, n_per_class), rep(1L, n_per_class))

# Validate every anchor is SPD before saving.
for (i in seq_along(example_anchors)) {
  C <- example_anchors[[i]]
  stopifnot(isSymmetric(C, tol = 1e-8))
  ev <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  stopifnot(all(ev > 0))
}

# ---------------------------------------------------------------------------
# 2. Save lazy-load (data/) and file-based (inst/extdata/) copies.
# ---------------------------------------------------------------------------
save(example_anchors,
     file = file.path(pkg_root, "data", "example_anchors.rda"),
     compress = "xz")

saveRDS(example_anchors,
        file = file.path(pkg_root, "inst", "extdata", "example_anchors.rds"))
saveRDS(example_labels,
        file = file.path(pkg_root, "inst", "extdata", "example_labels.rds"))

# ---------------------------------------------------------------------------
# 3. Minimal calibration manifest CSV. Schema must match the columns
#    documented in R/manifest.R::read_calibration_manifest.
#
#    One demo cell:
#      dataset = "Example", subject = 1, budget_per_class = 4, resample_id = 0
#    Five method codes (R0/E0/G0/G1/G2), each split into class 0 and class 1.
#    Synthetic ratio n_aug = 3 per anchor. Real trial ids are 0-based to
#    follow the Python protocol convention; replay_from_manifest converts
#    them to 1-based R indices internally.
# ---------------------------------------------------------------------------
methods <- c("R0", "E0", "G0", "G1", "G2")
method_seeds <- c(R0 = 1000L, E0 = 1001L, G0 = 1002L,
                  G1 = 1003L, G2 = 1004L)
real_ids_A <- "[0, 1, 2, 3]"   # anchors_A occupy positions 1..4 in R
real_ids_B <- "[4, 5, 6, 7]"   # anchors_B occupy positions 5..8 in R

manifest_rows <- list()
row_idx <- 1L
for (m in methods) {
  for (cls in 0:1) {
    manifest_rows[[row_idx]] <- data.frame(
      dataset = "Example",
      subject = 1L,
      budget_per_class = 4L,
      resample_id = 0L,
      method_code = m,
      class_label = cls,
      n_real_trials = 4L,
      real_trial_ids = if (cls == 0) real_ids_A else real_ids_B,
      real_seed = 42L,
      method_seed = method_seeds[[m]],
      n_synthetic_per_anchor = if (m == "R0") 0L else 3L,
      target_session_ids_used = "False",
      status = "ok",
      notes = if (m == "R0") "real_only_no_synthetic" else "demo_cell",
      stringsAsFactors = FALSE
    )
    row_idx <- row_idx + 1L
  }
}
example_manifest <- do.call(rbind, manifest_rows)

write.csv(example_manifest,
          file = file.path(pkg_root, "inst", "extdata",
                           "example_manifest.csv"),
          row.names = FALSE)
# NB: keep the default quote = TRUE. The real_trial_ids column holds
# JSON-like strings such as "[0, 1, 2, 3]" whose internal commas would
# otherwise be parsed as field separators by read.csv on round-trip.

cat(sprintf("Wrote %d anchors (p = %d), %d manifest rows.\n",
            length(example_anchors), p, nrow(example_manifest)))
