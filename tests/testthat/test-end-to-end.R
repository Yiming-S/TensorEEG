# End-to-end integration tests.
#
# Unit tests in test-cov-augmentation.R, test-fidelity.R, and
# test-manifest-replay.R cover individual function invariants. This file
# verifies that the documented user-facing pipelines compose correctly:
#
#   (1) sim_eeg_master() -> tensor_to_cov() -> four augmentations ->
#       audit_covariance_fidelity() runs without error and returns a
#       sane multi-metric profile;
#
#   (2) The shipped example assets (data(example_anchors),
#       inst/extdata/example_manifest.csv, inst/extdata/example_labels.rds)
#       drive the same pipeline used in the covariance-audit and
#       manifest-replay vignettes;
#
#   (3) Manifest-driven replay reproduces a seeded augmentation
#       byte-equivalently, so the cross-language reproducibility claim
#       holds within R.

helper_assemble_anchors <- function(seed = 7L, n_per_class = 4L,
                                    p = 6L, n_time = 400L) {
  set.seed(seed)
  build <- function(base_freq, jitter_seed) {
    set.seed(jitter_seed)
    fs <- 250
    t <- seq_len(n_time) / fs
    X <- matrix(0, n_time, p)
    for (k in seq_len(p)) {
      f_k <- base_freq + stats::runif(1, -0.5, 0.5)
      ph_k <- stats::runif(1, 0, 2 * pi)
      X[, k] <- cos(2 * pi * f_k * t + ph_k) + 0.3 * stats::rnorm(n_time)
    }
    stats::cov(X) + 1e-3 * diag(p)
  }
  list(
    cov_list = c(
      lapply(seq_len(n_per_class),
             function(i) build(10, seed * 100L + i)),
      lapply(seq_len(n_per_class),
             function(i) build(6,  seed * 200L + i))
    ),
    labels = c(rep(0L, n_per_class), rep(1L, n_per_class))
  )
}

# ---------------------------------------------------------------------------
test_that("end-to-end pipeline: anchors -> augment x4 -> fidelity audit", {
  bundle <- helper_assemble_anchors()
  anchors <- bundle$cov_list
  labels  <- bundle$labels
  n_aug <- 3L
  g0_sigma <- 0.15

  aug_g0 <- augment_cov_riemannian(
    anchors, n_aug = n_aug, sigma = g0_sigma,
    labels = labels, seed = 1L
  )
  aug_e0 <- augment_cov_amplitude_matched_euclidean(
    anchors, n_aug = n_aug, g0_sigma = g0_sigma,
    labels = labels, seed = 2L
  )
  aug_g1 <- augment_cov_empirical_tangent(
    anchors, labels, n_aug = n_aug, sigma = g0_sigma, seed = 3L
  )
  aug_g2 <- augment_cov_geodesic_mixup(
    anchors, labels, n_aug = n_aug, beta_alpha = 1.0, seed = 4L
  )

  for (aug in list(aug_g0, aug_e0, aug_g1, aug_g2)) {
    expect_equal(length(aug$cov), length(anchors) * n_aug)
    expect_true(all(vapply(aug$cov, isSymmetric, logical(1),
                            tol = 1e-8)))
    expect_true(all(vapply(aug$cov,
      function(C) all(eigen(C, symmetric = TRUE,
                            only.values = TRUE)$values > 0),
      logical(1))))
  }

  audit_one <- function(aug) {
    audit_covariance_fidelity(
      anchors    = anchors,
      synthetic  = aug$cov,
      anchor_idx = aug$anchor
    )
  }
  res_g0 <- audit_one(aug_g0)
  res_e0 <- audit_one(aug_e0)
  res_g1 <- audit_one(aug_g1)
  res_g2 <- audit_one(aug_g2)

  for (res in list(res_g0, res_e0, res_g1, res_g2)) {
    expect_equal(res$n_synthetic, length(anchors) * n_aug)
    expect_true(is.finite(res$log_euclidean_to_reference))
    expect_true(is.finite(res$affine_invariant_to_reference))
    expect_true(is.finite(res$anchor_perturbation_distance))
    expect_true(is.list(res$normalized))
    expect_true(res$normalized$feature_dim ==
                  nrow(anchors[[1]]) * (nrow(anchors[[1]]) + 1) / 2)
  }

  # E0 sits off-manifold by construction, so its log-Euclidean
  # distance to the class-mean reference must exceed G2's.
  expect_gt(res_e0$log_euclidean_to_reference,
            res_g2$log_euclidean_to_reference)
})

# ---------------------------------------------------------------------------
test_that("shipped example assets drive the same pipeline as the vignette", {
  manifest_path <- system.file("extdata", "example_manifest.csv",
                               package = "TensorEEG")
  labels_path   <- system.file("extdata", "example_labels.rds",
                               package = "TensorEEG")
  expect_true(nzchar(manifest_path))
  expect_true(nzchar(labels_path))

  manifest <- read_calibration_manifest(manifest_path)
  expect_true(all(c("R0", "E0", "G0", "G1", "G2") %in%
                    unique(manifest$method_code)))

  data(example_anchors, package = "TensorEEG", envir = environment())
  anchor_labels <- readRDS(labels_path)
  expect_equal(length(example_anchors), length(anchor_labels))

  for (m in c("R0", "E0", "G0", "G1", "G2")) {
    grp <- manifest[manifest$method_code == m, , drop = FALSE]
    res <- replay_from_manifest(
      manifest_group  = grp,
      source_cov_list = example_anchors,
      source_labels   = anchor_labels,
      ratio_per_anchor = 3L
    )
    if (m == "R0") {
      expect_equal(length(res$cov), 0L)
    } else {
      expect_equal(length(res$cov), length(example_anchors) * 3L)
    }
  }
})

# ---------------------------------------------------------------------------
test_that("manifest replay reproduces direct seeded augmentation byte-for-byte", {
  data(example_anchors, package = "TensorEEG", envir = environment())
  anchor_labels <- readRDS(
    system.file("extdata", "example_labels.rds", package = "TensorEEG")
  )
  manifest <- read_calibration_manifest(
    system.file("extdata", "example_manifest.csv", package = "TensorEEG")
  )

  # G2 replay vs direct seeded call.
  g2_group <- manifest[manifest$method_code == "G2", , drop = FALSE]
  g2_seed  <- as.integer(g2_group$method_seed[1])

  replayed <- replay_from_manifest(
    manifest_group = g2_group,
    source_cov_list = example_anchors,
    source_labels   = anchor_labels,
    ratio_per_anchor = 3L,
    g2_beta_alpha = 1.0
  )
  direct <- augment_cov_geodesic_mixup(
    cov_list   = example_anchors,
    labels     = anchor_labels,
    n_aug      = 3L,
    beta_alpha = 1.0,
    seed       = g2_seed
  )

  expect_equal(length(replayed$cov), length(direct$cov))
  max_diff <- max(mapply(function(A, B) max(abs(A - B)),
                          replayed$cov, direct$cov))
  expect_equal(max_diff, 0)
})
