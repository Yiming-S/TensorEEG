helper_spd <- function(p, seed) {
  set.seed(seed)
  A <- matrix(stats::rnorm(p * p), p, p)
  A %*% t(A) + p * diag(p)
}

# ---------------------------------------------------------------------
# Distance metrics
# ---------------------------------------------------------------------

test_that("cov_affine_invariant_distance is zero for identical SPD matrices", {
  C <- helper_spd(5L, 1L)
  expect_equal(cov_affine_invariant_distance(C, C), 0, tolerance = 1e-10)
})

test_that("cov_affine_invariant_distance is positive for distinct matrices", {
  C1 <- helper_spd(4L, 11L)
  C2 <- helper_spd(4L, 12L)
  expect_true(cov_affine_invariant_distance(C1, C2) > 0)
})

test_that("cov_eigenvalue_correlation is 1 for identical spectra", {
  C <- helper_spd(6L, 2L)
  expect_equal(cov_eigenvalue_correlation(C, C), 1, tolerance = 1e-12)
})

test_that("cov_trace_ratio is 1 for identical matrices", {
  C <- helper_spd(4L, 3L)
  expect_equal(cov_trace_ratio(C, C), 1, tolerance = 1e-12)
})

test_that("cov_condition_ratio is 1 for identical matrices", {
  C <- helper_spd(4L, 4L)
  expect_equal(cov_condition_ratio(C, C), 1, tolerance = 1e-10)
})

test_that("cov_anchor_perturbation_distance equals 0 when synthetic = anchor", {
  anchors <- list(helper_spd(4L, 5L), helper_spd(4L, 6L))
  d <- cov_anchor_perturbation_distance(
    anchors,
    synthetic = anchors,
    anchor_idx = c(1L, 2L)
  )
  expect_equal(d, 0, tolerance = 1e-10)
})

# ---------------------------------------------------------------------
# Six-metric audit wrapper
# ---------------------------------------------------------------------

test_that("audit_covariance_fidelity returns all expected fields", {
  anchors <- lapply(1:5, function(k) helper_spd(6L, 100L + k))
  syn <- lapply(1:5, function(k) helper_spd(6L, 200L + k))
  out <- audit_covariance_fidelity(anchors, syn,
                                    anchor_idx = c(1L, 2L, 3L, 4L, 5L))
  expect_named(out, c("n_synthetic", "log_euclidean_to_reference",
                       "affine_invariant_to_reference",
                       "eigenvalue_correlation", "trace_ratio",
                       "condition_number_ratio",
                       "anchor_perturbation_distance",
                       "normalized"))
  expect_equal(out$n_synthetic, 5L)
  expect_true(is.finite(out$log_euclidean_to_reference))
  expect_true(is.finite(out$affine_invariant_to_reference))
  expect_true(out$normalized$feature_dim == 21L)   # p=6 -> 21
  expect_equal(out$normalized$log_euclidean,
               out$log_euclidean_to_reference / sqrt(21))
})

test_that("audit_covariance_fidelity handles empty synthetic stack", {
  anchors <- list(helper_spd(4L, 1L))
  out <- audit_covariance_fidelity(anchors, list(), integer(0))
  expect_equal(out$n_synthetic, 0L)
  expect_true(is.na(out$log_euclidean_to_reference))
})

# ---------------------------------------------------------------------
# Mechanism check: E0 should be farther from class-mean than G0 on
# average for the same anchor set under matched amplitude. This is the
# manuscript's main mechanism finding (Shen & Degras 2026, Table 6).
# ---------------------------------------------------------------------

test_that("E0 has larger AI distance to class-mean than G0 on average", {
  anchors <- lapply(1:12, function(k) helper_spd(5L, 300L + k))
  labels <- rep(c(0L, 1L), each = 6L)
  # Suppress the budget-ratio warning: this test verifies a mechanism
  # claim, not the warning itself; the warning is exercised separately
  # in test-audit-family.R.
  g0 <- suppressWarnings(augment_cov_riemannian(
    anchors, n_aug = 3L, sigma = 0.10, labels = labels, seed = 1L))
  e0 <- augment_cov_amplitude_matched_euclidean(anchors, n_aug = 3L,
                                                 g0_sigma = 0.10,
                                                 labels = labels, seed = 2L)
  ref <- 0.5 * (Reduce("+", anchors) / length(anchors) +
                  t(Reduce("+", anchors) / length(anchors)))
  d_g0 <- mean(vapply(g0$cov,
                      function(C) cov_affine_invariant_distance(C, ref),
                      numeric(1)))
  d_e0 <- mean(vapply(e0$cov,
                      function(C) cov_affine_invariant_distance(C, ref),
                      numeric(1)))
  expect_true(d_e0 > d_g0)
})
