helper_random_spd <- function(p, seed = 1L) {
  set.seed(seed)
  A <- matrix(stats::rnorm(p * p), p, p)
  A %*% t(A) + p * diag(p)
}

helper_anchor_list <- function(n, p, seed = 1L) {
  lapply(seq_len(n), function(k) helper_random_spd(p, seed + k))
}

# ---------------------------------------------------------------------
# E0 amplitude-matched Euclidean perturbation
# ---------------------------------------------------------------------

test_that("augment_cov_amplitude_matched_euclidean returns SPD outputs", {
  anchors <- helper_anchor_list(8L, 6L, seed = 11L)
  res <- augment_cov_amplitude_matched_euclidean(
    anchors, n_aug = 3L, g0_sigma = 0.10,
    labels = rep(c(0L, 1L), each = 4L), seed = 42L)
  expect_length(res$cov, 24L)
  for (C in res$cov) {
    expect_equal(dim(C), c(6L, 6L))
    expect_equal(C, t(C), tolerance = 1e-10)
    expect_true(min(eigen(C, symmetric = TRUE,
                          only.values = TRUE)$values) > 0)
  }
  # rho diagnostic: should approach 1 (within tolerance) on this size.
  expect_true(abs(res$diagnostic$rho - 1) < 0.5)
})

test_that("augment_cov_amplitude_matched_euclidean inherits anchor labels", {
  anchors <- helper_anchor_list(6L, 4L, seed = 21L)
  labels <- c(0L, 0L, 0L, 1L, 1L, 1L)
  res <- augment_cov_amplitude_matched_euclidean(
    anchors, n_aug = 2L, g0_sigma = 0.10, labels = labels, seed = 7L)
  expect_equal(res$labels, labels[res$anchor])
})

# ---------------------------------------------------------------------
# G1 empirical tangent Gaussian
# ---------------------------------------------------------------------

test_that("augment_cov_empirical_tangent returns SPD outputs", {
  anchors <- helper_anchor_list(10L, 5L, seed = 31L)
  labels <- c(rep(0L, 5L), rep(1L, 5L))
  res <- augment_cov_empirical_tangent(anchors, labels,
                                       n_aug = 2L, sigma = 0.05,
                                       seed = 7L)
  expect_length(res$cov, 20L)
  for (C in res$cov) {
    expect_true(min(eigen(C, symmetric = TRUE,
                          only.values = TRUE)$values) > 0)
  }
  expect_equal(res$labels, labels[res$anchor])
})

test_that("augment_cov_empirical_tangent rejects mismatched labels", {
  anchors <- helper_anchor_list(4L, 4L, seed = 41L)
  expect_error(augment_cov_empirical_tangent(anchors, c(0L, 1L),
                                              n_aug = 1L),
               "must have length")
})

# ---------------------------------------------------------------------
# G2 log-Euclidean geodesic mixup
# ---------------------------------------------------------------------

test_that("augment_cov_geodesic_mixup picks same-class partners", {
  anchors <- helper_anchor_list(8L, 4L, seed = 51L)
  labels <- c(rep(0L, 4L), rep(1L, 4L))
  res <- augment_cov_geodesic_mixup(anchors, labels, n_aug = 3L,
                                    beta_alpha = 1.0, seed = 13L)
  expect_length(res$cov, 24L)
  expect_length(res$trace$partner, 24L)
  expect_length(res$trace$alpha, 24L)
  for (k in seq_along(res$cov)) {
    expect_equal(labels[res$anchor[k]], labels[res$trace$partner[k]])
    expect_true(min(eigen(res$cov[[k]], symmetric = TRUE,
                          only.values = TRUE)$values) > 0)
  }
})

# ---------------------------------------------------------------------
# A0 Riemannian alignment
# ---------------------------------------------------------------------

test_that("augment_cov_alignment_riemannian outputs SPD anchors", {
  source_cov <- helper_anchor_list(6L, 4L, seed = 61L)
  target_cov <- helper_anchor_list(8L, 4L, seed = 71L)
  res <- augment_cov_alignment_riemannian(source_cov, target_cov)
  expect_length(res$cov, 6L)
  for (C in res$cov) {
    expect_equal(C, t(C), tolerance = 1e-8)
    expect_true(min(eigen(C, symmetric = TRUE,
                          only.values = TRUE)$values) > 0)
  }
})

# ---------------------------------------------------------------------
# Budget-ratio warning on G0
# ---------------------------------------------------------------------

test_that("augment_cov_riemannian warns under high synthetic-to-real ratio", {
  anchors <- helper_anchor_list(20L, 4L, seed = 81L)   # below threshold 40
  expect_warning(
    augment_cov_riemannian(anchors, n_aug = 3L, sigma = 0.05, seed = 1L),
    "Synthetic-to-real ratio"
  )
})

test_that("augment_cov_riemannian does not warn under safe ratio", {
  anchors <- helper_anchor_list(60L, 4L, seed = 91L)   # above threshold 40
  expect_silent(
    augment_cov_riemannian(anchors, n_aug = 3L, sigma = 0.05, seed = 1L)
  )
})
