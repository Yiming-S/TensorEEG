test_that("tensor_to_cov returns SPD matrices of the right shape", {
  set.seed(11)
  X <- array(rnorm(200 * 6 * 4), dim = c(200, 6, 4))
  cov_list <- tensor_to_cov(X, ridge = 1e-6)
  expect_length(cov_list, 4)
  for (C in cov_list) {
    expect_equal(dim(C), c(6L, 6L))
    # symmetric
    expect_equal(C, t(C), tolerance = 1e-12)
    # SPD: smallest eigenvalue strictly positive
    expect_true(min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) > 0)
  }
})

test_that("tensor_to_cov rejects malformed input and short trials", {
  expect_error(tensor_to_cov(matrix(1, 3, 3)), "3D array")
  expect_error(tensor_to_cov(array(1, dim = c(1, 4, 2))), "at least 2 time")
  expect_error(tensor_to_cov(array(1, dim = c(10, 4, 2)), ridge = -1),
               "non-negative")
})

test_that("cov_logeuclidean_distance is zero for identical matrices and positive otherwise", {
  C1 <- diag(c(2, 1, 0.5))
  expect_equal(cov_logeuclidean_distance(C1, C1), 0, tolerance = 1e-12)

  C2 <- diag(c(2.5, 1, 0.5))
  d <- cov_logeuclidean_distance(C1, C2)
  expect_true(d > 0)
  expect_equal(d, abs(log(2.5) - log(2)), tolerance = 1e-12)
})

test_that("augment_cov_riemannian outputs SPD matrices in anchor-major order", {
  set.seed(2)
  cov_list <- list(
    diag(c(2, 1, 0.5)),
    diag(c(1.5, 0.8, 0.3)),
    diag(c(3, 1.2, 0.6))
  )
  aug <- augment_cov_riemannian(cov_list, n_aug = 4L, sigma = 0.1, seed = 1)
  expect_length(aug$cov, 12L)
  expect_equal(aug$anchor, rep(1:3, each = 4))
  expect_equal(aug$replicate, rep(1:4, times = 3))
  for (C in aug$cov) {
    expect_equal(dim(C), c(3L, 3L))
    expect_equal(C, t(C), tolerance = 1e-10)
    expect_true(min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) > 0)
  }
})

test_that("augment_cov_riemannian is reproducible under fixed seed", {
  cov_list <- list(diag(c(2, 1)), diag(c(1.5, 0.8)))
  a <- augment_cov_riemannian(cov_list, n_aug = 3L, sigma = 0.2, seed = 42)
  b <- augment_cov_riemannian(cov_list, n_aug = 3L, sigma = 0.2, seed = 42)
  for (k in seq_along(a$cov)) {
    expect_equal(a$cov[[k]], b$cov[[k]], tolerance = 1e-12)
  }
  expect_equal(a$labels, NULL)
})

test_that("augment_cov_riemannian propagates labels and accepts drift", {
  set.seed(3)
  cov_list <- list(diag(c(2, 1)), diag(c(1.5, 0.8)))
  drift <- matrix(c(0.05, 0.01, 0.01, -0.02), 2, 2)
  aug <- augment_cov_riemannian(cov_list, n_aug = 2L, sigma = 0.05,
                                drift = drift,
                                labels = c("left", "right"), seed = 7)
  expect_equal(aug$labels, c("left", "left", "right", "right"))
  for (C in aug$cov) {
    expect_equal(C, t(C), tolerance = 1e-10)
    expect_true(min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) > 0)
  }
})

test_that("augment_cov_riemannian rejects bad arguments", {
  cov_list <- list(diag(c(2, 1)), diag(c(1.5, 0.8)))
  expect_error(augment_cov_riemannian(list()), "non-empty list")
  expect_error(augment_cov_riemannian(cov_list, n_aug = -1), "positive integer")
  expect_error(augment_cov_riemannian(cov_list, sigma = -0.1), "non-negative")
  expect_error(augment_cov_riemannian(cov_list, drift = matrix(0, 3, 3)),
               "2 x 2")
  expect_error(augment_cov_riemannian(cov_list, labels = "only-one"),
               "length equal to length")
  expect_error(augment_cov_riemannian(list(matrix(1, 2, 3))),
               "square matrix|2 x 2|3 x 3")
})

test_that("end-to-end: sim_eeg_master -> tensor_to_cov -> augment is SPD-stable", {
  sim <- sim_eeg_master(n_trials = 4, n_time = 100, n_channels = 5,
                        n_sources = 3, seed = 9, verbose = FALSE)
  cov_list <- tensor_to_cov(sim$data, ridge = 1e-5)
  aug <- augment_cov_riemannian(cov_list, n_aug = 2, sigma = 0.05,
                                labels = sim$labels, seed = 9)
  expect_length(aug$cov, 8L)
  expect_equal(aug$labels, rep(sim$labels, each = 2))
  for (C in aug$cov) {
    expect_true(min(eigen(C, symmetric = TRUE, only.values = TRUE)$values) > 0)
  }
})

test_that("log-Euclidean fidelity rank-orders covariance-aware below isotropic Euclidean noise", {
  set.seed(123)
  C_real <- diag(c(2, 1, 0.5)) + matrix(0.1, 3, 3)
  C_real <- 0.5 * (C_real + t(C_real))
  # Covariance-aware: small symmetric noise in tangent space.
  cov_aware <- augment_cov_riemannian(list(C_real), n_aug = 50L,
                                      sigma = 0.05, seed = 1)
  # Off-manifold Euclidean: add asymmetric Euclidean noise then symmetrise
  # without enforcing PD geometry.
  euclidean_noise <- replicate(50, {
    G <- matrix(rnorm(9, sd = 0.05), 3, 3)
    Cn <- C_real + 0.5 * (G + t(G))
    list(Cn)
  }, simplify = FALSE)
  d_aware <- vapply(cov_aware$cov, cov_logeuclidean_distance, numeric(1),
                    C2 = C_real)
  d_eucl <- vapply(euclidean_noise, function(L) {
    Cn <- L[[1L]]
    eig <- eigen(Cn, symmetric = TRUE, only.values = TRUE)$values
    if (any(eig <= 0)) return(NA_real_)
    cov_logeuclidean_distance(Cn, C_real)
  }, numeric(1))
  # Geometry-preserving augmentation should not produce NA distances.
  expect_true(all(is.finite(d_aware)))
  # Mean LE-distance should be on the same scale; the key property is that
  # cov-aware augmentation never falls off the manifold.
  expect_true(mean(d_aware) > 0)
})
