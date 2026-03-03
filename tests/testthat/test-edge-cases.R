test_that("single-source drift rotations are identity", {
  rots <- generate_drift_rotations(n_sources = 1, n_trials = 3)

  expect_length(rots, 3)
  for (r in rots) {
    expect_equal(dim(r), c(1, 1))
    expect_equal(r[1, 1], 1)
  }
})

test_that("single-source VAR2 simulation preserves matrix shape", {
  vp <- setup_var2_system(n_sources = 1, fs = 128, target_freqs = 10)
  x <- sim_source_var2(n_time = 50, n_sources = 1, var_params = vp)

  expect_equal(dim(x), c(50, 1))
  expect_true(all(is.finite(x)))
})

test_that("sim_eeg_master supports single-source simulation", {
  sim <- sim_eeg_master(
    n_trials = 4,
    n_time = 80,
    n_channels = 8,
    n_sources = 1,
    fs = 128,
    target_freqs = 10,
    seed = 1,
    verbose = FALSE
  )

  expect_equal(dim(sim$data), c(80, 8, 4))
  expect_equal(nrow(sim$audit), 4)
  expect_true(all(is.finite(sim$data)))
})

test_that("geometry generation handles single-channel edge case", {
  geo <- generate_geometry_mixing(n_channels = 1, n_sources = 2)

  expect_equal(dim(geo$A_base), c(1, 2))
  expect_equal(dim(geo$L_sym), c(1, 1))
  expect_true(all(is.finite(geo$A_base)))
})
