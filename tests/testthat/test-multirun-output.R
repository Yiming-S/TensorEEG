test_that("multi-run wrapper returns consistent structure", {
  sess <- sim_multirun_session(
    n_runs = 2,
    trials_per_run = 3,
    gap_trials = 1,
    n_time = 60,
    n_channels = 8,
    n_sources = 3,
    fs = 128,
    seed = 123,
    verbose = FALSE
  )

  expect_length(sess$x, 6)
  expect_equal(length(sess$y), 6)
  expect_equal(length(sess$run), 6)
  expect_equal(sort(unique(as.integer(sess$run))), c(1, 2))
  expect_equal(sess$fs, 128)
})
