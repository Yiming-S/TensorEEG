test_that("sim_eeg_master rejects non-integer trials", {
  expect_error(
    sim_eeg_master(
      n_trials = 2.5,
      n_time = 40,
      n_channels = 8,
      n_sources = 3,
      fs = 128,
      verbose = FALSE
    ),
    "n_trials must be a positive integer"
  )
})

test_that("sim_eeg_master rejects invalid class labels", {
  expect_error(
    sim_eeg_master(
      n_trials = 3,
      n_time = 40,
      n_channels = 8,
      n_sources = 3,
      fs = 128,
      class_labels = c(0, 1, 0.2),
      verbose = FALSE
    ),
    "class_labels must only contain integer 0 and 1 values"
  )
})

test_that("sim_multirun_session validates run arguments", {
  expect_error(
    sim_multirun_session(n_runs = 0, verbose = FALSE),
    "n_runs must be a positive integer"
  )
  expect_error(
    sim_multirun_session(n_runs = 2, trials_per_run = 3, gap_trials = -1, verbose = FALSE),
    "gap_trials must be a non-negative integer"
  )
})

test_that("calc_ac_power validates matrix input", {
  expect_error(
    calc_ac_power(c(1, 2, 3), fs = 250),
    "X must be a numeric matrix"
  )
})
