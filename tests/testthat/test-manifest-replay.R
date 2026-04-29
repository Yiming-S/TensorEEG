helper_spd <- function(p, seed) {
  set.seed(seed)
  A <- matrix(stats::rnorm(p * p), p, p)
  A %*% t(A) + p * diag(p)
}

helper_write_manifest <- function(path, method, method_seed,
                                   real_ids_class0,
                                   real_ids_class1) {
  json_list <- function(v) {
    paste0("[", paste(as.integer(v), collapse = ", "), "]")
  }
  rows <- data.frame(
    dataset = c("DummySet", "DummySet"),
    subject = c(1L, 1L),
    budget_per_class = c(length(real_ids_class0),
                          length(real_ids_class0)),
    resample_id = c(0L, 0L),
    method_code = c(method, method),
    class_label = c(0L, 1L),
    n_real_trials = c(length(real_ids_class0),
                       length(real_ids_class1)),
    real_trial_ids = c(json_list(real_ids_class0),
                        json_list(real_ids_class1)),
    real_seed = c(123L, 123L),
    method_seed = c(method_seed, method_seed),
    n_synthetic_per_anchor = c(if (method == "R0") 0L else 3L,
                                if (method == "R0") 0L else 3L),
    target_session_ids_used = c("false", "false"),
    status = c("ok", "ok"),
    notes = c("", ""),
    stringsAsFactors = FALSE
  )
  write.csv(rows, path, row.names = FALSE)
}

# ---------------------------------------------------------------------
# read_calibration_manifest
# ---------------------------------------------------------------------

test_that("read_calibration_manifest parses real_trial_ids JSON", {
  tmp <- tempfile(fileext = ".csv")
  helper_write_manifest(tmp, "G0", method_seed = 42L,
                         real_ids_class0 = c(0L, 3L, 5L),
                         real_ids_class1 = c(1L, 4L, 7L))
  df <- read_calibration_manifest(tmp)
  expect_equal(nrow(df), 2L)
  expect_type(df$real_trial_ids, "list")
  expect_equal(df$real_trial_ids[[1L]], c(0L, 3L, 5L))
  expect_equal(df$real_trial_ids[[2L]], c(1L, 4L, 7L))
  unlink(tmp)
})

test_that("read_calibration_manifest fails on missing columns", {
  tmp <- tempfile(fileext = ".csv")
  write.csv(data.frame(dataset = "x"), tmp, row.names = FALSE)
  expect_error(read_calibration_manifest(tmp), "missing required")
  unlink(tmp)
})

# ---------------------------------------------------------------------
# replay_from_manifest dispatches the correct method
# ---------------------------------------------------------------------

test_that("replay_from_manifest reproduces G0 from manifest seed", {
  tmp <- tempfile(fileext = ".csv")
  helper_write_manifest(tmp, "G0", method_seed = 42L,
                         real_ids_class0 = c(0L, 1L, 2L),
                         real_ids_class1 = c(3L, 4L, 5L))
  manifest <- read_calibration_manifest(tmp)

  source_cov <- lapply(1:6, function(k) helper_spd(4L, 100L + k))
  source_labels <- c(0L, 0L, 0L, 1L, 1L, 1L)
  res <- replay_from_manifest(manifest, source_cov, source_labels,
                               ratio_per_anchor = 3L, g0_sigma = 0.10)
  expect_length(res$cov, 18L)
  expect_equal(res$labels, source_labels[res$anchor])
  unlink(tmp)
})

test_that("replay_from_manifest returns empty stack for R0 method", {
  tmp <- tempfile(fileext = ".csv")
  helper_write_manifest(tmp, "R0", method_seed = 42L,
                         real_ids_class0 = c(0L, 1L),
                         real_ids_class1 = c(2L, 3L))
  manifest <- read_calibration_manifest(tmp)

  source_cov <- lapply(1:4, function(k) helper_spd(3L, 200L + k))
  source_labels <- c(0L, 0L, 1L, 1L)
  res <- replay_from_manifest(manifest, source_cov, source_labels)
  expect_length(res$cov, 0L)
  expect_equal(length(res$anchors_used), 4L)
  unlink(tmp)
})

test_that("replay_from_manifest rejects A0 (transductive)", {
  tmp <- tempfile(fileext = ".csv")
  helper_write_manifest(tmp, "A0", method_seed = 42L,
                         real_ids_class0 = c(0L, 1L),
                         real_ids_class1 = c(2L, 3L))
  manifest <- read_calibration_manifest(tmp)

  source_cov <- lapply(1:4, function(k) helper_spd(3L, 300L + k))
  expect_error(
    replay_from_manifest(manifest, source_cov, c(0L, 0L, 1L, 1L)),
    "A0.*transductive"
  )
  unlink(tmp)
})

# ---------------------------------------------------------------------
# generate_drift_rotations: fBm option
# ---------------------------------------------------------------------

test_that("generate_drift_rotations produces orthogonal matrices under fbm", {
  set.seed(7L)
  R_list <- generate_drift_rotations(n_sources = 5L, n_trials = 10L,
                                      sigma_eps = 0.05,
                                      process = "fbm", hurst = 0.85)
  expect_length(R_list, 10L)
  for (R in R_list) {
    expect_equal(R %*% t(R), diag(nrow(R)), tolerance = 1e-8)
  }
})

test_that("generate_drift_rotations rejects invalid hurst", {
  expect_error(generate_drift_rotations(3L, 5L, process = "fbm",
                                         hurst = 0),
               "in \\(0, 1\\)")
  expect_error(generate_drift_rotations(3L, 5L, process = "fbm",
                                         hurst = 1.5),
               "in \\(0, 1\\)")
})
