#' EEG Topomap Interpolation Helper
#'
#' @param x Numeric vector. Sensor x coordinates.
#' @param y Numeric vector. Sensor y coordinates.
#' @param v Numeric vector. Sensor values.
#' @param n_grid Integer. Grid resolution.
#' @param sigma Numeric (optional). Kernel width.
#'
#' @return List with interpolated grid components.
.interpolate_topomap <- function(x, y, v, n_grid = 90, sigma = NULL) {
  x_seq <- seq(-1, 1, length.out = n_grid)
  y_seq <- seq(-1, 1, length.out = n_grid)
  gx <- rep(x_seq, each = n_grid)
  gy <- rep(y_seq, times = n_grid)
  mask <- gx^2 + gy^2 <= 1
  
  if(is.null(sigma)) {
    d_mat <- as.matrix(stats::dist(cbind(x, y)))
    sigma <- stats::median(d_mat[d_mat > 0], na.rm = TRUE)
    if(!is.finite(sigma) || sigma <= 0) sigma <- 0.2
  }
  
  gxy <- cbind(gx[mask], gy[mask])
  d2 <- outer(gxy[, 1], x, "-")^2 + outer(gxy[, 2], y, "-")^2
  w <- exp(-d2 / (2 * sigma^2))
  w_sum <- rowSums(w)
  w_sum[w_sum == 0] <- 1
  z_vals <- as.numeric((w %*% v) / w_sum)
  
  z_mat <- matrix(NA_real_, nrow = n_grid, ncol = n_grid)
  z_mat[mask] <- z_vals
  
  list(x = x_seq, y = y_seq, z = z_mat, sigma = sigma)
}

#' Plot EEG Topomap from Simulated Trial
#'
#' @description
#' Visualizes a scalp topography from one simulated trial at a selected time point
#' or over a selected time window.
#'
#' @param sim_res List output from \code{\link{sim_eeg_master}}.
#' @param trial Integer. Trial index.
#' @param time_idx Integer (optional). Sample index for topomap snapshot.
#' @param time_window_ms Numeric vector length 2 (optional). Time window in milliseconds.
#'   If provided, channel values are averaged in this window.
#' @param palette Character vector of colors for interpolation heatmap.
#' @param n_grid Integer. Grid resolution for interpolation.
#' @param sigma Numeric (optional). Interpolation kernel width.
#' @param show_sensors Logical. Whether to draw sensor positions.
#' @param main Character (optional). Plot title.
#' @param save_path Character (optional). If provided, saves plot to PNG.
#'
#' @return Invisible list containing selected trial/time metadata and topography values.
#' @export
plot_topomap <- function(sim_res,
                         trial = 1,
                         time_idx = NULL,
                         time_window_ms = NULL,
                         palette = grDevices::hcl.colors(32, "Spectral", rev = TRUE),
                         n_grid = 90,
                         sigma = NULL,
                         show_sensors = TRUE,
                         main = NULL,
                         save_path = NULL) {
  if(is.null(sim_res$data) || is.null(sim_res$geometry$coords_sens)) {
    stop("sim_res must contain $data and $geometry$coords_sens.")
  }
  
  X <- sim_res$data
  coords <- sim_res$geometry$coords_sens
  fs <- sim_res$params$fs
  
  dims <- dim(X)
  if(length(dims) != 3) stop("sim_res$data must be [time x channel x trial].")
  n_time <- dims[1]
  n_channels <- dims[2]
  n_trials <- dims[3]
  
  trial <- as.integer(trial[1])
  if(trial < 1 || trial > n_trials) stop("trial index out of range.")
  
  if(nrow(coords) != n_channels || ncol(coords) < 2) {
    stop("coords_sens must have one row per channel and at least 2 columns.")
  }
  
  if(!is.null(save_path)) {
    grDevices::png(filename = save_path, width = 900, height = 850, res = 140)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  if(!is.null(time_window_ms)) {
    if(length(time_window_ms) != 2) stop("time_window_ms must have length 2.")
    time_window_ms <- sort(as.numeric(time_window_ms))
    idx_start <- max(1L, min(n_time, round(time_window_ms[1] / 1000 * fs)))
    idx_end <- max(1L, min(n_time, round(time_window_ms[2] / 1000 * fs)))
    if(idx_end < idx_start) {
      tmp <- idx_start
      idx_start <- idx_end
      idx_end <- tmp
    }
    topo_vec <- colMeans(X[idx_start:idx_end, , trial])
    time_label <- sprintf("%d-%d ms",
                          round((idx_start - 1) / fs * 1000),
                          round((idx_end - 1) / fs * 1000))
  } else {
    if(is.null(time_idx)) time_idx <- round(n_time / 2)
    time_idx <- as.integer(time_idx[1])
    time_idx <- max(1L, min(n_time, time_idx))
    topo_vec <- as.numeric(X[time_idx, , trial])
    time_label <- sprintf("%d ms", round((time_idx - 1) / fs * 1000))
  }
  
  x <- coords[, 1]
  y <- coords[, 2]
  interp <- .interpolate_topomap(x, y, topo_vec, n_grid = n_grid, sigma = sigma)
  
  if(is.null(main)) {
    main <- sprintf("Topomap Trial %d (%s)", trial, time_label)
  }
  
  graphics::image(interp$x, interp$y, interp$z,
                  col = palette, axes = FALSE, asp = 1,
                  xlab = "", ylab = "", main = main, useRaster = TRUE)
  
  graphics::contour(interp$x, interp$y, interp$z,
                    add = TRUE, drawlabels = FALSE, nlevels = 8,
                    col = grDevices::adjustcolor("black", alpha.f = 0.35))
  
  theta <- seq(0, 2 * pi, length.out = 300)
  graphics::lines(cos(theta), sin(theta), lwd = 2)
  graphics::lines(c(-0.08, 0, 0.08), c(1.0, 1.1, 1.0), lwd = 2)
  graphics::lines(c(-1.02, -1.08, -1.02), c(0.16, 0, -0.16), lwd = 2)
  graphics::lines(c(1.02, 1.08, 1.02), c(0.16, 0, -0.16), lwd = 2)
  
  if(show_sensors) {
    graphics::points(x, y, pch = 16, cex = 0.45, col = "gray20")
  }
  
  invisible(list(
    trial = trial,
    time_label = time_label,
    sigma = interp$sigma,
    topography = topo_vec
  ))
}

#' Plot Run-Wise Drift Diagnostics
#'
#' @description
#' Computes a per-trial drift metric and visualizes the shift across runs from
#' \code{\link{sim_multirun_session}} output.
#'
#' @param session_res List output from \code{\link{sim_multirun_session}}.
#' @param metric Character. One of \code{"bandpower"}, \code{"cov_distance"}, \code{"mean_abs"}.
#' @param band Numeric vector length 2. Frequency band for \code{metric = "bandpower"}.
#' @param channel Integer (optional). If set, metric is computed on one channel.
#'   Otherwise, channel-averaged series is used.
#' @param main Character (optional). Plot title.
#' @param save_path Character (optional). If provided, saves plot to PNG.
#'
#' @return Invisible data frame with trial-level metric values.
#' @export
plot_run_drift <- function(session_res,
                           metric = c("bandpower", "cov_distance", "mean_abs"),
                           band = c(8, 30),
                           channel = NULL,
                           main = NULL,
                           save_path = NULL) {
  metric <- match.arg(metric)
  
  if(is.null(session_res$x) || is.null(session_res$run) || is.null(session_res$fs)) {
    stop("session_res must contain $x, $run, and $fs.")
  }
  
  x_list <- session_res$x
  runs <- as.factor(session_res$run)
  fs <- session_res$fs
  n_trials <- length(x_list)
  
  if(length(runs) != n_trials) {
    stop("length(session_res$run) must match length(session_res$x).")
  }
  
  if(n_trials == 0) stop("session_res$x is empty.")
  
  n_channels <- ncol(x_list[[1]])
  if(!is.null(channel)) {
    channel <- as.integer(channel[1])
    if(channel < 1 || channel > n_channels) stop("channel index out of range.")
  }
  
  if(!is.null(save_path)) {
    grDevices::png(filename = save_path, width = 1100, height = 800, res = 140)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  trial_metric <- rep(NA_real_, n_trials)
  
  if(metric == "bandpower") {
    band <- sort(as.numeric(band))
    for(i in seq_len(n_trials)) {
      X <- x_list[[i]]
      x_series <- if(is.null(channel)) rowMeans(X) else X[, channel]
      sp <- stats::spectrum(x_series, plot = FALSE)
      f_hz <- sp$freq * fs
      idx <- which(f_hz >= band[1] & f_hz <= band[2])
      trial_metric[i] <- if(length(idx) > 0) {
        mean(10 * log10(pmax(sp$spec[idx], .Machine$double.eps)))
      } else {
        NA_real_
      }
    }
    y_label <- sprintf("Bandpower %.1f-%.1f Hz (dB)", band[1], band[2])
  }
  
  if(metric == "mean_abs") {
    for(i in seq_len(n_trials)) {
      X <- x_list[[i]]
      x_series <- if(is.null(channel)) rowMeans(X) else X[, channel]
      trial_metric[i] <- mean(abs(x_series))
    }
    y_label <- "Mean Absolute Amplitude"
  }
  
  if(metric == "cov_distance") {
    bf <- signal::butter(4, 1/(fs/2), type = "high")
    run_levels <- levels(runs)
    ref_idx <- which(runs == run_levels[1])
    
    cov_ref <- matrix(0, n_channels, n_channels)
    for(j in ref_idx) {
      X_hp <- apply(x_list[[j]], 2, function(x) signal::filtfilt(bf, x))
      cov_ref <- cov_ref + stats::cov(X_hp)
    }
    cov_ref <- cov_ref / max(1, length(ref_idx))
    
    for(i in seq_len(n_trials)) {
      X_hp <- apply(x_list[[i]], 2, function(x) signal::filtfilt(bf, x))
      cov_i <- stats::cov(X_hp)
      trial_metric[i] <- sqrt(sum((cov_i - cov_ref)^2))
    }
    y_label <- "Covariance Distance to Run 1"
  }
  
  drift_df <- data.frame(
    trial = seq_len(n_trials),
    run = runs,
    value = trial_metric
  )
  
  if(is.null(main)) {
    main <- sprintf("Run Drift Diagnostic (%s)", metric)
  }
  
  box_cols <- grDevices::hcl.colors(nlevels(runs), "Set 2")
  graphics::boxplot(value ~ run, data = drift_df,
                    col = box_cols, xlab = "Run", ylab = y_label, main = main)
  
  run_mean <- tapply(drift_df$value, drift_df$run, mean, na.rm = TRUE)
  graphics::lines(seq_along(run_mean), run_mean, type = "b", pch = 19, lwd = 2, col = "black")
  
  invisible(drift_df)
}

#' Plot Simulated Artifact Components
#'
#' @description
#' Visualizes fast artifacts and low-frequency drift generated by
#' \code{\link{sim_artifacts}}.
#'
#' @param artifact_res List output from \code{\link{sim_artifacts}} containing
#'   \code{N_fast} and \code{N_drift}.
#' @param fs Numeric. Sampling frequency.
#' @param channel Integer. Channel index to inspect.
#' @param fmax Numeric. Maximum frequency for PSD subplots.
#' @param save_path Character (optional). If provided, saves plot to PNG.
#'
#' @return Invisible list with RMS summary.
#' @export
plot_artifact_components <- function(artifact_res,
                                     fs,
                                     channel = 1,
                                     fmax = 60,
                                     save_path = NULL) {
  if(is.null(artifact_res$N_fast) || is.null(artifact_res$N_drift)) {
    stop("artifact_res must contain N_fast and N_drift.")
  }
  
  N_fast <- artifact_res$N_fast
  N_drift <- artifact_res$N_drift
  
  if(!all(dim(N_fast) == dim(N_drift))) {
    stop("N_fast and N_drift must have the same dimensions.")
  }
  
  n_time <- nrow(N_fast)
  n_channels <- ncol(N_fast)
  if(channel < 1 || channel > n_channels) stop("channel index out of range.")
  
  if(!is.null(save_path)) {
    grDevices::png(filename = save_path, width = 1300, height = 900, res = 140)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  t_vec <- (seq_len(n_time) - 1) / fs
  x_fast <- N_fast[, channel]
  x_drift <- N_drift[, channel]
  
  graphics::plot(t_vec, x_fast, type = "l", col = "firebrick",
                 xlab = "s", ylab = "Amplitude", main = sprintf("Fast Artifact (Ch%d)", channel))
  
  graphics::plot(t_vec, x_drift, type = "l", col = "steelblue",
                 xlab = "s", ylab = "Amplitude", main = sprintf("Drift Artifact (Ch%d)", channel))
  
  sp_fast <- stats::spectrum(x_fast, plot = FALSE)
  f_fast <- sp_fast$freq * fs
  psd_fast <- 10 * log10(pmax(sp_fast$spec, .Machine$double.eps))
  graphics::plot(f_fast, psd_fast, type = "l",
                 xlim = c(0, min(fmax, max(f_fast))), col = "firebrick",
                 xlab = "Hz", ylab = "dB", main = "Fast Artifact PSD")
  
  sp_drift <- stats::spectrum(x_drift, plot = FALSE)
  f_drift <- sp_drift$freq * fs
  psd_drift <- 10 * log10(pmax(sp_drift$spec, .Machine$double.eps))
  graphics::plot(f_drift, psd_drift, type = "l",
                 xlim = c(0, min(fmax, max(f_drift))), col = "steelblue",
                 xlab = "Hz", ylab = "dB", main = "Drift Artifact PSD")
  
  metrics <- list(
    channel = channel,
    rms_fast = sqrt(mean(x_fast^2)),
    rms_drift = sqrt(mean(x_drift^2))
  )
  
  invisible(metrics)
}

#' Unified Visualization Dashboard Entry
#'
#' @description
#' Dispatches to the corresponding visualization function by type.
#'
#' @param obj Simulation or artifact object.
#' @param type Character. One of \code{"validation"}, \code{"topomap"},
#'   \code{"run_drift"}, or \code{"artifact"}.
#' @param ... Additional arguments passed to target plotting function.
#'
#' @return Output of the selected plotting function.
#' @export
plot_sim_dashboard <- function(obj,
                               type = c("validation", "topomap", "run_drift", "artifact"),
                               ...) {
  type <- match.arg(type)
  
  if(type == "validation") {
    return(validate_sim_eeg(obj, ...))
  }
  if(type == "topomap") {
    return(plot_topomap(obj, ...))
  }
  if(type == "run_drift") {
    return(plot_run_drift(obj, ...))
  }
  if(type == "artifact") {
    return(plot_artifact_components(obj, ...))
  }
  
  invisible(NULL)
}
