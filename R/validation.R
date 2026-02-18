#' Validate Simulation Output via Visual Inspection
#'
#' @description
#' Generates a diagnostic dashboard (2x2 grid) to visually verify the spectral, spatial,
#' and temporal properties of the generated EEG tensor. This function is essential for
#' confirming that the simulation parameters (SNR, frequency peaks, spatial smoothness)
#' were correctly realized in the output.
#'
#' @details
#' The function produces four key validation plots:
#' \enumerate{
#'   \item \strong{Temporal Contrast (Top-Left):} Overlays single-trial traces from Class 0 vs. Class 1
#'     (after 1Hz High-Pass filtering). This verifies the amplitude difference defined by the class logic
#'     (e.g., Class 1 should show higher task amplitude).
#'   \item \strong{Power Spectral Density (Top-Right):} Computes the PSD of a Class 1 trial.
#'     Verifies the presence of the 20Hz Task peak (red dashed line) and the specific background
#'     oscillation frequencies (grey dotted lines) defined in the VAR(2) model.
#'   \item \strong{Spatial Covariance (Bottom-Left):} Displays the sensor covariance matrix.
#'     A "smooth" diagonal structure confirms that the Laplacian regularization successfully
#'     simulated volume conduction (neighboring sensors are highly correlated).
#'   \item \strong{SNR Audit (Bottom-Right):} Plots the "Realized Neural SNR" across all trials.
#'     Verifies that the closed-loop gain control system successfully maintained the target SNR
#'     (e.g., 5dB) despite random fluctuations in signal power.
#' }
#'
#' @param sim_res List. The simulation object returned by \code{\link{sim_eeg_master}}, containing
#'   \code{$data}, \code{$params}, \code{$labels}, and \code{$audit}.
#' @param class0_trial Integer (optional). Trial index to use as Class 0 representative.
#'   If NULL, the first available class 0 trial is used.
#' @param class1_trial Integer (optional). Trial index to use as Class 1 representative.
#'   If NULL, the first available class 1 trial is used.
#' @param channel Integer. Channel index used for temporal and PSD inspection.
#' @param fmax Numeric. Maximum frequency (Hz) for PSD plot.
#' @param snr_target_db Numeric (optional). Target SNR line for audit plot.
#'   If NULL, median realized SNR is used.
#' @param save_path Character (optional). If provided, saves the dashboard to PNG file.
#' @param return_metrics Logical. If TRUE, returns a metrics list instead of only plotting.
#'
#' @return Invisible metrics list. If \code{return_metrics = TRUE}, the metrics list
#'   is returned explicitly.
#'
#' @importFrom stats spectrum cov
#' @importFrom graphics par plot lines legend abline image
#' @importFrom signal butter filtfilt
#' @export
validate_sim_eeg <- function(sim_res,
                             class0_trial = NULL,
                             class1_trial = NULL,
                             channel = 1,
                             fmax = 60,
                             snr_target_db = NULL,
                             save_path = NULL,
                             return_metrics = FALSE) {
  
  if(is.null(sim_res$data) || length(dim(sim_res$data)) != 3) {
    stop("sim_res$data must be a 3D array [time x channel x trial].")
  }
  if(is.null(sim_res$params$fs)) {
    stop("sim_res$params$fs is required.")
  }
  
  X_tensor <- sim_res$data
  fs <- sim_res$params$fs
  target_freqs <- sim_res$params$target_freqs
  dims <- dim(X_tensor)
  n_time <- dims[1]
  n_channels <- dims[2]
  n_trials <- dims[3]
  
  if(!is.numeric(channel) || length(channel) != 1 || channel < 1 || channel > n_channels) {
    stop("channel must be a valid index in [1, n_channels].")
  }
  channel <- as.integer(channel)
  
  if(!is.null(save_path)) {
    grDevices::png(filename = save_path, width = 1400, height = 1000, res = 140)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  
  # Validation Filter: 1Hz HPF (Standard EEG Preproc)
  bf_val <- signal::butter(4, 1/(fs/2), type="high")
  
  # Plot 1: Temporal (Class Comparison)
  # Average Trial for Class 0 vs Class 1
  cls <- sim_res$labels
  if(is.null(cls)) cls <- rep(0, n_trials)
  
  if(is.null(class0_trial)) {
    idx0_candidates <- which(cls == 0)
    idx0 <- if(length(idx0_candidates) > 0) idx0_candidates[1] else 1L
  } else {
    idx0 <- as.integer(class0_trial[1])
  }
  
  if(is.null(class1_trial)) {
    idx1_candidates <- which(cls == 1)
    idx1 <- if(length(idx1_candidates) > 0) idx1_candidates[1] else min(2L, n_trials)
  } else {
    idx1 <- as.integer(class1_trial[1])
  }
  
  idx0 <- max(1L, min(n_trials, idx0))
  idx1 <- max(1L, min(n_trials, idx1))
  
  x_c0 <- signal::filtfilt(bf_val, X_tensor[, channel, idx0])
  x_c1 <- signal::filtfilt(bf_val, X_tensor[, channel, idx1])
  t_vec <- (1:length(x_c0)) / fs
  
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  plot(t_vec, x_c1, type="l", col="red",
       main=sprintf("Class Contrast (Ch%d)", channel), ylab="uV", xlab="s")
  lines(t_vec, x_c0, col="blue")
  legend("topright", c("Class 1 (Task High)", "Class 0 (Task Low)"), col=c("red", "blue"), lty=1, cex=0.8)
  
  # Plot 2: PSD (Check for 20Hz Task Peak + Alpha Background)
  spec1 <- spectrum(x_c1, plot=FALSE)
  freq1 <- spec1$freq * fs
  psd1 <- 10 * log10(pmax(spec1$spec, .Machine$double.eps))
  fmax_plot <- min(max(freq1), fmax)
  
  plot(freq1, psd1, type="l", xlim=c(0, fmax_plot),
       main="PSD (Class 1)", xlab="Hz", ylab="dB")
  if(!is.null(target_freqs)) {
    abline(v = target_freqs, col=rgb(0,0,0,0.2), lty=3)
  }
  abline(v = 20, col="red", lty=2) # Task Freq
  legend("topright", c("Task (20Hz)", "Bg Freqs"), col=c("red", "black"), lty=c(2,3), bty="n", cex=0.8)
  
  # Plot 3: Spatial Covariance (1Hz HP)
  X_hp_mat <- apply(X_tensor[,,idx1], 2, function(x) signal::filtfilt(bf_val, x))
  cov_mat <- cov(X_hp_mat)
  image(cov_mat[,ncol(cov_mat):1], main="Spatial Cov (1Hz HP)", axes=FALSE)
  
  # Plot 4: Audit
  audit <- as.data.frame(sim_res$audit)
  snr_values <- NULL
  if("Realized_SNR_Neural" %in% names(audit)) {
    snr_values <- as.numeric(audit$Realized_SNR_Neural)
  }
  
  if(!is.null(snr_values)) {
    if(is.null(snr_target_db)) snr_target_db <- stats::median(snr_values, na.rm = TRUE)
    y_lim <- range(c(snr_values, snr_target_db), na.rm = TRUE)
    if(!all(is.finite(y_lim)) || diff(y_lim) < 1e-6) y_lim <- c(0, 10)
    
    plot(snr_values, type='b', ylim=y_lim,
         main="SNR Audit", ylab="dB", xlab="Trial")
    abline(h=snr_target_db, col="green", lty=2)
  } else {
    plot.new()
    text(0.5, 0.5, "No audit SNR data found")
  }
  
  par(mfrow=c(1,1))
  
  spec0 <- spectrum(x_c0, plot=FALSE)
  freq0 <- spec0$freq * fs
  dom_f0 <- freq0[which.max(spec0$spec)]
  dom_f1 <- freq1[which.max(spec1$spec)]
  
  metrics <- list(
    class0_trial = idx0,
    class1_trial = idx1,
    channel = channel,
    dominant_freq_class0_hz = dom_f0,
    dominant_freq_class1_hz = dom_f1,
    snr_mean_db = if(!is.null(snr_values)) mean(snr_values, na.rm = TRUE) else NA_real_,
    snr_sd_db = if(!is.null(snr_values)) stats::sd(snr_values, na.rm = TRUE) else NA_real_,
    snr_target_db = snr_target_db
  )
  
  if(return_metrics) return(metrics)
  invisible(metrics)
}
