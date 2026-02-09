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
#'
#' @return NULL. The function is called for its side effect of generating plots.
#'
#' @importFrom stats spectrum cov
#' @importFrom graphics par plot lines legend abline image
#' @importFrom signal butter filtfilt
#' @export
validate_sim_eeg <- function(sim_res) {
  
  X_tensor <- sim_res$data
  fs <- sim_res$params$fs
  target_freqs <- sim_res$params$target_freqs
  
  # Validation Filter: 1Hz HPF (Standard EEG Preproc)
  bf_val <- signal::butter(4, 1/(fs/2), type="high")
  
  # Plot 1: Temporal (Class Comparison)
  # Average Trial for Class 0 vs Class 1
  cls <- sim_res$labels
  idx0 <- which(cls == 0)[1]; idx1 <- which(cls == 1)[1]
  
  x_c0 <- signal::filtfilt(bf_val, X_tensor[, 1, idx0])
  x_c1 <- signal::filtfilt(bf_val, X_tensor[, 1, idx1])
  t_vec <- (1:length(x_c0)) / fs
  
  par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
  
  plot(t_vec, x_c1, type="l", col="red", main="Class Contrast (Ch1)", ylab="uV", xlab="s")
  lines(t_vec, x_c0, col="blue")
  legend("topright", c("Class 1 (Task High)", "Class 0 (Task Low)"), col=c("red", "blue"), lty=1, cex=0.8)
  
  # Plot 2: PSD (Check for 20Hz Task Peak + Alpha Background)
  spec <- spectrum(x_c1, plot=FALSE)
  plot(spec$freq * fs, 10*log10(spec$spec), type="l", xlim=c(0, 60),
       main="PSD (Class 1)", xlab="Hz", ylab="dB")
  abline(v = target_freqs, col=rgb(0,0,0,0.2), lty=3)
  abline(v = 20, col="red", lty=2) # Task Freq
  legend("topright", c("Task (20Hz)", "Bg Freqs"), col=c("red", "black"), lty=c(2,3), bty="n", cex=0.8)
  
  # Plot 3: Spatial Covariance (1Hz HP)
  X_hp_mat <- apply(X_tensor[,,1], 2, function(x) signal::filtfilt(bf_val, x))
  cov_mat <- cov(X_hp_mat)
  image(cov_mat[,ncol(cov_mat):1], main="Spatial Cov (1Hz HP)", axes=FALSE)
  
  # Plot 4: Audit
  audit <- as.data.frame(sim_res$audit)
  plot(audit$Realized_SNR_Neural, type='b', ylim=c(0, 10), 
       main="SNR Audit", ylab="dB", xlab="Trial")
  abline(h=5, col="green", lty=2)
  
  par(mfrow=c(1,1))
}