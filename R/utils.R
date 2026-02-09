#' Compute Radial Basis Function (RBF) Kernel
#'
#' @description
#' Calculates the pair-wise Gaussian RBF kernel matrix between two sets of data points.
#' This is used in the spatial module to determine the raw mixing weights based on
#' the Euclidean distance between sensors and sources.
#'
#' @details
#' The kernel is defined as:
#' \deqn{K(x, y) = \exp(-\gamma ||x - y||^2)}
#' where \eqn{\gamma} is determined by the \code{sigma} and \code{standard_scale} parameters.
#'
#' @param x Numeric matrix. First set of coordinates (e.g., Sensor positions).
#'   Rows represent points, columns represent dimensions (x, y, z).
#' @param y Numeric matrix. Second set of coordinates (e.g., Source positions).
#'   Rows represent points, columns represent dimensions.
#' @param sigma Numeric. The bandwidth parameter (standard deviation) of the kernel.
#'   Controls the smoothness/decay of the weights.
#' @param standard_scale Logical. Controls the scaling of the exponent.
#'   If \code{TRUE} (default), \eqn{\gamma = 1 / (2\sigma^2)}.
#'   If \code{FALSE}, \eqn{\gamma = 1 / \sigma^2}.
#'
#' @return A numeric matrix of dimensions \code{nrow(x)} by \code{nrow(y)},
#'   containing the computed kernel weights.
#'
#' @export
rbf_kernel <- function(x, y, sigma, standard_scale = TRUE) {
  x <- as.matrix(x); y <- as.matrix(y)
  stopifnot(ncol(x) == ncol(y), is.numeric(sigma), length(sigma) == 1L, sigma > 0)
  
  xx <- rowSums(x * x)
  yy <- rowSums(y * y)
  
  D2 <- outer(xx, yy, "+") - 2 * tcrossprod(x, y) 
  D2 <- pmax(D2, 0) 
  
  gamma <- if (standard_scale) 1 / (2 * sigma^2) else 1 / (sigma^2)
  exp(-gamma * D2)
}


#' Calculate Effective AC Power (High-Pass Filtered)
#'
#' @description
#' Calculates the mean squared power of a signal matrix after removing low-frequency
#' drifts and DC offsets. This metric is critical for the Closed-Loop SNR calibration
#' system, ensuring that high-amplitude drifts do not skew the signal-to-noise ratio calculations.
#'
#' @details
#' The function applies a 4th-order Butterworth high-pass filter with a cutoff
#' frequency of 0.1 Hz using forward-backward filtering (\code{filtfilt}) to prevent
#' phase distortion.
#'
#' @param X Numeric matrix (Time x Channels). The input EEG data trial.
#' @param fs Numeric. The sampling frequency in Hz.
#'
#' @return A single numeric value representing the average AC power across all
#'   channels and time points. Returns 0 if input contains NAs or Infinite values.
#'
#' @importFrom signal butter filtfilt
#' @export
calc_ac_power <- function(X, fs) {
  if(any(is.na(X)) || any(is.infinite(X))) return(0)
  
  # Fix B1: Use 0.1 Hz to capture drift AC components better
  bf <- signal::butter(4, 0.1 / (fs/2), type = "high")
  X_filt <- apply(X, 2, function(x) signal::filtfilt(bf, x))
  
  P_ac <- sum(X_filt^2) / length(X_filt)
  
  if(is.na(P_ac)) return(0)
  return(P_ac)
}