#' Setup Structured VAR(2) System
#'
#' @description
#' Constructs the coefficient matrices for a Vector Autoregressive model of order 2 (VAR(2)).
#' This system generates the background brain oscillations (e.g., Alpha, Beta rhythms).
#'
#' @details
#' The function designs a stable VAR(2) process where:
#' \enumerate{
#'   \item \strong{Diagonal Elements:} Are tuned to generate specific oscillation frequencies
#'     (\code{target_freqs}) by placing poles on the complex plane with magnitude \code{r_pole} (0.95).
#'   \item \strong{Off-Diagonal Elements:} Represent functional connectivity (coupling) between sources.
#'     These are initialized randomly but scaled down iteratively (\code{gamma}) until the
#'     spectral radius of the companion matrix is < 0.995, ensuring system stability.
#' }
#' The process equation is: \eqn{\mathbf{s}(t) = \mathbf{\Phi}_1 \mathbf{s}(t-1) + \mathbf{\Phi}_2 \mathbf{s}(t-2) + \epsilon(t)}
#'
#' @param n_sources Integer. Number of sources in the VAR model.
#' @param fs Numeric. Sampling frequency in Hz.
#' @param target_freqs Numeric vector (length = n_sources). The peak frequency for each source.
#'
#' @return A list containing the VAR(2) coefficient matrices:
#' \item{Phi1}{Matrix (n_sources x n_sources). Lag-1 coefficients.}
#' \item{Phi2}{Matrix (n_sources x n_sources). Lag-2 coefficients.}
#'
#' @importFrom stats rbinom rnorm
#' @export
setup_var2_system <- function(n_sources, fs, target_freqs) {
  if(!.is_whole_number(n_sources) || n_sources < 1) {
    stop("n_sources must be a positive integer.")
  }
  if(!is.numeric(fs) || length(fs) != 1L || !is.finite(fs) || fs <= 0) {
    stop("fs must be a single positive finite number.")
  }
  n_sources <- as.integer(round(n_sources))
  
  if(!is.numeric(target_freqs) || length(target_freqs) != n_sources || any(!is.finite(target_freqs))) {
    stop("target_freqs must be a finite numeric vector of length n_sources.")
  }
  target_freqs <- as.numeric(target_freqs)
  nyquist <- fs / 2
  if(any(target_freqs <= 0) || any(target_freqs >= nyquist)) {
    stop(sprintf("target_freqs values must be in (0, fs/2 = %.3f).", nyquist))
  }
  
  Phi1 <- matrix(0, n_sources, n_sources)
  Phi2 <- matrix(0, n_sources, n_sources)
  r_pole <- 0.95 
  
  for(i in seq_len(n_sources)) {
    omega <- 2 * pi * target_freqs[i] / fs
    Phi1[i, i] <- 2 * r_pole * cos(omega)
    Phi2[i, i] <- -r_pole^2
  }
  
  coupling_mask <- matrix(stats::rbinom(n_sources^2, 1, 0.2), n_sources, n_sources)
  diag(coupling_mask) <- 0
  couplings_base <- matrix(stats::rnorm(n_sources^2, 0, 0.05), n_sources, n_sources) * coupling_mask
  
  is_stable <- FALSE
  gamma <- 1.0 
  
  while(!is_stable && gamma > 0.01) {
    Phi1_curr <- diag(diag(Phi1)) + gamma * couplings_base
    
    top <- cbind(Phi1_curr, Phi2)
    bot <- cbind(diag(n_sources), matrix(0, n_sources, n_sources))
    Comp <- rbind(top, bot)
    
    rho <- max(Mod(eigen(Comp, only.values = TRUE)$values))
    
    if(rho < 0.995) {
      is_stable <- TRUE
      Phi1 <- Phi1_curr 
    } else {
      gamma <- gamma * 0.9 
    }
  }
  
  if(!is_stable) {
    Phi1 <- diag(diag(Phi1), n_sources, n_sources)
  }
  
  return(list(Phi1 = Phi1, Phi2 = Phi2))
}

#' Generate VAR(2) Time Series
#'
#' @description
#' Simulates the time series data for background neural activity using the pre-calculated
#' VAR(2) coefficients.
#'
#' @details
#' Iterates the difference equation:
#' \deqn{s(t) = \Phi_1 s(t-1) + \Phi_2 s(t-2) + \mathcal{N}(0, 1)}
#' Includes a burn-in period of 500 samples to remove transient effects and ensures
#' the output is zero-centered (mean removal).
#'
#' @param n_time Integer. Number of time points to generate.
#' @param n_sources Integer. Number of sources.
#' @param var_params List. The output from \code{\link{setup_var2_system}}, containing \code{Phi1} and \code{Phi2}.
#'
#' @return Numeric matrix (n_time x n_sources). The generated source time series.
#'
#' @importFrom stats rnorm
#' @export
sim_source_var2 <- function(n_time, n_sources, var_params) {
  if(!.is_whole_number(n_time) || n_time < 1) {
    stop("n_time must be a positive integer.")
  }
  if(!.is_whole_number(n_sources) || n_sources < 1) {
    stop("n_sources must be a positive integer.")
  }
  if(!is.list(var_params) || is.null(var_params$Phi1) || is.null(var_params$Phi2)) {
    stop("var_params must be a list containing Phi1 and Phi2.")
  }
  n_time <- as.integer(round(n_time))
  n_sources <- as.integer(round(n_sources))
  
  Phi1 <- as.matrix(var_params$Phi1)
  Phi2 <- as.matrix(var_params$Phi2)
  if(!all(dim(Phi1) == c(n_sources, n_sources)) || !all(dim(Phi2) == c(n_sources, n_sources))) {
    stop("Phi1 and Phi2 must both be matrices of size n_sources x n_sources.")
  }
  if(any(!is.finite(Phi1)) || any(!is.finite(Phi2))) {
    stop("Phi1 and Phi2 must contain only finite values.")
  }
  
  burn_in <- 500
  noise <- matrix(stats::rnorm((n_time + burn_in) * n_sources), n_time + burn_in, n_sources)
  S_temp <- matrix(0, n_time + burn_in, n_sources)
  
  for(t in 3:(n_time + burn_in)) {
    val <- as.numeric(Phi1 %*% S_temp[t-1, ] + Phi2 %*% S_temp[t-2, ] + noise[t, ])
    if(any(is.infinite(val))) val[is.infinite(val)] <- sign(val[is.infinite(val)]) * 1e10
    S_temp[t, ] <- val
  }
  
  result <- S_temp[(burn_in + 1):(n_time + burn_in), , drop = FALSE]
  result <- sweep(result, 2, colMeans(result), "-")
  return(result)
}

#' Generate Task Source (Gabor Wavelet)
#'
#' @description
#' Generates the task-related Event-Related Potential (ERP) source activity using a
#' time-warped Gabor wavelet.
#'
#' @details
#' The task signal is modeled as a 20Hz Gabor wavelet (Gaussian-windowed cosine).
#' It incorporates trial-specific variability via:
#' \itemize{
#'   \item \strong{Latency Jitter (\code{tau_ms}):} Shifts the peak time.
#'   \item \strong{Time Warping (\code{gamma}):} Stretches or compresses the waveform duration.
#' }
#' This simulates the P300/Beta-rebound components often seen in BCI paradigms.
#'
#' @param n_time Integer. Total number of time points.
#' @param n_sources Integer. Total number of sources.
#' @param fs Numeric. Sampling frequency in Hz.
#' @param tau_ms Numeric. Latency shift in milliseconds (default: 0).
#' @param gamma Numeric. Time scaling factor (default: 1.0).
#' @param active_idx Integer vector. Indices of the sources that are "task-active" (default: 1:3).
#'   Other sources will remain silent (zeros).
#'
#' @return Numeric matrix (n_time x n_sources). The task-related source activity.
#'
#' @export
sim_source_task <- function(n_time, n_sources, fs,
                            tau_ms = 0, gamma = 1, active_idx = 1:3) {
  if(!.is_whole_number(n_time) || n_time < 1) {
    stop("n_time must be a positive integer.")
  }
  if(!.is_whole_number(n_sources) || n_sources < 1) {
    stop("n_sources must be a positive integer.")
  }
  if(!is.numeric(fs) || length(fs) != 1L || !is.finite(fs) || fs <= 0) {
    stop("fs must be a single positive finite number.")
  }
  if(!is.numeric(tau_ms) || length(tau_ms) != 1L || !is.finite(tau_ms)) {
    stop("tau_ms must be a single finite number.")
  }
  if(!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma <= 0) {
    stop("gamma must be a single positive finite number.")
  }
  if(!is.null(active_idx) && !(is.numeric(active_idx) || is.integer(active_idx))) {
    stop("active_idx must be a numeric or integer vector of source indices.")
  }
  
  n_time <- as.integer(round(n_time))
  n_sources <- as.integer(round(n_sources))
  
  T_duration <- n_time / fs
  t_vec <- (seq_len(n_time) - 1) / fs
  center_time <- T_duration / 2
  S <- matrix(0, n_time, n_sources)
  tau_s <- tau_ms / 1000 
  if(is.null(active_idx)) {
    active_idx <- integer(0)
  }
  is_valid_idx <- is.finite(active_idx) &
    abs(active_idx - round(active_idx)) < .Machine$double.eps^0.5
  active_idx <- as.integer(active_idx[is_valid_idx])
  active_idx <- unique(active_idx[is.finite(active_idx) & active_idx >= 1L & active_idx <= n_sources])
  if(length(active_idx) == 0L) {
    return(S)
  }
  
  gabor <- function(t) {
    # Fix A1: Use 20 Hz carrier (Beta band) instead of 5Hz
    # Matches target freqs better for PSD validation
    exp(- (t^2) / (2 * 0.05^2)) * cos(2 * pi * 20 * t)
  }
  
  for(i in active_idx) {
    t_prime <- (t_vec - center_time - tau_s) / gamma
    S[, i] <- gabor(t_prime)
  }
  return(S) 
}
