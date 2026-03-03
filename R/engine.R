#' Master Controller for Physics-Constrained Tensor EEG Simulation
#'
#' @description
#' Orchestrates the entire simulation pipeline to generate a single session of EEG data.
#' This function integrates:
#' \enumerate{
#'   \item \strong{Spatial Geometry:} Generates sensor positions and a Laplacian-smoothed mixing matrix.
#'   \item \strong{Manifold Dynamics:} Computes trial-specific rotation matrices to simulate non-stationarity.
#'   \item \strong{Source Activity:} Generates task-related ERPs (Gabor) and background oscillations (VAR2).
#'   \item \strong{Artifact Injection:} Adds EOG, EMG, and Drift noise with precise SNR control.
#'   \item \strong{Closed-Loop SNR:} Automatically adjusts component amplitudes to meet target dB levels.
#' }
#'
#' @details
#' The simulation implements a class-dependent logic:
#' \itemize{
#'   \item \strong{Class 1:} Characterized by stronger task activation (amplitude 1.5x) and
#'     weaker background activity (amplitude 0.8x), simulating Event-Related Desynchronization (ERD) or strong ERPs.
#'   \item \strong{Class 0:} Baseline task activation (amplitude 1.0x) and normal background (1.0x).
#' }
#' The Signal-to-Noise Ratio (SNR) is calibrated using "Effective AC Power" (High-Pass > 0.1Hz)
#' to ensure that high-amplitude low-frequency drifts do not bias the mixing ratios.
#'
#' @param n_trials Integer. Total number of trials to simulate.
#' @param n_time Integer. Number of time points per trial.
#' @param n_channels Integer. Number of EEG channels.
#' @param n_sources Integer. Number of latent sources.
#' @param fs Numeric. Sampling frequency in Hz.
#' @param snr_neural_db Numeric. Target SNR (dB) for the Neural signal relative to the Background brain activity.
#' @param snr_artifact_db Numeric. Target SNR (dB) for the Total Neural signal relative to Fast Artifacts (EOG/EMG).
#' @param drift_power_ratio Numeric. The ratio of Drift power relative to Neural power.
#'   A value of 0.5 means drift power is half the neural power (in the AC band).
#' @param target_freqs Numeric vector (optional). Specific peak frequencies for the VAR(2) background sources.
#'   If NULL, frequencies are randomly drawn from 8-20 Hz.
#' @param class_labels Integer vector (optional). A vector of length \code{n_trials} containing 0s and 1s.
#'   If NULL, defaults to alternating 0/1.
#' @param seed Integer (optional). If provided, sets RNG seed for reproducible simulation.
#' @param verbose Logical. If TRUE, prints simulation progress messages.
#'
#' @return A list containing the simulation results:
#' \describe{
#'   \item{\code{data}}{Numeric Array (n_time x n_channels x n_trials). The final simulated EEG tensor.}
#'   \item{\code{geometry}}{List. Contains sensor coordinates and the base mixing matrix.}
#'   \item{\code{audit}}{Data frame. Trial-by-trial log of class labels and realized SNR values.}
#'   \item{\code{labels}}{Integer vector. The class labels used for generation.}
#'   \item{\code{params}}{List. Metadata including sampling rate and source frequencies.}
#' }
#'
#' @importFrom stats rnorm rlnorm
#' @export
sim_eeg_master <- function(n_trials = 20, 
                           n_time = 500, 
                           n_channels = 64, 
                           n_sources = 10,
                           fs = 250,
                           snr_neural_db = 5,
                           snr_artifact_db = 0,
                           drift_power_ratio = 0.5, 
                           target_freqs = NULL,
                           class_labels = NULL,
                           seed = NULL,
                           verbose = TRUE) {
  is_whole_number <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      abs(x - round(x)) < .Machine$double.eps^0.5
  }
  if(!is_whole_number(n_trials) || n_trials < 1) {
    stop("n_trials must be a positive integer.")
  }
  if(!is_whole_number(n_time) || n_time < 1) {
    stop("n_time must be a positive integer.")
  }
  if(!is_whole_number(n_channels) || n_channels < 1) {
    stop("n_channels must be a positive integer.")
  }
  if(!is_whole_number(n_sources) || n_sources < 1) {
    stop("n_sources must be a positive integer.")
  }
  if(!is.numeric(fs) || length(fs) != 1L || !is.finite(fs) || fs <= 0) {
    stop("fs must be a single positive finite number.")
  }
  if(!is.numeric(snr_neural_db) || length(snr_neural_db) != 1L || !is.finite(snr_neural_db)) {
    stop("snr_neural_db must be a single finite number.")
  }
  if(!is.numeric(snr_artifact_db) || length(snr_artifact_db) != 1L || !is.finite(snr_artifact_db)) {
    stop("snr_artifact_db must be a single finite number.")
  }
  if(!is.numeric(drift_power_ratio) || length(drift_power_ratio) != 1L || !is.finite(drift_power_ratio) || drift_power_ratio < 0) {
    stop("drift_power_ratio must be a single non-negative finite number.")
  }
  if(!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
  
  n_trials <- as.integer(round(n_trials))
  n_time <- as.integer(round(n_time))
  n_channels <- as.integer(round(n_channels))
  n_sources <- as.integer(round(n_sources))
  
  if(!is.null(seed)) {
    if(!is_whole_number(seed)) {
      stop("seed must be a single finite integer.")
    }
    set.seed(as.integer(round(seed)))
  }
  
  if(!is.numeric(target_freqs) && !is.null(target_freqs)) {
    stop("target_freqs must be NULL or a numeric vector of length n_sources.")
  }
  nyquist <- fs / 2
  if(is.null(target_freqs)) {
    lower <- min(8, nyquist * 0.25)
    upper <- min(20, nyquist * 0.9)
    if(lower <= 0 || upper <= 0 || lower >= upper) {
      stop("fs is too low for default target_freqs; please provide target_freqs in (0, fs/2).")
    }
    target_freqs <- stats::runif(n_sources, lower, upper)
  } else {
    if(length(target_freqs) != n_sources) {
      stop("target_freqs must have length n_sources.")
    }
    target_freqs <- as.numeric(target_freqs)
    if(any(!is.finite(target_freqs)) || any(target_freqs <= 0) || any(target_freqs >= nyquist)) {
      stop("target_freqs must contain finite values in (0, fs/2).")
    }
  }
  
  if(is.null(class_labels)) {
    class_labels <- rep(c(0L, 1L), length.out = n_trials)
  } else {
    if(!(is.numeric(class_labels) || is.integer(class_labels))) {
      stop("class_labels must be a numeric/integer vector containing only 0 and 1.")
    }
    if(length(class_labels) != n_trials) {
      stop("class_labels must have length n_trials.")
    }
    if(any(!is.finite(class_labels)) || any(abs(class_labels - round(class_labels)) >= .Machine$double.eps^0.5)) {
      stop("class_labels must only contain integer 0 and 1 values.")
    }
    class_labels <- as.integer(class_labels)
    if(any(is.na(class_labels)) || !all(class_labels %in% c(0L, 1L))) {
      stop("class_labels must only contain 0 and 1.")
    }
  }
  
  geo <- generate_geometry_mixing(n_channels, n_sources)
  A_base <- geo$A_base
  coords <- geo$coords_sens
  rotations <- generate_drift_rotations(n_sources, n_trials)
  var_system <- setup_var2_system(n_sources, fs, target_freqs)
  
  X_Tensor <- array(0, dim = c(n_time, n_channels, n_trials))
  audit_list <- vector("list", n_trials)
  
  if(verbose) {
    message(sprintf("Simulating %d trials...", n_trials))
  }
  for(k in seq_len(n_trials)) {
    if(verbose && (k %% 10L == 0L || k == n_trials)) {
      message(sprintf("  trial %d/%d", k, n_trials))
    }
    
    A_k <- A_base %*% rotations[[k]]
    
    # Source Generation with Class Logic (Fix B2)
    # Class 1: Stronger Task, Weaker Background (ERD)
    cls <- class_labels[k]
    task_amp <- if(cls == 1) 1.5 else 1.0
    bg_amp   <- if(cls == 1) 0.8 else 1.0
    
    tau <- rnorm(1, 0, 20) 
    gamma <- rlnorm(1, 0, 0.1) 
    S_task <- sim_source_task(n_time, n_sources, fs, tau_ms=tau, gamma=gamma) * task_amp
    S_bg <- sim_source_var2(n_time, n_sources, var_system) * bg_amp
    
    X_task_pure <- S_task %*% t(A_k)
    X_bg_pure   <- S_bg %*% t(A_k)
    
    # Neural SNR
    P_task <- calc_ac_power(X_task_pure, fs)
    P_bg   <- calc_ac_power(X_bg_pure, fs)
    
    if(P_bg > 1e-9) { 
      g_bg <- sqrt(P_task / (P_bg * 10^(snr_neural_db / 10)))
    } else { g_bg <- 0 }
    
    X_neural <- X_task_pure + g_bg * X_bg_pure
    P_neural <- calc_ac_power(X_neural, fs)
    
    # Artifacts
    arts <- sim_artifacts(n_time, n_channels, fs, coords)
    N_fast  <- arts$N_fast
    N_drift_raw <- arts$N_drift
    
    # Fast Artifact SNR
    P_fast <- calc_ac_power(N_fast, fs)
    if(P_fast > 1e-9) {
      g_art <- sqrt(P_neural / (P_fast * 10^(snr_artifact_db / 10)))
    } else { g_art <- 0 }
    
    # Drift Budget Control (Fix B3: Using 0.1Hz HP Power)
    P_drift_ac <- calc_ac_power(N_drift_raw, fs)
    if(P_drift_ac > 1e-9) {
      g_drift <- sqrt((P_neural * drift_power_ratio) / P_drift_ac)
    } else { g_drift <- 0 }
    
    # Final Assemble
    X_final <- X_neural + g_art * N_fast + g_drift * N_drift_raw
    X_final[is.na(X_final)] <- 0
    X_Tensor[,,k] <- X_final
    
    snr_denom <- g_bg^2 * P_bg
    realized_snr <- if(P_task > 1e-12 && snr_denom > 1e-12) {
      10 * log10(P_task / snr_denom)
    } else {
      NA_real_
    }
    
    audit_list[[k]] <- data.frame(
      trial = k,
      class = cls,
      Realized_SNR_Neural = realized_snr
    )
  }
  
  audit_df <- do.call(rbind, audit_list)
  
  return(list(
    data = X_Tensor,
    geometry = geo,
    audit = audit_df,
    labels = class_labels,
    params = list(fs = fs, target_freqs = target_freqs, seed = seed)
  ))
}
