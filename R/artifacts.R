#' Generate Synthetic Physiological Artifacts (EOG, EMG, Drift)
#'
#' @description
#' Generates a mixture of common EEG artifacts to simulate realistic noise contamination.
#' This includes ocular artifacts (blinks), muscle artifacts (EMG bursts), and low-frequency
#' electrode drifts.
#'
#' @details
#' The function simulates three distinct types of noise:
#' \enumerate{
#'   \item \strong{EOG (Eye Blinks):} Modeled as Bell-shaped cosine pulses (Hanning-like) roughly 300ms in duration.
#'     They occur according to a Poisson process (~15/min) and are projected spatially from a virtual
#'     frontal source (approximate Fpz location) to the sensors based on distance.
#'   \item \strong{EMG (Muscle Activity):} Modeled as short bursts (~100ms) of high-frequency noise (>20Hz).
#'     These are spatially localized to random clusters of 3 neighboring electrodes and occur
#'     sporadically (Poisson process, ~5/min).
#'   \item \strong{Drift:} Modeled as a non-stationary random walk (integrated Gaussian noise)
#'     superimposed with a linear trend. This simulates electrode impedance shifts and DC offsets.
#'     The drift is normalized to a maximum amplitude of 100 units.
#' }
#'
#' @param n_time Integer. Number of time points to generate.
#' @param n_channels Integer. Number of EEG channels.
#' @param fs Numeric. Sampling frequency in Hz.
#' @param coords_sens Matrix (n_channels x 3). XYZ coordinates of the sensors, used to calculate
#'   spatial projections for EOG and spatial clustering for EMG.
#'
#' @return A list containing two matrices (n_time x n_channels):
#' \describe{
#'   \item{\code{N_fast}}{Contains high-frequency artifacts (EOG + EMG).}
#'   \item{\code{N_drift}}{Contains low-frequency drift artifacts.}
#' }
#'
#' @importFrom stats rpois rnorm dist
#' @importFrom signal butter filtfilt
#' @export
sim_artifacts <- function(n_time, n_channels, fs, coords_sens) {
  
  N_fast <- matrix(0, n_time, n_channels)
  
  # 1. EOG
  fpz <- c(0.9, 0, 0.4) 
  dists <- sqrt(rowSums((coords_sens - matrix(fpz, nrow=n_channels, ncol=3, byrow=T))^2))
  proj_eog <- exp(-dists^2 / 0.1)
  
  n_blink <- round(0.3 * fs) 
  t_blink <- seq(-pi, pi, length.out=n_blink)
  blink_shape <- (1 + cos(t_blink))/2 
  
  n_events_eog <- rpois(1, (n_time/fs) * (15/60)) 
  if(n_events_eog > 0) {
    possible_starts <- 1:(n_time - n_blink)
    if(length(possible_starts) > 0) {
      onsets <- sample(possible_starts, min(n_events_eog, length(possible_starts)))
      ts_eog <- numeric(n_time)
      for(t_start in onsets) {
        ts_eog[t_start:(t_start+n_blink-1)] <- ts_eog[t_start:(t_start+n_blink-1)] + blink_shape
      }
      N_fast <- N_fast + (matrix(ts_eog, ncol=1) %*% matrix(proj_eog, nrow=1)) * 50
    }
  }
  
  # 2. EMG
  dist_SS <- as.matrix(dist(coords_sens))
  n_bursts <- rpois(1, (n_time/fs) * (5/60)) 
  
  if(n_bursts > 0) {
    bf_emg <- signal::butter(4, 20/(fs/2), type="high")
    
    for(i in 1:n_bursts) {
      center_ch <- sample(1:n_channels, 1)
      cluster <- order(dist_SS[center_ch,])[1:3]
      dur_samps <- round(0.1 * fs)
      start_t <- sample(1:(n_time - dur_samps), 1)
      
      buffer_len <- fs 
      noise_buffer <- matrix(rnorm(buffer_len * length(cluster)), buffer_len, length(cluster))
      noise_hp_long <- apply(noise_buffer, 2, function(x) signal::filtfilt(bf_emg, x))
      crop_start <- floor((buffer_len - dur_samps)/2)
      noise_burst <- noise_hp_long[crop_start:(crop_start+dur_samps-1), ]
      
      N_fast[start_t:(start_t+dur_samps-1), cluster] <- N_fast[start_t:(start_t+dur_samps-1), cluster] + noise_burst * 20
    }
  }
  
  # B. Drift
  drift <- apply(matrix(rnorm(n_time*n_channels), n_time, n_channels), 2, cumsum)
  trend <- seq(0, 1, length.out = n_time)
  trend_mat <- matrix(trend, n_time, n_channels) * matrix(rnorm(n_channels, 0, 50), n_time, n_channels, byrow=T)
  drift <- drift + trend_mat
  
  max_exc <- max(abs(drift))
  if(max_exc > 0) drift <- drift / max_exc * 100
  
  return(list(N_fast = N_fast, N_drift = drift))
}