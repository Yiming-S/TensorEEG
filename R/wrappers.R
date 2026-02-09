#' Simulate a Full Session with Multiple Runs (Cross-Run Drift)
#'
#' @description
#' Simulates a complete EEG recording session consisting of multiple experimental "runs"
#' separated by rest periods. This is critical for validating algorithms that must handle
#' non-stationarity (drift) across different blocks of time.
#'
#' @details
#' \strong{Drift Injection Strategy:}
#' To simulate the natural evolution of brain signals and electrode impedance over a long session,
#' this function generates a single continuous "Giant Simulation" internally. It then:
#' \enumerate{
#'   \item Extracts \code{trials_per_run} trials for the active run.
#'   \item Discards the next \code{gap_trials} trials to simulate the "break" or "rest" period.
#'     During this gap, the underlying manifold drift (rotation matrix) continues to evolve,
#'     ensuring that the next run starts with a slightly different distribution than the previous one.
#' }
#'
#' @param n_runs Integer. Number of experimental runs (blocks) to simulate (default: 3).
#' @param trials_per_run Integer. Number of trials per run (default: 20).
#' @param gap_trials Integer. Number of "virtual" trials to simulate but discard between runs (default: 15).
#'   Higher values create larger distributional shifts between runs.
#' @param ... Additional arguments passed to \code{\link{sim_eeg_master}} (e.g., \code{n_time}, \code{n_sources}).
#'
#' @return A list formatted similarly to Physionet datasets, containing:
#' \describe{
#'   \item{\code{x}}{List of matrices. The EEG data for all valid trials.}
#'   \item{\code{y}}{Factor. Class labels for all valid trials.}
#'   \item{\code{run}}{Factor. Run ID for each trial (useful for Leave-One-Run-Out CV).}
#'   \item{\code{fs}}{Numeric. Sampling frequency.}
#'   \item{\code{ch_names}}{Character vector. Channel names.}
#' }
#'
#' @export
sim_multirun_session <- function(n_runs = 3, 
                                 trials_per_run = 20, 
                                 gap_trials = 15,
                                 ...) { # Pass other args (n_time, n_sources) to master
  
  # 1. Calculate Total Trials needed
  # Structure: [Run 1] --gap-- [Run 2] --gap-- [Run 3]
  total_active_trials <- n_runs * trials_per_run
  total_gap_trials    <- (n_runs - 1) * gap_trials
  grand_total         <- total_active_trials + total_gap_trials
  
  cat(sprintf("--- Simulating Multi-Run Session ---\n"))
  cat(sprintf("Config: %d Runs, %d Trials/Run, %d Gap Trials (Drift Injection)\n", 
              n_runs, trials_per_run, gap_trials))
  
  # 2. Run ONE Giant Simulation (Maintains Geometry & Drift Continuity)
  # We pass '...' to sim_eeg_master to allow user control (e.g., target_freqs)
  giant_sim <- sim_eeg_master(n_trials = grand_total, ...)
  
  # 3. Slice and Dice (Extract Runs, Discard Gaps)
  final_x_list <- list()
  final_y_vec  <- c()
  run_indicator <- c()
  
  current_idx <- 1
  
  for(r in 1:n_runs) {
    # Define start and end for this run
    end_idx <- current_idx + trials_per_run - 1
    
    cat(sprintf("Extracting Run %d: Indices [%d - %d]\n", r, current_idx, end_idx))
    
    # Extract data slices
    for(k in current_idx:end_idx) {
      # Extract matrix [Time x Channel]
      final_x_list[[length(final_x_list) + 1]] <- giant_sim$data[,,k]
    }
    
    # Extract labels
    final_y_vec <- c(final_y_vec, giant_sim$labels[current_idx:end_idx])
    
    # Record Run ID
    run_indicator <- c(run_indicator, rep(r, trials_per_run))
    
    # Move pointer past the Run AND the Gap
    current_idx <- end_idx + 1 + gap_trials
  }
  
  # 4. Construct Final Physionet-Style List
  formatted_data <- list(
    x = final_x_list,           # List of [500x64] matrices
    y = factor(final_y_vec),    # Class labels
    run = factor(run_indicator),# New: Run ID for cross-validation
    fs = giant_sim$params$fs,
    ch_names = paste0("Ch", 1:dim(giant_sim$data)[2])
  )
  
  return(formatted_data)
}

