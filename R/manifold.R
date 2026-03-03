#' Generate Source Rotation Matrix (Manifold Drift)
#'
#' @description
#' Simulates the non-stationary evolution of source orientations across trials.
#' Instead of adding additive noise to the mixing matrix, this function models drift
#' as a rotation on the source manifold, preserving the orthogonality and geometry
#' of the sources while allowing their projection to change over time.
#'
#' @details
#' The drift is modeled as a geodesic random walk on the Lie Group of rotation matrices \eqn{SO(n)}.
#' \enumerate{
#'   \item \strong{Generator:} A random skew-symmetric matrix \eqn{\mathbf{\Omega}_{base}} is generated
#'     to define a fixed axis of rotation in the high-dimensional source space.
#'   \item \strong{Dynamics:} A scalar angle parameter \eqn{\theta_k} evolves according to a
#'     mean-reverting Ornstein-Uhlenbeck (AR(1)) process:
#'     \deqn{\theta_k = \alpha \theta_{k-1} + \epsilon, \quad \epsilon \sim \mathcal{N}(0, \sigma^2)}
#'   \item \strong{Mapping:} The rotation matrix for trial \eqn{k} is obtained via the Matrix Exponential map:
#'     \deqn{\mathbf{R}_k = \exp(\theta_k \mathbf{\Omega}_{base})}
#' }
#' This ensures that \eqn{\mathbf{R}_k} is always an orthogonal matrix, preventing the
#' source mixing from degenerating (i.e., preserving rank).
#'
#' @param n_sources Integer. The number of sources (dimension of the rotation matrix).
#' @param n_trials Integer. The number of trials to simulate.
#' @param alpha_ou Numeric. The autoregressive coefficient for the OU process (default: 0.95).
#'   Values closer to 1 imply strong memory (slow drift); values closer to 0 imply white noise (jitter).
#' @param sigma_eps Numeric. The standard deviation of the innovation noise (default: 0.05).
#'   Controls the magnitude/speed of the drift.
#'
#' @return A list of length \code{n_trials}, where each element is an
#'   \code{n_sources x n_sources} numeric rotation matrix.
#'
#' @importFrom expm expm
#' @importFrom stats rnorm
#' @export
generate_drift_rotations <- function(n_sources, n_trials, 
                                     alpha_ou = 0.95, sigma_eps = 0.05) {
  is_whole_number <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      abs(x - round(x)) < .Machine$double.eps^0.5
  }
  if(!is_whole_number(n_sources) || n_sources < 1) {
    stop("n_sources must be a positive integer.")
  }
  if(!is_whole_number(n_trials) || n_trials < 1) {
    stop("n_trials must be a positive integer.")
  }
  if(!is.numeric(alpha_ou) || length(alpha_ou) != 1L || !is.finite(alpha_ou) || abs(alpha_ou) >= 1) {
    stop("alpha_ou must be a single finite number in (-1, 1).")
  }
  if(!is.numeric(sigma_eps) || length(sigma_eps) != 1L || !is.finite(sigma_eps) || sigma_eps < 0) {
    stop("sigma_eps must be a single non-negative finite number.")
  }
  
  n_sources <- as.integer(round(n_sources))
  n_trials <- as.integer(round(n_trials))
  
  if(n_sources == 1L) {
    return(rep(list(matrix(1, 1, 1)), n_trials))
  }
  
  G <- matrix(stats::rnorm(n_sources^2), n_sources, n_sources)
  Omega_base <- (G - t(G)) / 2
  frob_norm <- sqrt(sum(Omega_base^2))
  if(frob_norm <= .Machine$double.eps) {
    return(rep(list(diag(n_sources)), n_trials))
  }
  Omega_base <- Omega_base / frob_norm
  
  theta <- numeric(n_trials)
  theta[1] <- 0
  noise <- stats::rnorm(n_trials, 0, sigma_eps)
  
  for(k in 2:n_trials) {
    theta[k] <- alpha_ou * theta[k-1] + noise[k]
  }
  
  R_list <- vector("list", n_trials)
  for(k in seq_len(n_trials)) {
    R_list[[k]] <- expm::expm(theta[k] * Omega_base)
  }
  
  return(R_list)
}
