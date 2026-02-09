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
  G <- matrix(rnorm(n_sources^2), n_sources, n_sources)
  Omega_base <- (G - t(G)) / 2
  frob_norm <- sqrt(sum(Omega_base^2))
  Omega_base <- Omega_base / frob_norm
  
  theta <- numeric(n_trials)
  theta[1] <- 0
  noise <- rnorm(n_trials, 0, sigma_eps)
  
  for(k in 2:n_trials) {
    theta[k] <- alpha_ou * theta[k-1] + noise[k]
  }
  
  R_list <- list()
  for(k in 1:n_trials) {
    R_list[[k]] <- expm::expm(theta[k] * Omega_base)
  }
  
  return(R_list)
}