#' Generate Virtual Head Geometry and Mixing Matrix
#'
#' @description
#' Generates a virtual EEG recording setup including 3D sensor coordinates (Fibonacci Grid),
#' random source locations, and a physically constrained mixing matrix. This function
#' serves as the static geometric foundation for the simulation.
#'
#' @details
#' This function simulates the physics of volume conduction without requiring a
#' full BEM/FEM head model. It employs a graph-theoretical approach:
#' \enumerate{
#'   \item \strong{Sensor Geometry:} A Fibonacci Lattice is used to generate evenly
#'     distributed sensor positions on the upper hemisphere (\eqn{z \ge 0}).
#'   \item \strong{Source Geometry:} Sources are randomly placed within a sphere
#'     of radius 0.8 to ensure they are located "deep" inside the brain volume.
#'   \item \strong{Physics Constraint:} A raw Radial Basis Function (RBF) lead field
#'     is smoothed using a Normalized Graph Laplacian to simulate skull conductivity.
#'     The smoothing operation is defined as Tikhonov regularization:
#'     \deqn{\mathbf{A}_{smooth} = (\mathbf{I} + \lambda \mathcal{L}_{sym})^{-1} \mathbf{A}_{raw}}
#' }
#' The resulting matrix is column-normalized to ensure unit energy per source.
#'
#' @param n_channels Integer. Number of scalp electrodes to simulate (default: 64).
#' @param n_sources Integer. Number of latent brain sources (default: 10).
#' @param sigma_geo Numeric. The width (sigma) of the initial RBF kernel.
#'   Controls the raw electrical field decay before smoothing.
#' @param lambda_smooth Numeric. The regularization parameter for the graph Laplacian.
#'   Higher values create smoother, more biologically plausible topographies
#'   (stronger volume conduction effect).
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{coords_sens}}{Matrix (n_channels x 3). XYZ coordinates of sensors on the unit hemisphere.}
#'   \item{\code{coords_src}}{Matrix (n_sources x 3). XYZ coordinates of sources inside the sphere.}
#'   \item{\code{A_base}}{Matrix (n_channels x n_sources). The final column-normalized, smoothed mixing matrix.}
#'   \item{\code{L_sym}}{Matrix (n_channels x n_channels). The normalized graph Laplacian used for the smoothing operator.}
#' }
#'
#' @importFrom stats dist
#' @export
generate_geometry_mixing <- function(n_channels = 64, n_sources = 10, 
                                     sigma_geo = 0.3, lambda_smooth = 0.5) {
  
  # 1.1 Sensor Coordinates: Fibonacci Grid on Hemisphere
  idx <- 0:(n_channels - 1)
  # Correct mapping: z uniform [0, 1] preserves area
  z <- stats::runif(n_channels, 0, 1)
  theta <- (sqrt(5) * pi * idx) %% (2 * pi)
  r_xy <- sqrt(1 - z^2)
  x <- r_xy * cos(theta)
  y <- r_xy * sin(theta)
  coords_sens <- cbind(x, y, z)
  
  # 1.2 Source Coordinates: Random inside Sphere
  coords_src <- matrix(0, n_sources, 3)
  count <- 0
  while(count < n_sources) {
    pt <- stats::runif(3, -0.8, 0.8)
    if(sum(pt^2) < 0.8^2) {
      count <- count + 1
      coords_src[count, ] <- pt
    }
  }
  
  # 1.3 Raw Weights
  A_raw <- rbf_kernel(x = coords_sens, y = coords_src, 
                      sigma = sigma_geo, standard_scale = TRUE)
  
  # 1.4 Normalized Laplacian Smoothing (Weighted Adjacency)
  dist_SS <- as.matrix(dist(coords_sens))
  k_nn <- 4
  W <- matrix(0, n_channels, n_channels)
  sigma_w <- mean(dist_SS) * 0.5 
  
  for(i in 1:n_channels) {
    nbs <- order(dist_SS[i,])[2:(k_nn+1)]
    W[i, nbs] <- exp(-dist_SS[i, nbs]^2 / (2 * sigma_w^2))
    W[nbs, i] <- W[i, nbs] 
  }
  
  deg <- rowSums(W)
  deg[deg == 0] <- 1 
  D_inv_sqrt <- diag(1 / sqrt(deg))
  L_sym <- diag(n_channels) - D_inv_sqrt %*% W %*% D_inv_sqrt
  
  H <- solve(diag(n_channels) + lambda_smooth * L_sym)
  A_smooth <- H %*% A_raw
  
  col_norms <- sqrt(colSums(A_smooth^2))
  col_norms[col_norms == 0] <- 1
  A_base <- sweep(A_smooth, 2, col_norms, "/")
  
  return(list(coords_sens = coords_sens, coords_src = coords_src, 
              A_base = A_base, L_sym = L_sym))
}
