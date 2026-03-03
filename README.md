# TensorEEG: Physics-Constrained EEG Simulation on Manifolds

## Overview

**TensorEEG** is an R library designed to generate **physically consistent** and **mathematically rigorous** synthetic EEG data in the form of 3rd-order tensors ($\mathcal{X} \in \mathbb{R}^{T \times C \times K}$).

Unlike simple additive noise models, TensorEEG constructs a generative process rooted in:
1.  **Volume Conduction Physics:** Using Normalized Graph Laplacian smoothing on spherical manifolds.
2.  **Manifold Dynamics:** Modeling trial-to-trial variability as geodesic random walks on the Stiefel/Rotation manifold ($SO(n)$).
3.  **Structured Temporal Processes:** Combining VAR(2) background oscillations with time-warped ERP components.

This library is specifically engineered for benchmarking **Tensor Decomposition (CP/PARAFAC2)** algorithms, **Riemannian Geometry** classifiers, and **Source Localization** methods.

---

## Key Features

### 1. Physically Constrained Mixing (Spatial Mode)
We model the volume conduction effect using a **Fibonacci Grid** spherical mapping. The mixing matrix is regularized via **Tikhonov smoothing** over a **Normalized Graph Laplacian**:

$$
\mathbf{A}_{smooth} = (\mathbf{I} + \lambda \mathcal{L}_{sym})^{-1} \mathbf{A}_{raw}
$$

This ensures that the generated topographies exhibit realistic spatial coherence and rank-deficiency common in real EEG recordings.

### 2. Trial-wise Manifold Drift (Trial Mode)
To simulate non-stationarity across trials (e.g., fatigue, electrode shifts), we model the source orientation drift using **Matrix Exponentials** on the tangent space of the mixing matrix:

$$
\mathbf{A}_k = \mathbf{A}_{base} \exp(\theta_k \mathbf{\Omega}_{base})
$$

where $\theta_k$ follows a mean-reverting Ornstein-Uhlenbeck process. This produces data that strictly adheres to **PARAFAC2** structures (shared covariance eigenvalues, evolving eigenvectors).

### 3. Closed-Loop SNR Calibration
TensorEEG implements an **Effective AC Power** ($P_{AC}$) metric using a 4th-order Butterworth high-pass filter ($f_c = 0.1$ Hz). This prevents high-amplitude low-frequency drifts from skewing Signal-to-Noise Ratio (SNR) calculations, ensuring precise control over neural vs. artifact energy.

---

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("Yiming-S/TensorEEG")
```

## Quick Start

```r
library(TensorEEG)

sim <- sim_eeg_master(
  n_trials = 20,
  n_time = 500,
  n_channels = 64,
  n_sources = 10,
  seed = 42,
  verbose = FALSE
)

str(sim$data)
```

## References

1. Harshman, R. A. (1972). *PARAFAC2: Mathematical and technical notes*. UCLA Working Papers in Phonetics.
2. Lütkepohl, H. (2005). *New Introduction to Multiple Time Series Analysis*. Springer.
3. Uhlenbeck, G. E., & Ornstein, L. S. (1930). On the theory of the Brownian motion. *Physical Review*, 36(5), 823-841.
4. Higham, N. J. (2008). *Functions of Matrices: Theory and Computation*. SIAM.
5. Chung, F. R. K. (1997). *Spectral Graph Theory*. American Mathematical Society.
6. Tikhonov, A. N., & Arsenin, V. Y. (1977). *Solutions of Ill-posed Problems*. Winston & Sons.
7. Nunez, P. L., & Srinivasan, R. (2006). *Electric Fields of the Brain: The Neurophysics of EEG* (2nd ed.). Oxford University Press.
8. Barachant, A., Bonnet, S., Congedo, M., & Jutten, C. (2012). Multiclass brain-computer interface classification by Riemannian geometry. *IEEE Transactions on Biomedical Engineering*, 59(4), 920-928.
