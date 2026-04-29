# ----------------------------------------------------------------------------
# data-raw/generate_cover.R
#
# Regenerates the README cover figure at man/figures/cover.png.
# The cover is a 2x2 paper-figure-style hero showing the four functions
# the package supports end-to-end:
#
#   (top-left)     Simulate              -- 6-channel SimEEG trace
#   (top-right)    SPD covariance        -- 6 x 6 anchor heatmap
#   (bottom-left)  Augmentation footprint -- E0/G0/G1/G2 vs real anchors
#                                            in tangent-space PC1-PC2
#   (bottom-right) Fidelity audit         -- dimension-normalised LE
#                                            distance bars (E0 > G0 > G1 > G2)
#
# Run from the package root:
#   Rscript data-raw/generate_cover.R
# ----------------------------------------------------------------------------

# Source the package in-place so this script runs without a prior install.
src_files <- list.files("R", pattern = "[.]R$", full.names = TRUE)
for (f in src_files) source(f)

# Load the bundled anchors and labels.
e <- new.env()
load("data/example_anchors.rda", envir = e)
example_anchors <- e$example_anchors
labels <- readRDS("inst/extdata/example_labels.rds")
p <- nrow(example_anchors[[1]])

# ---- Panel data --------------------------------------------------------

# Augmentations (denser n_aug to make the scatter cloud visible).
n_aug <- 8L
aug_e0 <- augment_cov_amplitude_matched_euclidean(
  example_anchors, n_aug = n_aug, g0_sigma = 0.15,
  labels = labels, seed = 1001L
)
aug_g0 <- augment_cov_riemannian(
  example_anchors, n_aug = n_aug, sigma = 0.15,
  labels = labels, seed = 1002L
)
aug_g1 <- augment_cov_empirical_tangent(
  example_anchors, labels, n_aug = n_aug, sigma = 0.15, seed = 1003L
)
aug_g2 <- augment_cov_geodesic_mixup(
  example_anchors, labels, n_aug = n_aug, beta_alpha = 1.0, seed = 1004L
)

# Fidelity audits.
audit_one <- function(aug) {
  audit_covariance_fidelity(example_anchors, aug$cov, aug$anchor)
}
metrics <- list(
  E0 = audit_one(aug_e0),
  G0 = audit_one(aug_g0),
  G1 = audit_one(aug_g1),
  G2 = audit_one(aug_g2)
)

# Tangent-space embedding for the augmentation-footprint scatter:
# vech of the symmetric matrix logarithm, then PCA across all points.
spd_logm_local <- function(C) {
  ev <- eigen(C, symmetric = TRUE)
  ev$vectors %*% diag(log(pmax(ev$values, 1e-12))) %*% t(ev$vectors)
}
triu_mask <- upper.tri(matrix(0, p, p), diag = TRUE)
weight <- ifelse(row(matrix(0, p, p))[triu_mask] ==
                   col(matrix(0, p, p))[triu_mask], 1, sqrt(2))
vech_log <- function(C) {
  Z <- spd_logm_local(C)
  Z[triu_mask] * weight
}

all_cov <- c(example_anchors, aug_e0$cov, aug_g0$cov, aug_g1$cov, aug_g2$cov)
all_kind <- c(rep("real", length(example_anchors)),
              rep("E0",   length(aug_e0$cov)),
              rep("G0",   length(aug_g0$cov)),
              rep("G1",   length(aug_g1$cov)),
              rep("G2",   length(aug_g2$cov)))
M <- t(vapply(all_cov, vech_log, numeric(p * (p + 1L) / 2L)))
pc <- stats::prcomp(M, center = TRUE, scale. = FALSE)
xy <- pc$x[, 1:2]

# Synthetic 6-channel trace for the Simulate panel.
sim <- sim_eeg_master(
  n_trials = 1, n_time = 250, n_channels = 6, n_sources = 4,
  fs = 250, snr_neural_db = 5, snr_artifact_db = 0,
  drift_power_ratio = 0.4,
  target_freqs = c(10, 11, 10, 12),  # one per source
  verbose = FALSE
)
trace <- sim$data[, , 1]

# ---- Render ------------------------------------------------------------

png("man/figures/cover.png",
    width = 1600, height = 1100, res = 180,
    type = "cairo-png", bg = "white")
op <- par(mfrow = c(2, 2), mar = c(3.4, 3.6, 2.4, 1.2),
          mgp = c(2.1, 0.5, 0), tcl = -0.3,
          family = "Helvetica")

# ---- Panel 1: Simulate -------------------------------------------------
fs <- 250
tt <- seq_len(nrow(trace)) / fs
ch_colors <- grDevices::hcl.colors(6, "Plasma", rev = TRUE)
plot(NA, xlim = range(tt), ylim = c(0.3, 6.7),
     xlab = "Time (s)", ylab = "",
     yaxt = "n", main = "Simulate",
     bty = "n", cex.main = 1.3)
for (k in 1:6) {
  s <- scale(trace[, k])[, 1]
  graphics::lines(tt, s / 3 + k, col = ch_colors[k], lwd = 1.4)
  graphics::text(min(tt) - 0.06, k, sprintf("ch%d", k),
                 pos = 2, cex = 0.7, xpd = NA, col = ch_colors[k])
}

# ---- Panel 2: SPD covariance heatmap -----------------------------------
C <- example_anchors[[1]]
# Normalise to [0, 1] for the colour scale so the legend stays comparable.
C_plot <- C[nrow(C):1, ]
graphics::image(seq_len(p), seq_len(p), t(C_plot),
                col = grDevices::hcl.colors(64, "Viridis"),
                xlab = "Channel j", ylab = "Channel i",
                main = "SPD covariance (anchor 1)",
                axes = FALSE, cex.main = 1.3)
graphics::axis(1, at = 1:p, labels = 1:p, cex.axis = 0.7)
graphics::axis(2, at = 1:p, labels = p:1, cex.axis = 0.7, las = 2)
graphics::box()

# ---- Panel 3: Augmentation footprint ------------------------------------
kind_colors <- c(real = "#1b1b1b",
                 E0   = "#d95f02",
                 G0   = "#7570b3",
                 G1   = "#1b9e77",
                 G2   = "#3a9bdc")
plot(xy[, 1], xy[, 2],
     col   = kind_colors[all_kind],
     pch   = ifelse(all_kind == "real", 19, 1),
     cex   = ifelse(all_kind == "real", 1.6, 0.85),
     lwd   = ifelse(all_kind == "real", 1.0, 1.4),
     xlab  = "Tangent PC1",
     ylab  = "Tangent PC2",
     main  = "Augmentation footprint",
     bty   = "n",
     cex.main = 1.3)
graphics::abline(h = 0, v = 0, col = "grey85", lty = 3)
graphics::legend("topright",
       legend = c("real anchor", "E0  off-manifold",
                  "G0  isotropic", "G1  empirical",
                  "G2  geodesic"),
       col = kind_colors[c("real", "E0", "G0", "G1", "G2")],
       pch = c(19, 1, 1, 1, 1),
       cex = 0.65, bty = "n", inset = c(0, 0))

# ---- Panel 4: Fidelity audit -------------------------------------------
methods_lbl <- c("E0", "G0", "G1", "G2")
le_dist <- vapply(metrics, function(m) m$normalized$log_euclidean,
                  numeric(1))
bar_colors <- kind_colors[methods_lbl]
bp <- graphics::barplot(le_dist,
                        names.arg = methods_lbl,
                        col = bar_colors,
                        border = NA,
                        main = "Fidelity audit",
                        ylab = "Normalised LE distance to class mean",
                        ylim = c(0, max(le_dist) * 1.15),
                        cex.main = 1.3, cex.names = 1.0,
                        cex.lab = 0.95)
graphics::text(bp, le_dist + max(le_dist) * 0.04,
               sprintf("%.2f", le_dist), cex = 0.85)

par(op)
dev.off()

cat("Wrote man/figures/cover.png ",
    file.info("man/figures/cover.png")$size, "bytes\n")
cat("Fidelity ranking (normalised LE distance):\n")
print(round(le_dist, 4))
