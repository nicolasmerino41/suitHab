#!/usr/bin/env Rscript

library(ggplot2)
library(patchwork)
library(viridis)

set.seed(42)

NX <- 60
NY <- 60

# ============================================================
# Helper: rescale 0–1
# ============================================================

rescale01 <- function(x) (x - min(x)) / (max(x) - min(x))

# ============================================================
# 1) RANDOM FIELD (white noise)
# ============================================================

E_random <- matrix(rnorm(NX * NY), NX, NY)

# ============================================================
# 2) AUTOCORRELATED FIELD
#    Gaussian smoothing via convolution kernel
# ============================================================

# Build Gaussian kernel
make_kernel <- function(size = 15, sigma = 4) {
  ax <- seq(-(size-1)/2, (size-1)/2)
  kernel <- outer(ax, ax, function(x, y) exp(-(x^2 + y^2)/(2*sigma^2)))
  kernel / sum(kernel)
}

kernel <- make_kernel(size = 15, sigma = 4)

# Convolution using fft (fast and clean)
convolve2d <- function(mat, kernel) {
  fft_mat <- fft(mat)
  fft_kernel <- fft(kernel, dim = dim(mat))
  Re(fft(fft_mat * fft_kernel, inverse = TRUE) / length(mat))
}

# Pad kernel to grid size
kernel_pad <- matrix(0, NX, NY)
kdim <- dim(kernel)
kernel_pad[1:kdim[1], 1:kdim[2]] <- kernel

E_auto_raw <- matrix(rnorm(NX * NY), NX, NY)
E_auto <- convolve2d(E_auto_raw, kernel_pad)

# ============================================================
# Rescale both to common range
# ============================================================

E_random <- rescale01(E_random)
E_auto   <- rescale01(E_auto)

common_range <- range(c(E_random, E_auto))

# ============================================================
# Build dataframe
# ============================================================

make_df <- function(mat, label) {
  expand.grid(x = 1:NX, y = 1:NY) |>
    transform(value = as.vector(mat),
              type = label)
}

df <- rbind(
  make_df(E_auto, "Autocorrelated"),
  make_df(E_random, "Random")
)

# ============================================================
# Plot
# ============================================================

p <- ggplot(df, aes(x = x, y = y, fill = value)) +
  geom_raster(interpolate = TRUE) +
  facet_wrap(~type, nrow = 1) +
  scale_fill_viridis(
    option = "viridis",
    limits = common_range,
    name = "Environmental value"
  ) +
  coord_equal(expand = FALSE) +
  theme_minimal(base_size = 15) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 16),
    strip.background = element_rect(fill = "grey92", color = "black"),
    legend.position = "right"
  )

ggsave(
  "SI_environmental_grids_clean.png",
  p,
  width = 11,
  height = 5,
  dpi = 600
)