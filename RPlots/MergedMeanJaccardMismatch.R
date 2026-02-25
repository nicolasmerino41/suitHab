#!/usr/bin/env Rscript

library(tidyverse)
library(viridis)
library(ggh4x)
# ============================================================
# 0) PATHS
# ============================================================

OUTDIR <- "RPlots/Plots/MeanJaccardMismatch_heatmaps/data"

N_CONNECT <- 15
N_CORR    <- 15

Cvals <- seq(0.005, 0.15, length.out = N_CONNECT)
Rvals <- seq(0.0, 1.0, length.out = N_CORR)

ENVKINDS <- c("random", "autocorr")
NETFAMS  <- c("Random", "Modular", "Heavytail", "Cascade")

REGIMES <- c(
  "Narrow niche + low variance",
  "Narrow niche + high variance",
  "Broad niche + low variance",
  "Broad niche + high variance"
)

METRIC_LABELS <- c(
  dSrel = "Relative richness loss\n(1 − S[AB] / S[A])",
  mean_jaccard_mismatch = "Mean Jaccard mismatch",
  frac_affected = "Fraction of consumers affected",
  realized_overlap = "Realized prey-support overlap",
  achieved_r = "Achieved niche correlation",
  Creal = "Realized connectance (L / S²)"
)

# ============================================================
# 1) MATRIX → LONG
# ============================================================

read_matrix_long <- function(file, env, net, reg, metric) {
  
  M <- as.matrix(read.table(file, sep = "\t"))
  
  df <- expand.grid(
    r_index = 1:nrow(M),
    c_index = 1:ncol(M)
  )
  
  df$value <- as.vector(M)
  df$Connectance <- Cvals[df$c_index]
  df$NicheCorr   <- Rvals[df$r_index]
  
  df$Environment <- env
  df$Network     <- net
  df$Regime      <- REGIMES[reg]
  df$Metric      <- metric
  
  df
}

# ============================================================
# 2) LOAD DATA
# ============================================================

all_data <- list()

for (env in ENVKINDS) {
  for (net in NETFAMS) {
    for (ri in 1:4) {
      for (metric in names(METRIC_LABELS)) {
        
        file <- file.path(
          OUTDIR,
          paste0("mat_", env, "_", net, "_reg", ri, "_", metric, ".tsv")
        )
        
        if (file.exists(file)) {
          all_data[[length(all_data) + 1]] <-
            read_matrix_long(file, env, net, ri, metric)
        }
      }
    }
  }
}

heat_df <- bind_rows(all_data)

heat_df$Network <- factor(
  heat_df$Network,
  levels = c("Random", "Modular", "Heavytail", "Cascade")
)

heat_df$Regime <- factor(
  heat_df$Regime,
  levels = REGIMES
)

# ------------------------------------------------------------
# 2b) FILL MISSING CELLS (robust horizontal interpolation)
# ------------------------------------------------------------
heat_df <- heat_df %>%
  arrange(Environment, Network, Regime, Metric,
          NicheCorr, Connectance) %>%
  group_by(Environment, Network, Regime, Metric, NicheCorr) %>%
  mutate(
    value = approx(
      x = Connectance,
      y = value,
      xout = Connectance,
      rule = 2
    )$y
  ) %>%
  ungroup()

# ============================================================
# 3) BALANCED THEME
# ============================================================

theme_heat_balanced <- function() {
  theme_classic(base_family = "Arial", base_size = 13) +
    theme(
      plot.margin = margin(12, 18, 12, 18),
      strip.text = element_text(size = 13),
      strip.background = element_blank(),
      panel.spacing.x = unit(1.2, "lines"),
      panel.spacing.y = unit(1.6, "lines"),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.4
      ),
      axis.title.x = element_text(
        face = "bold",
        size = 15,
        margin = margin(t = 14)
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 15,
        margin = margin(r = 14)
      ),
      axis.text = element_text(size = 12, color = "black"),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      legend.title = element_text(
        face = "bold",
        size = 13,
        margin = margin(b = 10)
      ),
      legend.text = element_text(size = 11),
      legend.position = "right"
    )
}

# ============================================================
# 4) HELPERS
# ============================================================
connectance_axis_labels <- function(x) {
  x_round <- round(x / 0.005) * 0.005
  
  sapply(x_round, function(val) {
    third_decimal <- round(val * 1000) %% 10
    
    if (third_decimal == 0) {
      sprintf("%.2f", val)
    } else {
      sprintf("%.3f", val)
    }
  })
}

make_fill_scale <- function(horizontal_legend = TRUE) {
  
  if (horizontal_legend) {
    list(
      scale = scale_fill_viridis_c(
        option = "viridis",
        oob = scales::squish,
        guide = guide_colorbar(
          direction = "horizontal",
          title.position = "left",
          label.position = "bottom",
          barwidth  = unit(8, "cm"),
          barheight = unit(0.5, "cm"),
          title.theme = element_text(
            face = "bold",
            size = 13,
            hjust = 0.5,
            vjust = 0.9
          ),
          frame.colour = "black",
          ticks.colour = "black"
        )
      ),
      legend_pos = "bottom"
    )
  } else {
    list(
      scale = scale_fill_viridis_c(
        option = "viridis",
        oob = scales::squish,
        guide = guide_colorbar(
          direction = "vertical",
          title.position = "top",
          title.hjust = 0.5,
          label.position = "right",
          barwidth  = unit(0.7, "cm"),
          barheight = unit(6.5, "cm"),
          frame.colour = "black",
          ticks.colour = "black"
        )
      ),
      legend_pos = "right"
    )
  }
}

apply_legend_theme <- function(p, horizontal_legend = TRUE, legend_pos = "bottom") {
  
  if (horizontal_legend) {
    p + theme(
      legend.position = legend_pos,
      legend.title = element_text(
        face = "bold",
        size = 13,
        hjust = 0,
        vjust = 1.2,
        margin = margin(r = 10)
      ),
      legend.text = element_text(size = 11),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
  } else {
    p + theme(
      legend.position = legend_pos,
      legend.title = element_text(
        face = "bold",
        size = 13,
        margin = margin(b = 12)
      ),
      legend.text = element_text(size = 11),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
  }
}

# ============================================================
# 5) MIXED 4-ROW HEATMAP (Random/Heavytail × random/autocorr)
# ============================================================
plot_heatmap_metric_mixed4 <- function(metric_name,
                                       horizontal_legend = TRUE) {
  
  df <- heat_df %>%
    filter(
      Metric == metric_name,
      Environment %in% c("random", "autocorr"),
      Network %in% c("Random", "Heavytail")
    ) %>%
    mutate(
      GridStructure = factor(
        case_when(
          Environment == "random"   ~ "Random environment",
          Environment == "autocorr" ~ "Autocorrelated environment"
        ),
        levels = c("Random environment", "Autocorrelated environment")
      ),
      NetworkStructure = factor(
        case_when(
          Network == "Random"    ~ "Uniform network",
          Network == "Heavytail" ~ "Heavytail network"
        ),
        levels = c("Uniform network", "Heavytail network")
      ),
      Regime = factor(Regime, levels = REGIMES)
    ) %>%
    arrange(GridStructure, NetworkStructure, Regime, NicheCorr, Connectance)
  
  fill_cfg <- make_fill_scale(horizontal_legend = horizontal_legend)
  
  p <- ggplot(df,
              aes(x = Connectance,
                  y = NicheCorr,
                  fill = value)) +
    geom_tile() +
    ggh4x::facet_nested(
      rows = vars(GridStructure, NetworkStructure),
      cols = vars(Regime),
      scales = "fixed"
    ) +
    scale_x_continuous(
      breaks = seq(0.005, 0.15, length.out = 5),
      labels = connectance_axis_labels,
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(
      x = "Connectance",
      y = "Niche correlation",
      fill = METRIC_LABELS[[metric_name]]
    ) +
    fill_cfg$scale +
    theme_heat_balanced() +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y = element_text(size = 12, lineheight = 1.05),
      strip.text.y.left = element_text(angle = 0, hjust = 0),
      ggh4x.facet.nestline = element_line(linewidth = 0.6, colour = "black")
    )
  
  p <- apply_legend_theme(
    p,
    horizontal_legend = horizontal_legend,
    legend_pos = fill_cfg$legend_pos
  )
  
  return(p)
}
# ============================================================
# 6) EXPORT
# ============================================================

# Main target (your requested mixed panel)
p_mix <- plot_heatmap_metric_mixed4(
  metric_name = "mean_jaccard_mismatch",
  horizontal_legend = TRUE
)

ggsave(
  filename = file.path(OUTDIR, "../heatmap_mixed4_mean_jaccard_mismatch.png"),
  plot = p_mix,
  width = 12,
  height = 12,
  dpi = 600
)

# Optional diagnostics (uncomment if you want the same mixed layout for these too)
# for (metric in c("achieved_r", "Creal")) {
#   p_diag <- plot_heatmap_metric_mixed4(metric_name = metric, horizontal_legend = TRUE)
#   ggsave(
#     filename = file.path(OUTDIR, paste0("../heatmap_mixed4_", metric, ".png")),
#     plot = p_diag,
#     width = 12,
#     height = 12,
#     dpi = 600
#   )
# }