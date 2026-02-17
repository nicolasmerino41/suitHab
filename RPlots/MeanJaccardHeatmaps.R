#!/usr/bin/env Rscript

library(tidyverse)
library(viridis)

# ============================================================
# 0) PATHS
# ============================================================

OUTDIR <- "output_biotic_divergence_"

N_CONNECT <- 12
N_CORR    <- 12

Cvals <- seq(0.01, 0.1, length.out = N_CONNECT)
Rvals <- seq(0.0, 0.9, length.out = N_CORR)

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
  mean_jaccard_mismatch = "Mean spatial mismatch\n(1 − Jaccard index)",
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
          all_data[[length(all_data)+1]] <-
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

# ============================================================
# 3) BALANCED THEME
# ============================================================

theme_heat_balanced <- function() {
  theme_classic(base_family = "Arial", base_size = 13) +
    theme(
      plot.margin = margin(12, 18, 12, 18),
      strip.text = element_text(size = 12), #face = "bold"),
      strip.background = element_blank(),
      panel.spacing.x = unit(1.2, "lines"),
      panel.spacing.y = unit(1.6, "lines"),
      panel.border = element_rect(color = "black",
                                  fill = NA,
                                  linewidth = 0.4),
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
# 4) PLOT FUNCTION (FIX: FILTER BY ENVIRONMENT)
# ============================================================
plot_heatmap_metric <- function(metric_name,
                                env_name,
                                fixed_limits = TRUE,
                                horizontal_legend = FALSE) {
  
  df <- heat_df %>%
    filter(Metric == metric_name,
           Environment == env_name)
  
  # Legend configuration
  if (horizontal_legend) {
    
    fill_scale <- scale_fill_viridis_c(
      option = "viridis",
      guide = guide_colorbar(
        direction = "horizontal",   # ← THIS IS THE FIX
        title.position = "top",
        title.hjust = 0.5,
        barwidth  = unit(8, "cm"),
        barheight = unit(0.5, "cm"),
        frame.colour = "black",
        ticks.colour = "black"
      )
    )
    
    legend_pos <- "bottom"
    
  } else {
    
    fill_scale <- scale_fill_viridis_c(
      option = "viridis",
      guide = guide_colorbar(
        direction = "vertical",     # ← explicit
        title.position = "top",
        title.hjust = 0.5,
        label.position = "right",
        barwidth  = unit(0.7, "cm"),
        barheight = unit(6.5, "cm"),
        frame.colour = "black",
        ticks.colour = "black"
      )
    )
    
    legend_pos <- "right"
  }
  
  p <- ggplot(df,
              aes(x = Connectance,
                  y = NicheCorr,
                  fill = value)) +
    geom_tile() +
    facet_grid(Network ~ Regime) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(
      x = "Connectance",
      y = "Niche correlation",
      fill = METRIC_LABELS[[metric_name]]
    ) +
    fill_scale +
    theme_heat_balanced()
  
  # APPLY legend position LAST (critical)
  p <- p + theme(
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
  
  return(p)
}

# ============================================================
# 5) EXPORT (NOW SEPARATED BY ENVIRONMENT)
# ============================================================
metrics_main <- c(
  # "dSrel",
  "mean_jaccard_mismatch"
  # "frac_affected",
  # "realized_overlap"
)

for (env in ENVKINDS) {
  
  for (metric in metrics_main) {
    
    p <- plot_heatmap_metric(metric,
                             env,
                             fixed_limits = TRUE,
                             horizontal_legend = TRUE)
    
    
    ggsave(
      filename = file.path(
        OUTDIR,
        paste0("ggplot_heatmap_", env, "_", metric, ".png")
      ),
      plot = p,
      width = 12,
      height = 10,
      dpi = 600
    )
  }
}

metrics_diag <- c("achieved_r", "Creal")

# for (env in ENVKINDS) {
#   
#   for (metric in metrics_diag) {
#     
#     p <- plot_heatmap_metric(metric, env, fixed_limits = FALSE, horizontal_legend = TRUE)
#     
#     ggsave(
#       filename = file.path(
#         OUTDIR,
#         paste0("ggplot_heatmap_", env, "_", metric, ".png")
#       ),
#       plot = p,
#       width = 14,
#       height = 10,
#       dpi = 600
#     )
#   }
# }

cat("Balanced heatmaps exported separately for random and autocorrelated environments.\n")

