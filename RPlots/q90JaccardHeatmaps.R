#!/usr/bin/env Rscript

library(tidyverse)
library(viridis)

# ============================================================
# 0) PATHS
# ============================================================

OUTDIR <- "output_jaccard_tail_60x60_"

N_CONNECT <- 15
N_CORR    <- 15

Cvals <- seq(0.005, 0.1, length.out = N_CONNECT)
Rvals <- seq(0.0, 0.9, length.out = N_CORR)

ENVKINDS <- c("random", "autocorr")
NETFAMS  <- c("random", "modular", "heavytail", "cascade")

REGIMES <- c(
  "Narrow niche + low variance",
  "Narrow niche + high variance",
  "Broad niche + low variance",
  "Broad niche + high variance"
)

METRIC_LABELS <- c(
  mismatch_q90 = "90th percentile of spatial mismatch\n(1 − Jaccard index)",
  mismatch_frac_gt = "Fraction of consumers with strong mismatch\n(1 − J > 0.8)"
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
  levels = c("random", "modular", "heavytail", "cascade")
)

heat_df$Regime <- factor(
  heat_df$Regime,
  levels = REGIMES
)

# ============================================================
# 3) FINAL POLISHED THEME (IDENTICAL STYLE)
# ============================================================

theme_heat_final <- function() {
  theme_classic(base_family = "Arial", base_size = 14) +
    theme(
      plot.margin = margin(20, 30, 20, 30),
      
      strip.text = element_text(size = 15, face = "bold"),
      strip.background = element_blank(),
      
      panel.spacing.x = unit(1.5, "lines"),
      panel.spacing.y = unit(2, "lines"),
      
      panel.border = element_rect(color = "black",
                                  fill = NA,
                                  linewidth = 0.4),
      
      axis.title.x = element_text(
        face = "bold",
        size = 16,
        margin = margin(t = 16)
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 16,
        margin = margin(r = 16)
      ),
      
      axis.text = element_text(size = 13, color = "black"),
      
      legend.position = "right",
      legend.title = element_text(
        face = "bold",
        size = 14,
        margin = margin(b = 14)
      ),
      legend.text = element_text(size = 12),
      legend.key = element_blank()
    )
}

# ============================================================
# 4) HEATMAP FUNCTION (SQUARED PANELS)
# ============================================================

plot_heatmap_metric <- function(metric_name,
                                env_name,
                                fixed_limits = TRUE) {
  
  df <- heat_df %>%
    filter(Metric == metric_name,
           Environment == env_name)
  
  p <- ggplot(df,
              aes(x = Connectance,
                  y = NicheCorr,
                  fill = value)) +
    
    geom_tile() +
    
    facet_grid(Network ~ Regime) +
    
    coord_fixed() +   # square panels
    
    scale_x_continuous(
      breaks = c(0.005, 0.04, 0.07, 0.10),
      labels = sprintf("%.2f", c(0.005, 0.04, 0.07, 0.10)),
      expand = c(0,0)
    ) +
    
    scale_y_continuous(expand = c(0,0)) +
    
    labs(
      x = "Connectance",
      y = "Niche correlation",
      fill = METRIC_LABELS[[metric_name]]
    ) +
    
    scale_fill_viridis_c(
      option = "viridis",
      limits = if (fixed_limits) c(0,1) else NULL,
      guide = guide_colorbar(
        direction = "vertical",
        title.position = "top",
        title.hjust = 0.5,
        barwidth  = unit(0.9, "cm"),
        barheight = unit(7, "cm"),
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +
    
    theme_heat_final()
  
  return(p)
}

# ============================================================
# 5) EXPORT
# ============================================================

for (env in ENVKINDS) {
  
  for (metric in names(METRIC_LABELS)) {
    
    p <- plot_heatmap_metric(metric,
                             env,
                             fixed_limits = TRUE)
    
    ggsave(
      filename = file.path(
        OUTDIR,
        paste0("ggplot_heatmap_", env, "_", metric, ".png")
      ),
      plot = p,
      width = 16,
      height = 11,
      dpi = 600
    )
  }
}

cat("Tail heatmaps (q90 and frac) exported separately for both environments.\n")

