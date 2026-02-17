#!/usr/bin/env Rscript

library(tidyverse)

# ============================================================
# PATH
# ============================================================

OUTROOT <- "FebruaryRestart/outputs_merged_all_metrics"

TRAITS <- c("lt50", "ltmin", "ltmax", "ctmin")

# ============================================================
# FUNCTION: scatter + lm + stats + n
# ============================================================

plot_scatter_lm <- function(df, xvar, yvar,
                            xlab, ylab,
                            title,
                            outfile) {
  
  df <- df %>%
    filter(is.finite(.data[[xvar]]),
           is.finite(.data[[yvar]]))
  
  n <- nrow(df)
  if (n < 3) return(NULL)
  
  x <- df[[xvar]]
  y <- df[[yvar]]
  
  # Linear model
  fit <- lm(y ~ x)
  r2  <- summary(fit)$r.squared
  r   <- cor(x, y)
  
  intercept <- coef(fit)[1]
  slope     <- coef(fit)[2]
  
  eq_text <- paste0(
    "y = ", round(intercept, 2),
    ifelse(slope >= 0, " + ", " - "),
    round(abs(slope), 2), "x\n",
    "R² = ", round(r2, 2), "\n",
    "r = ", round(r, 2), "\n",
    "n = ", n
  )
  
  p <- ggplot(df, aes(x = .data[[xvar]],
                      y = .data[[yvar]])) +
    
    geom_point(size = 2) +
    
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 1) +
    
    annotate("text",
             x = min(x),
             y = max(y),
             label = eq_text,
             hjust = 0,
             vjust = 1,
             size = 5) +
    
    labs(x = xlab,
         y = ylab,
         title = title) +
    
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  ggsave(outfile,
         p,
         width = 8,
         height = 6,
         dpi = 600)
  
  return(p)
}

# ============================================================
# LOOP OVER THERMAL METRICS
# ============================================================

for (trait in TRAITS) {
  
  cat("Processing:", trait, "\n")
  
  folder <- file.path(OUTROOT, trait)
  
  edge_file <- file.path(folder, "edge_table.csv")
  node_file <- file.path(folder, "predator_level.csv")
  
  if (!file.exists(edge_file)) next
  
  edges <- read_csv(edge_file, show_col_types = FALSE)
  nodes <- read_csv(node_file, show_col_types = FALSE)
  
  # ----------------------------------------------------------
  # 1) EDGE: predator vs prey
  # ----------------------------------------------------------
  
  plot_scatter_lm(
    edges,
    "pred_trait",
    "prey_trait",
    paste("Predator", trait),
    paste("Prey", trait),
    paste("Edge-level:", trait),
    file.path(folder,
              paste0("R_edge_pred_vs_prey_", trait, ".png"))
  )
  
  # ----------------------------------------------------------
  # 2) EDGE: abs difference vs predator
  # ----------------------------------------------------------
  
  plot_scatter_lm(
    edges,
    "pred_trait",
    "absdiff",
    paste("Predator", trait),
    paste("|Δ|", trait),
    paste("Edge-level |Δ|:", trait),
    file.path(folder,
              paste0("R_edge_pred_vs_absdiff_", trait, ".png"))
  )
  
  # ----------------------------------------------------------
  # 3) NODE: predator vs mean prey
  # ----------------------------------------------------------
  
  if (file.exists(node_file)) {
    plot_scatter_lm(
      nodes,
      "pred_trait",
      "mean_prey_trait",
      paste("Predator", trait),
      paste("Mean prey", trait),
      paste("Node-level:", trait),
      file.path(folder,
                paste0("R_node_pred_vs_meanprey_", trait, ".png"))
    )
  }
}

cat("All R scatter plots generated.\n")

