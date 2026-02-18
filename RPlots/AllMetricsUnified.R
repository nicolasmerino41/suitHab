#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(ggplot2)

# ============================================================
# PATH
# ============================================================

OUTROOT <- "RPlots/Plots/outputs_merged_all_metrics"

TRAITS <- c("lt50", "ltmin", "ltmax", "ctmin", "ctmax")

# ============================================================
# NATURE-STYLE THEME (MASTER STYLE)
# ============================================================

theme_nature <- function() {
  theme_classic(base_family = "Arial", base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.background = element_blank(),
      panel.border = element_blank()
    )
}

# ============================================================
# SCATTER FUNCTION WITH LM + R2 + r + n
# ============================================================

plot_scatter <- function(df,
                         xvar,
                         yvar,
                         colorvar = NULL,
                         xlab,
                         ylab,
                         title,
                         outfile) {
  
  df <- df %>%
    filter(is.finite(.data[[xvar]]),
           is.finite(.data[[yvar]]))
  
  n <- nrow(df)
  if (n < 3) return(NULL)
  
  x <- df[[xvar]]
  y <- df[[yvar]]
  
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
  
  if (!is.null(colorvar)) {
    p <- ggplot(df,
                aes(x = .data[[xvar]],
                    y = .data[[yvar]],
                    color = .data[[colorvar]])) +
      geom_point(size = 2.5, alpha = 0.9)
  } else {
    p <- ggplot(df,
                aes(x = .data[[xvar]],
                    y = .data[[yvar]])) +
      geom_point(size = 2.5, alpha = 0.9,
                 color = "black")
  }
  
  p <- p +
    geom_smooth(method = "lm",
                se = FALSE,
                color = "black",
                linewidth = 0.8) +
    annotate("text",
             x = min(x),
             y = max(y),
             label = eq_text,
             hjust = 0,
             vjust = 1,
             size = 4.8) +
    labs(x = xlab,
         y = ylab,
         title = title) +
    theme_nature()
  
  ggsave(outfile,
         p,
         width = 9,
         height = 6,
         dpi = 600)
  
  return(p)
}

# ============================================================
# MAIN LOOP
# ============================================================
for (trait in TRAITS) {
  
  cat("Processing:", trait, "\n")
  
  folder <- file.path(OUTROOT, trait)
  
  edge_file    <- file.path(folder, "edge_table.csv")
  pred_file    <- file.path(folder, "predator_level.csv")
  summary_file <- file.path(folder, "summary.csv")
  null_file    <- file.path(folder, "null_vals.csv")
  
  if (!file.exists(edge_file)) next
  
  edges <- read_csv(edge_file, show_col_types = FALSE)
  if ("pred_ctmax" %in% names(edges)) {
    edges <- edges %>%
      rename(
        pred_trait = pred_ctmax,
        prey_trait = prey_ctmax,
        trait_source_pred = ctmax_source_pred
      )
  }
  pred_level <- if (file.exists(pred_file))
    read_csv(pred_file, show_col_types = FALSE) else NULL
  
  summary_df <- if (file.exists(summary_file))
    read_csv(summary_file, show_col_types = FALSE) else NULL
  
  # ------------------------------------------
  # COLUMN NAMES ADAPTATION
  # ------------------------------------------
  
  pred_col  <- "pred_trait"
  prey_col  <- "prey_trait"
  mean_col  <- "mean_prey_trait"
  source_col <- "trait_source_pred"
  
  # ------------------------------------------
  # EDGE: color by interaction source
  # ------------------------------------------
  
  plot_scatter(
    edges,
    pred_col,
    prey_col,
    "interaction_sources",
    paste("Predator", trait),
    paste("Prey", trait),
    paste(toupper(trait)),
    file.path(folder,
              paste0("edge_pred_vs_prey_color_by_interactions.png"))
  )
  
  # ------------------------------------------
  # EDGE: color by trait source
  # ------------------------------------------
  
  if (source_col %in% names(edges)) {
    plot_scatter(
      edges,
      pred_col,
      prey_col,
      source_col,
      paste("Predator", trait),
      paste("Prey", trait),
      paste(toupper(trait)),
      file.path(folder,
                paste0("edge_pred_vs_prey_color_by_trait_source.png"))
    )
  }
  
  # ------------------------------------------
  # EDGE: single colour
  # ------------------------------------------
  
  plot_scatter(
    edges,
    pred_col,
    prey_col,
    NULL,
    paste("Predator", trait),
    paste("Prey", trait),
    paste(toupper(trait)),
    file.path(folder,
              "edge_pred_vs_prey_single_color.png")
  )
  
  # ------------------------------------------
  # NODE-LEVEL
  # ------------------------------------------
  
  pred_level <- if (file.exists(pred_file))
    read_csv(pred_file, show_col_types = FALSE) else NULL
  
  # PATCH FOR OLD CTMAX NODE FILE
  if (!is.null(pred_level) && "pred_ctmax" %in% names(pred_level)) {
    pred_level <- pred_level %>%
      rename(
        pred_trait = pred_ctmax,
        mean_prey_trait = mean_prey_ctmax,
        trait_source_pred = ctmax_source_pred
      )
  }
    
    # single
    plot_scatter(
      pred_level,
      pred_col,
      mean_col,
      NULL,
      paste("Predator", trait),
      paste("Mean prey", trait),
      paste(toupper(trait)),
      file.path(folder,
                "node_pred_vs_meanprey_single_color.png")
    )
}
  
  # ------------------------------------------
  # NULL HISTOGRAM
  # ------------------------------------------

if (!is.null(summary_df) && file.exists(null_file)) {
  
  null_vals <- read_csv(null_file, show_col_types = FALSE)$null_vals
  obs <- summary_df$obs_mean_absdiff[1]
  
  p_null <- ggplot(data.frame(null_vals),
                   aes(x = null_vals)) +
    geom_histogram(bins = 30,
                   fill = "grey70",
                   color = "black") +
    geom_vline(xintercept = obs,
               color = "red",
               linewidth = 1.2) +
    labs(
      title = paste("Null distribution of mean |Δ", trait, "|"),
      x = paste("Mean |Δ", trait, "|"),
      y = "Count"
    ) +
    theme_nature()
  
  ggsave(file.path(folder, "null_mean_absdiff.png"),
         p_null, width = 9, height = 6, dpi = 600)
}

cat("ALL TRAITS COMPLETE.\n")

