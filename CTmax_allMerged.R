#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(ggplot2)

# -----------------------------
# 0) PATHS
# -----------------------------
project_root <- "."
OUTROOT <- file.path(project_root, "FebruaryRestart/outputs_merged_all_ctmax")

edge_file    <- file.path(OUTROOT, "edge_table.csv")
pred_file    <- file.path(OUTROOT, "predator_level.csv")
summary_file <- file.path(OUTROOT, "summary.csv")
null_file    <- file.path(OUTROOT, "null_vals.csv")

# -----------------------------
# 1) LOAD DATA
# -----------------------------
edges <- read_csv(edge_file, show_col_types = FALSE)
pred_level <- read_csv(pred_file, show_col_types = FALSE)
summary_df <- read_csv(summary_file, show_col_types = FALSE)

edges$interaction_sources <- as.factor(edges$interaction_sources)
edges$ctmax_source_pred   <- as.factor(edges$ctmax_source_pred)
pred_level$ctmax_source_pred <- as.factor(pred_level$ctmax_source_pred)

# -----------------------------
# 2) PUBLICATION THEME
# -----------------------------
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

# -----------------------------
# 3) EDGE-LEVEL: color by interaction source
# -----------------------------
p1 <- ggplot(edges,
             aes(x = pred_ctmax,
                 y = prey_ctmax,
                 color = interaction_sources)) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE,
              color = "black", linewidth = 0.8) +
  labs(
    title = "Predator vs Prey CTmax (Edge-level)",
    x = "Predator CTmax",
    y = "Prey CTmax",
    color = "Interaction source"
  ) +
  theme_nature()

ggsave(file.path(OUTROOT,
                 "edge_pred_vs_prey_color_by_interactions.png"),
       p1, width = 9, height = 6, dpi = 600)

# -----------------------------
# 4) EDGE-LEVEL: color by CTmax source
# -----------------------------
p2 <- ggplot(edges,
             aes(x = pred_ctmax,
                 y = prey_ctmax,
                 color = ctmax_source_pred)) +
  geom_point(size = 2.5, alpha = 0.9) +
  geom_smooth(method = "lm", se = FALSE,
              color = "black", linewidth = 0.8) +
  labs(
    title = "Predator vs Prey CTmax (Edge-level)",
    x = "Predator CTmax",
    y = "Prey CTmax",
    color = "CTmax source"
  ) +
  theme_nature()

ggsave(file.path(OUTROOT,
                 "edge_pred_vs_prey_color_by_ctmax_source.png"),
       p2, width = 9, height = 6, dpi = 600)

# -----------------------------
# 5) EDGE-LEVEL: single colour version
# -----------------------------
p2_single <- ggplot(edges,
                    aes(x = pred_ctmax,
                        y = prey_ctmax)) +
  geom_point(size = 2.5, alpha = 0.9,
             color = "black") +
  geom_smooth(method = "lm", se = FALSE,
              color = "black", linewidth = 0.8) +
  labs(
    title = "Predator vs Prey CTmax (Edge-level)",
    x = "Predator CTmax",
    y = "Prey CTmax"
  ) +
  theme_nature()

ggsave(file.path(OUTROOT,
                 "edge_pred_vs_prey_single_color.png"),
       p2_single, width = 9, height = 6, dpi = 600)

# -----------------------------
# 6) NODE-LEVEL
# -----------------------------
p3 <- ggplot(pred_level,
             aes(x = pred_ctmax,
                 y = mean_prey_ctmax,
                 color = ctmax_source_pred)) +
  geom_point(size = 3, alpha = 0.95) +
  geom_smooth(method = "lm", se = FALSE,
              color = "black", linewidth = 0.8) +
  labs(
    title = "Predator vs Mean Prey CTmax (Node-level)",
    x = "Predator CTmax",
    y = "Mean prey CTmax",
    color = "CTmax source"
  ) +
  theme_nature()

ggsave(file.path(OUTROOT,
                 "node_pred_vs_meanprey_color_by_ctmax_source.png"),
       p3, width = 9, height = 6, dpi = 600)

# -----------------------------
# 7) NULL HISTOGRAM
# -----------------------------
obs <- summary_df$obs_mean_absdiff[1]

if (file.exists(null_file)) {
  
  null_vals <- read_csv(null_file, show_col_types = FALSE)$null_vals
  
  p4 <- ggplot(data.frame(null_vals),
               aes(x = null_vals)) +
    geom_histogram(bins = 30,
                   fill = "grey70",
                   color = "black") +
    geom_vline(xintercept = obs,
               color = "red",
               linewidth = 1.2) +
    labs(
      title = "Null distribution of mean |ΔCTmax|",
      x = "Mean |ΔCTmax|",
      y = "Count"
    ) +
    theme_nature()
  
} else {
  
  p4 <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Null values not exported.",
             size = 6) +
    theme_void()
}

ggsave(file.path(OUTROOT,
                 "null_mean_absdiff.png"),
       p4, width = 9, height = 6, dpi = 600)

cat("DONE\n")
