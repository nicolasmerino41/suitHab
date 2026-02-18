#!/usr/bin/env Rscript

library(tidyverse)

# ============================================================
# SETTINGS
# ============================================================

OUTROOT <- "RPlots/Plots/outputs_merged_all_metrics"
TRAITS  <- c("ctmin", "lt50", "ctmax", "ltmax")

OUTDIR  <- file.path(OUTROOT, "CombinedFigures")
dir.create(OUTDIR, showWarnings = FALSE)

# ============================================================
# THEME
# ============================================================

theme_nature <- function() {
  theme_classic(base_family = "Arial", base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      legend.position = "none"
    )
}

# ============================================================
# LOAD + STANDARDISE
# ============================================================

load_trait_data <- function(trait) {
  
  folder <- file.path(OUTROOT, trait)
  edge_file <- file.path(folder, "edge_table.csv")
  node_file <- file.path(folder, "predator_level.csv")
  
  if (!file.exists(edge_file)) return(NULL)
  
  edges <- read_csv(edge_file, show_col_types = FALSE)
  
  if ("pred_ctmax" %in% names(edges)) {
    edges <- edges %>%
      rename(
        pred_trait = pred_ctmax,
        prey_trait = prey_ctmax
      )
  }
  
  edges <- edges %>% mutate(trait = trait)
  
  nodes <- NULL
  if (file.exists(node_file)) {
    nodes <- read_csv(node_file, show_col_types = FALSE)
    
    if ("pred_ctmax" %in% names(nodes)) {
      nodes <- nodes %>%
        rename(
          pred_trait = pred_ctmax,
          mean_prey_trait = mean_prey_ctmax
        )
    }
    
    nodes <- nodes %>% mutate(trait = trait)
  }
  
  list(edges = edges, nodes = nodes)
}

all_edges <- list()
all_nodes <- list()

for (tr in TRAITS) {
  dat <- load_trait_data(tr)
  if (!is.null(dat)) {
    all_edges[[tr]] <- dat$edges
    if (!is.null(dat$nodes))
      all_nodes[[tr]] <- dat$nodes
  }
}

edges_df <- bind_rows(all_edges)
nodes_df <- bind_rows(all_nodes)

edges_df$trait <- factor(edges_df$trait, levels = TRAITS)
nodes_df$trait <- factor(nodes_df$trait, levels = TRAITS)

# ============================================================
# FUNCTION: Compute annotation per trait
# ============================================================

compute_stats <- function(df, xvar, yvar) {
  
  df %>%
    filter(is.finite(.data[[xvar]]),
           is.finite(.data[[yvar]])) %>%
    group_by(trait) %>%
    summarise(
      n = n(),
      r = cor(.data[[xvar]], .data[[yvar]]),
      label = paste0("n = ", n,
                     "\nr = ", round(r, 2)),
      .groups = "drop"
    )
}

# ============================================================
# EDGE STATS
# ============================================================

edge_stats <- compute_stats(edges_df, "pred_trait", "prey_trait")

# ============================================================
# NODE STATS
# ============================================================

node_stats <- compute_stats(nodes_df, "pred_trait", "mean_prey_trait")

# ============================================================
# EDGE-LEVEL 2x2 FACET
# ============================================================

p_edge <- ggplot(edges_df,
                 aes(x = pred_trait,
                     y = prey_trait)) +
  
  geom_point(size = 2, alpha = 0.9, color = "black") +
  
  geom_smooth(method = "lm",
              se = FALSE,
              linewidth = 1,
              color = "black") +
  
  geom_text(data = edge_stats,
            aes(label = label),
            x = -Inf,
            y = Inf,
            hjust = -0.1,
            vjust = 1.1,
            inherit.aes = FALSE,
            size = 5) +
  
  facet_wrap(~ trait, ncol = 4, scales = "free") +
  
  labs(
    title = "Paiwise Predator vs Prey Thermal Limits",
    x = "Predator thermal limit (ºC)",
    y = "Prey thermal limit (ºC)"
  ) +
  
  theme_nature()

# ggsave(file.path(OUTDIR, "EDGE_2x2_combined.png"),
#        p_edge, width = 10, height = 8, dpi = 600)
ggsave(
  file.path(OUTDIR, "EDGE_1x4_combined.png"),
  p_edge,
  width = 16,     # wide
  height = 5.5,   # short
  dpi = 600
)

# ============================================================
# NODE-LEVEL 2x2 FACET
# ============================================================
p_node <- ggplot(nodes_df,
                 aes(x = pred_trait,
                     y = mean_prey_trait)) +
  
  geom_point(size = 2.5, alpha = 0.95, color = "black") +
  
  geom_smooth(method = "lm",
              se = FALSE,
              linewidth = 1,
              color = "black") +
  
  geom_text(data = node_stats,
            aes(label = label),
            x = -Inf,
            y = Inf,
            hjust = -0.1,
            vjust = 1.1,
            inherit.aes = FALSE,
            size = 5) +
  
  facet_wrap(~ trait, ncol = 4, scales = "free") +
  
  labs(
    title = "Predator vs Mean Prey Thermal Limits",
    x = "Predator thermal limit (ºC)",
    y = "Mean prey thermal limit (ºC)"
  ) +
  
  theme_nature()

ggsave(file.path(OUTDIR, "NODE_2x2_combined.png"),
       p_node, width = 10, height = 8, dpi = 600)
ggsave(
  file.path(OUTDIR, "NODE_1x4_combined.png"),
  p_node,
  width = 16,     # wide
  height = 5.5,   # short
  dpi = 600
)

cat("Combined 2x2 figures generated with n and Pearson r.\n")

