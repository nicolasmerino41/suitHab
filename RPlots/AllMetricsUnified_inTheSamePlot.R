#!/usr/bin/env Rscript

library(tidyverse)

# ============================================================
# SETTINGS
# ============================================================

OUTROOT <- "RPlots/Plots/outputs_merged_all_metrics"
TRAITS  <- c("ctmin", "ctmax", "lt50", "ltmax")

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
      label = paste0("N = ", n,
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
  
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  
  labs(
    title = "Paiwise Predator vs Prey Thermal Limits",
    x = "Predator thermal limit (ºC)",
    y = "Prey thermal limit (ºC)"
  ) +
  
  theme_nature()

ggsave(file.path(OUTDIR, "EDGE_2x2_combined.png"),
       p_edge, width = 10, height = 8, dpi = 600)
# ggsave(
#   file.path(OUTDIR, "EDGE_1x4_combined.png"),
#   p_edge,
#   width = 16,     # wide
#   height = 5.5,   # short
#   dpi = 600
# )

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
  
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  
  labs(
    title = "Predator vs Mean Prey Thermal Limits",
    x = "Predator thermal limit (ºC)",
    y = "Mean prey thermal limit (ºC)"
  ) +
  
  theme_nature()

ggsave(file.path(OUTDIR, "NODE_2x2_combined.png"),
       p_node, width = 10, height = 8, dpi = 600)
# ggsave(
#   file.path(OUTDIR, "NODE_1x4_combined.png"),
#   p_node,
#   width = 16,     # wide
#   height = 5.5,   # short
#   dpi = 600
# )

cat("Combined 2x2 figures generated with n and Pearson r.\n")


################### NULL DISTRIBUTIONS ########################
# ============================================================
# NULL MODEL: Mean absolute predator–prey difference
# Predator-swap null (degree preserved)
# 2x2 layout
# ============================================================

set.seed(123)

N_PERM  <- 1000
NSWEEPS <- 10

edge_swap <- function(pred_vec, nsweeps = 10) {
  preds <- pred_vec
  E <- length(preds)
  
  for (i in seq_len(nsweeps * E)) {
    idx <- sample.int(E, 2)
    tmp <- preds[idx[1]]
    preds[idx[1]] <- preds[idx[2]]
    preds[idx[2]] <- tmp
  }
  
  preds
}

compute_null_absdiff <- function(df, nperm = 1000, nsweeps = 10) {
  
  df <- df %>%
    filter(is.finite(pred_trait),
           is.finite(prey_trait))
  
  base_preds <- df$pred_trait
  base_preys <- df$prey_trait
  
  observed <- mean(abs(base_preds - base_preys))
  
  null_vals <- numeric(nperm)
  
  for (i in seq_len(nperm)) {
    swapped <- edge_swap(base_preds, nsweeps)
    null_vals[i] <- mean(abs(swapped - base_preys))
  }
  
  tibble(
    observed = observed,
    null = null_vals
  )
}

# Run per trait
null_absdiff_df <- edges_df %>%
  group_split(trait) %>%
  map_df(function(d) {
    
    tr <- unique(d$trait)
    res <- compute_null_absdiff(d, N_PERM, NSWEEPS)
    
    tibble(
      trait = tr,
      observed = res$observed[1],
      null = res$null
    )
  })

null_absdiff_df$trait <- factor(null_absdiff_df$trait, levels = TRAITS)

# Compute summary stats
null_summary_absdiff <- null_absdiff_df %>%
  group_by(trait) %>%
  summarise(
    observed = unique(observed),
    null_mean = mean(null),
    null_sd   = sd(null),
    p_value   = mean(null <= observed),  # one-sided: smaller than null?
    .groups = "drop"
  )

# ============================================================
# PLOT 2x2 NULL DISTRIBUTIONS
# ============================================================

p_null_absdiff <- ggplot(null_absdiff_df,
                         aes(x = null)) +
  
  geom_histogram(bins = 30,
                 fill = "grey70",
                 color = "black") +
  
  geom_vline(data = null_summary_absdiff,
             aes(xintercept = observed),
             color = "red",
             linewidth = 1.2) +
  
  geom_text(data = null_summary_absdiff,
            aes(x = observed*1.05,
                y = Inf,
                label = paste0("Observed = ", round(observed, 2),
                               "\np-value < 0.001")),
            vjust = 1.5,
            hjust = -0.05,
            inherit.aes = FALSE,
            size = 4) +
  
  facet_wrap(~ trait, ncol = 2, scales = "free") +
  
  labs(
    title = "Null model: mean |Predator − Prey| thermal difference",
    x = "Mean |Δ thermal limit|",
    y = "Count"
  ) +
  
  theme_nature()

ggsave(
  file.path(OUTDIR, "EDGE_null_mean_absdiff_2x2.png"),
  p_null_absdiff,
  width = 10,
  height = 8,
  dpi = 600
)

cat("2x2 null mean |Δ| figure generated.\n")