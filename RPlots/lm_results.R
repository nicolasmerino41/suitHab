#!/usr/bin/env Rscript

library(tidyverse)
library(broom)
library(writexl)

# ============================================================
# SETTINGS
# ============================================================

OUTROOT <- "RPlots/Plots/outputs_merged_all_metrics"
TRAITS  <- c("ctmin", "ctmax", "lt50", "ltmax")

results_list <- list()

# ============================================================
# LOOP OVER METRICS
# ============================================================

for (trait in TRAITS) {
  
  folder <- file.path(OUTROOT, trait)
  node_file <- file.path(folder, "predator_level.csv")
  
  if (!file.exists(node_file)) next
  
  df <- read_csv(node_file, show_col_types = FALSE)
  
  if (trait == "ctmax") {
    df <- df %>%
      rename(
        pred_trait      = pred_ctmax,
        mean_prey_trait = mean_prey_ctmax
      )
  }
  
  df <- df %>%
    filter(is.finite(pred_trait),
           is.finite(mean_prey_trait))
  
  if (nrow(df) < 3) next
  
  # Fit linear model
  model <- lm(mean_prey_trait ~ pred_trait, data = df)
  
  tidy_mod  <- tidy(model)
  glance_mod <- glance(model)
  
  slope_row <- tidy_mod %>% filter(term == "pred_trait")
  int_row   <- tidy_mod %>% filter(term == "(Intercept)")
  
  results_list[[trait]] <- tibble(
    Metric      = toupper(trait),
    N           = nrow(df),
    Formula     = "mean_prey_trait ~ pred_trait",
    Intercept   = int_row$estimate,
    SE_Intercept= int_row$std.error,
    Slope       = slope_row$estimate,
    SE_Slope    = slope_row$std.error,
    t_Slope     = slope_row$statistic,
    p_Slope     = slope_row$p.value,
    R_squared   = glance_mod$r.squared
  )
}

# Combine all results
results_table <- bind_rows(results_list)

# ============================================================
# EXPORT TO EXCEL
# ============================================================

write_xlsx(results_table, "Thermal_alignment_lm_results.xlsx")

cat("Excel table generated: Thermal_alignment_lm_results.xlsx\n")
