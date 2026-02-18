#!/usr/bin/env Rscript

library(tidyverse)
library(patchwork)
library(viridis)

# ============================================================
# SETTINGS
# ============================================================

OUTROOT <- "RPlots/Plots/outputs_merged_all_metrics"
TRAITS  <- c("ctmin", "lt50", "ctmax", "ltmax")
OUTPUT_FILE <- "NODE_combined_with_niche.png"

# ============================================================
# PUBLICATION THEME
# ============================================================

theme_nature <- function() {
  theme_classic(base_family = "Arial", base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text  = element_text(color = "black"),
      legend.position = "none"
    )
}

# ============================================================
# NODE-LEVEL PANEL FUNCTION
# ============================================================

make_node_plot <- function(trait) {
  
  folder <- file.path(OUTROOT, trait)
  node_file <- file.path(folder, "predator_level.csv")
  if (!file.exists(node_file)) return(NULL)
  
  df <- read_csv(node_file, show_col_types = FALSE)
  
  # Handle ctmax naming
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
  
  if (nrow(df) < 3) return(NULL)
  
  n  <- nrow(df)
  r  <- cor(df$pred_trait, df$mean_prey_trait)
  
  label_text <- paste0("n = ", n,
                       "\nr = ", round(r, 2))
  
  ggplot(df,
         aes(x = pred_trait,
             y = mean_prey_trait)) +
    
    geom_point(size = 2.5, alpha = 0.9, color = "black") +
    
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 1,
                color = "black") +
    
    annotate("text",
             x = -Inf,
             y = Inf,
             label = label_text,
             hjust = -0.1,
             vjust = 1.2,
             size = 6) +
    
    labs(
      title = toupper(trait),
      x = "Predator thermal limit (°C)",
      y = "Mean prey thermal limit (°C)"
    ) +
    
    theme_nature()
}

# ============================================================
# BUILD TOP ROW
# ============================================================

node_plots <- lapply(TRAITS, make_node_plot)
node_plots <- node_plots[!sapply(node_plots, is.null)]
top_row <- wrap_plots(node_plots, ncol = 4)

# ============================================================
# NICHE CARTOON (ULTRA SMOOTH GRADIENT)
# ============================================================

# High resolution for ultra-smooth gradient
x <- seq(-3.5, 3.5, length.out = 5000)
y <- dnorm(x)
ymax <- max(y)

niche_df <- data.frame(x = x, y = y)

dx <- diff(x)[1]

niche_tiles <- niche_df %>%
  mutate(
    height = y,
    width  = dx * 1.02
  )

p_niche <- ggplot(niche_tiles) +
  
  geom_tile(aes(x = x,
                y = height/2,
                height = height,
                width = width,
                fill = height),
            colour = NA) +
  
  geom_line(data = niche_df,
            aes(x = x, y = y),
            linewidth = 1.4,
            color = "black") +
  
  scale_fill_viridis(option = "viridis", guide = "none") +
  
  scale_y_continuous(
    limits = c(0, ymax * 1.08),   # prevent clipping
    expand = c(0, 0)
  ) +
  
  scale_x_continuous(expand = c(0, 0)) +
  
  labs(
    x = "Temperature",
    y = "Performance"
  ) +
  
  theme_classic(base_size = 15) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
  
    axis.title.x = element_text(
      margin = margin(t = 10)   # move x-title further from axis
    ),
  
    axis.title.y = element_text(
      margin = margin(r = -40)   # pull y-title closer to axis
    ),
  
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10)
  )


# ============================================================
# FUNCTION TO ADD BOXED ARROWS
# ============================================================

add_metric_label <- function(p, x_pos, y_text, label) {
  
  y_curve <- dnorm(x_pos)
  scaled  <- y_curve / ymax
  
  box_col <- viridis(100)[round(scaled * 99) + 1]
  
  p +
    annotate("segment",
             x = x_pos,
             xend = x_pos,
             y = y_text,
             yend = y_curve,
             linewidth = 1.5,
             arrow = arrow(length = unit(0.35, "cm"))) +
    
    annotate("label",
             x = x_pos,
             y = y_text,
             label = label,
             size = 7,
             fontface = "bold",
             fill = box_col,
             color = "white",
             label.size = 1.1,
             label.padding = unit(0.3, "lines"))
}

# ============================================================
# ADD METRIC ARROWS (POSITIONED MANUALLY)
# ============================================================

p_niche <- add_metric_label(p_niche, -2.2, 0.26, "CTmin")
p_niche <- add_metric_label(p_niche,  1.7,  0.36, "LT50")
p_niche <- add_metric_label(p_niche,  2.2,  0.26, "CTmax")
p_niche <- add_metric_label(p_niche,  2.7,  0.16, "LTmax")

# ============================================================
# FINAL LAYOUT (WITH SPACER)
# ============================================================

final_plot <-
  top_row /
  plot_spacer() /
  p_niche +
  plot_layout(heights = c(2.2, 0.12, 1.2))

ggsave(
  OUTPUT_FILE,
  final_plot,
  width = 18,
  height = 8.5,
  dpi = 600
)

cat("Final publication-ready figure generated.\n")

