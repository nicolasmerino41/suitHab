#!/usr/bin/env Rscript

library(tidyverse)
library(patchwork)
library(viridis)

# ============================================================
# SETTINGS
# ============================================================

OUTROOT <- "RPlots/Plots/outputs_merged_all_metrics"
TRAITS  <- c("ctmin", "ctmax", "lt50", "ltmax")
OUTPUT_FILE <- "RPlots/Plots/NODE_combined_with_niche_colored_background.png"

# ============================================================
# PUBLICATION THEME (background injected later)
# ============================================================
theme_nature <- function(bg_color = "white") {
  theme_classic(base_family = "Arial", base_size = 15) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.title.x = element_text(
        face = "bold",
        margin = margin(t = 10)
      ),
      
      axis.title.y = element_text(
        face = "bold",
        margin = margin(r = 10)
      ),
      axis.text  = element_text(color = "black"),
      legend.position = "none",
      # panel.background = element_rect(
      #   fill = scales::alpha(bg_color, 1.0),
      #   colour = NA
      # )
    )
}

# ============================================================
# NICHE CARTOON FIRST (compute colours here)
# ============================================================
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

# Function to compute viridis colour from x-position
get_label_colour <- function(x_pos) {
  y_curve <- dnorm(x_pos)
  scaled  <- y_curve / ymax
  viridis(100)[round(scaled * 99) + 1]
}

# Manually chosen x-positions
trait_positions <- c(
  ctmin = -2.2,
  lt50  =  2.4,
  ctmax =  1.9,
  ltmax =  2.9
)

# Compute exact label colours
trait_colors <- sapply(trait_positions, get_label_colour)

# ------------------------------------------------------------
# Build niche plot
# ------------------------------------------------------------

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
    limits = c(0, ymax * 1.08),
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
      margin = margin(t = 15), face = "bold"
    ),
    
    axis.title.y = element_text(
      margin = margin(r = 4),
      vjust = 0.5, face = "bold"
    ),
    
    plot.margin = margin(t = 30, r = 10, b = 10, l = 10)
  )

# Add boxed arrows using computed colours

add_metric_label <- function(p, x_pos, y_text, label, fill_col) {
  
  y_curve <- dnorm(x_pos)
  
  p +
    annotate("segment",
             x = x_pos,
             xend = x_pos,
             y = y_text,
             yend = y_curve,
             linewidth = 1.6,
             arrow = arrow(length = unit(0.35, "cm"))) +
    annotate("label",
             x = x_pos,
             y = y_text,
             label = label,
             size = 7,
             fontface = "bold",
             fill = fill_col,
             color = "white",
             label.size = 1.1,
             label.padding = unit(0.35, "lines"))
}

p_niche <- add_metric_label(p_niche, -2.2, 0.31, "CTmin", trait_colors["ctmin"])
p_niche <- add_metric_label(p_niche,  1.9,  0.35, "CTmax",  trait_colors["ctmax"])
p_niche <- add_metric_label(p_niche,  2.4,  0.28, "LT50", trait_colors["lt50"])
p_niche <- add_metric_label(p_niche,  2.9,  0.21, "LTmax", trait_colors["ltmax"])

# ============================================================
# NODE-LEVEL PANEL FUNCTION (uses computed background)
# ============================================================
title_map <- c(
  ctmin = "CTmin",
  ctmax = "CTmax",
  lt50  = "LT50",
  ltmax = "LTmax"
)

make_node_plot <- function(trait, bg_color) {
  
  folder <- file.path(OUTROOT, trait)
  node_file <- file.path(folder, "predator_level.csv")
  if (!file.exists(node_file)) return(NULL)
  
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
  
  if (nrow(df) < 3) return(NULL)
  
  n  <- nrow(df)
  r  <- cor(df$pred_trait, df$mean_prey_trait)
  
  label_text <- paste0("N = ", n,
                       "\nr = ", round(r, 2))
  
  ggplot(df,
         aes(x = pred_trait,
             y = mean_prey_trait)) +
    
    geom_point(size = 2.5, alpha = 0.9, color = bg_color) +
    
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 1,
                color = bg_color) +
    
    annotate("text",
             x = -Inf,
             y = Inf,
             label = label_text,
             hjust = -0.1,
             vjust = 1.2,
             size = 6,
             color = "black") +
    
    labs(
      title = title_map[[trait]],
      x = "Predator thermal limit (°C)",
      y = "Mean prey thermal limit (°C)"
    ) +
    
    theme_nature(bg_color)
}

# ============================================================
# BUILD TOP ROW USING MATCHED COLOURS
# ============================================================

node_plots <- lapply(TRAITS, function(trait) {
  make_node_plot(trait, trait_colors[trait])
})

node_plots <- node_plots[!sapply(node_plots, is.null)]
top_row <- wrap_plots(node_plots, ncol = 4)

# ============================================================
# FINAL LAYOUT
# ============================================================

final_plot <-
  top_row /
  plot_spacer() /
  p_niche +
  plot_layout(heights = c(1.5, 0.0, 0.5))

ggsave(
  OUTPUT_FILE,
  final_plot,
  width = 20,
  height = 8.5,
  dpi = 600
)

cat("Final figure generated with perfectly matched viridis backgrounds.\n")

