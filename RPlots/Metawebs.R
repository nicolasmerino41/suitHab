library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)

mean_degree <- c(
  96.1992, 253.1312, 21.3818, 19.0386, 27.4189, 17.5256, 23.7921, 14.7914,
  103.4709, 11.2174, 9.25, 6.4, 16.305, 19.2512, 65.9871, 3.3191, 7.9309, 4.0606,
  4.65909, 14.3603, 22.5812, 34.1791, 8.7313
)

type_of_network <- c(
  "Trophic, parasitic, mutualistic",
  "Trophic", "Trophic", "Trophic", "Parasitic", "Trophic", "Parasitic",
  "Trophic", "Trophic", "Trophic", "Trophic", "Trophic, parasitic",
  "Trophic", "Trophic", "Trophic", "Mutualistic", "Mutualistic", "Mutualistic",
  "Mutualistic", "Trophic, parasitic", "Trophic, parasitic", "Trophic",
  "Mutualistic"
)

n_nodes <- c(
  23020, 1151, 48, 244, 148, 215, 327, 278, 446, 23, 16, 100, 846, 211, 155, 47,
  4838, 66, 88, 1809, 499, 11365, 268
)

net_names <- c(
  "Switzerland", "TETRA-EU", "German Bight", "Barent Sea", "Northern European",
  "European Crop Pests", "Eurasian Rodent Flea", "Brazilian Ant-Tree", "KF metaweb",
  "Cedar Creek Arthropods", "Jena experiment Arthropods", "Forest Plant-Herb-Host",
  "Marine Anctarctic", "Adirondacks lake", "Florida Keys", "Pollen-insect Australia", "Global Plant-frugivore",
  "Southwest China", "Brazilian Forest", "California Rocky Intertidal", "Baja California Estuaries",
  "Marine Global", "Argentinian Pampas"
)

stopifnot(length(mean_degree) == length(type_of_network),
          length(mean_degree) == length(n_nodes),
          length(mean_degree) == length(net_names))

legend_order <- c(
  "Trophic",
  "Parasitic",
  "Mutualistic",
  "Trophic, parasitic",
  "Trophic, parasitic, mutualistic"
)

x_ref <- 6.25

df <- tibble(
  mean_degree = mean_degree,
  type = trimws(type_of_network),
  n_nodes = n_nodes,
  name = net_names
) |>
  mutate(type = factor(type, levels = legend_order)) |>
  arrange(mean_degree) |>
  mutate(
    rank = row_number(),
    xmin = min(mean_degree),
    label_x = mean_degree * 1.04
  )

m   <- median(df$mean_degree)
q25 <- quantile(df$mean_degree, 0.25)
q75 <- quantile(df$mean_degree, 0.75)

xmin <- min(df$mean_degree)
xmax <- max(df$mean_degree)

tickvals <- c(1, 2, 3.3191, 5, 10, 20, 50, 100, 200, 300)
tickvals <- tickvals[tickvals >= xmin & tickvals <= xmax]

pal <- brewer.pal(5, "Set1")
rank_spacing <- 1.4

left_limit <- min(df$mean_degree) * 0.9

df <- df |>
  mutate(
    rank = row_number() * rank_spacing,
    stem_start = left_limit,
    label_x = mean_degree * 1.12
  )

average_y <- 3 * rank_spacing

ggplot(df, aes(x = mean_degree, y = rank)) +
  
  # ---- Reference line FIRST (so it stays behind) ----
geom_vline(xintercept = x_ref,
           color = "red",
           linetype = "dashed",
           linewidth = 0.9) +
  
  # ---- Median line ----
geom_vline(xintercept = m,
           linetype = "dashed",
           linewidth = 0.6) +
  
  # Average label
  annotate("text",
           x = m * 1.05,
           y = average_y*0.85,
           label = "Global average", fontface = "bold",
           hjust = 0,
           size = 4.0) +
  
  # ---- Lollipop stems ----
geom_segment(aes(x = stem_start,
                 xend = mean_degree,
                 y = rank,
                 yend = rank),
             linewidth = 0.4,
             alpha = 0.25) +
  
  # ---- Points (better circles) ----
geom_point(aes(fill = type,
               size = n_nodes),
           shape = 21,
           color = "black",
           stroke = 0.4) +
  
  # ---- Labels ----
geom_text(aes(x = label_x,
              label = name,
              color = "black"),
          hjust = 0,
          size = 3.1,
          show.legend = FALSE) +
  
  # ---- Star + label moved to top ----
annotate("point",
         x = x_ref,
         y = max(df$rank) - rank_spacing*1.1,
         shape = 18,
         size = 6,
         color = "red") +
  
  annotate("text",
           x = x_ref -2.5,
           y = (max(df$rank) - rank_spacing)*0.92,
           label = "Simulation reference (6.25)",
           hjust = 0,
           vjust = -0.5,
           size = 3.2,
           fontface = "bold") +
  
  scale_x_log10(
    breaks = tickvals,
    labels = label_number(accuracy = 1)
  ) +
  
  scale_size_continuous(
    trans = "log10",
    range = c(2, 6),
    guide = "none"
  ) +
  
  scale_fill_viridis_d(
    option = "viridis",
    begin = 0.1,
    end = 0.9
  ) +
  
  scale_color_viridis_d(
    option = "viridis",
    begin = 0.1,
    end = 0.9,
    guide = "none"
  ) +

  coord_cartesian(
  xlim = c(left_limit, xmax * 1.45),

    ylim = c(-0, max(df$rank)),
    clip = "off"
  ) +
  
  labs(
    x = "Mean number of interactions per species",
    y = "Community Rank",
    fill = "Network type"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13, face = "bold"),
    
    axis.line = element_line(linewidth = 0.6),
    
    plot.margin = margin(10, 180, 10, 10)
  ) +
  
  theme(
    axis.text = element_text(color = "black")
  ) +
  
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 4,      # controls legend circle size
        shape = 21,
        stroke = 0.6
      )
    )
  )


