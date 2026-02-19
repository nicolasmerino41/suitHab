library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)

OUTROOT <- "RPlots/Plots"
TRAITS  <- c("ctmin", "lt50", "ctmax", "ltmax")
OUTPUT_FILE <- "RPlots/Plots/MetawebCompilation_PlusCartoon.png"

cartoon_xmin <- xmax * 0.62
cartoon_xmax <- xmax *1.2

cartoon_ymin <- 0.5
cartoon_ymax <- 20

# total height available
total_height <- cartoon_ymax - cartoon_ymin

# give each network equal height
net_height <- total_height * 0.35

# vertical gap for arrow
gap <- total_height * 0.12

# Bottom network box
bottom_ymin <- cartoon_ymin
bottom_ymax <- bottom_ymin + net_height

# Top network box
top_ymax <- cartoon_ymax
top_ymin <- top_ymax - net_height

make_cartoon_network <- function(connectance = 0.1,
                                 xmin, xmax,
                                 ymin, ymax) {
  
  set.seed(1)
  n <- 100
  
  nodes <- tibble(
    id = 1:n,
    x = runif(n),
    y = runif(n)
  )
  
  # scale to rectangle
  nodes <- nodes |>
    mutate(
      x = xmin + x * (xmax - xmin),
      y = ymin + y * (ymax - ymin)
    )
  
  edges <- expand.grid(i = 1:n, j = 1:n) |>
    filter(i < j) |>
    sample_frac(connectance) |>
    left_join(nodes, by = c("i" = "id")) |>
    rename(x1 = x, y1 = y) |>
    left_join(nodes, by = c("j" = "id")) |>
    rename(x2 = x, y2 = y)
  
  list(nodes = nodes, edges = edges)
}

net_sparse <- make_cartoon_network(
  connectance = 0.01,
  xmin = cartoon_xmin,
  xmax = cartoon_xmax,
  ymin = bottom_ymin,
  ymax = bottom_ymax
)

net_dense <- make_cartoon_network(
  connectance = 0.03,
  xmin = cartoon_xmin,
  xmax = cartoon_xmax,
  ymin = top_ymin,
  ymax = top_ymax
)


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
  "Switzerland", "European Tetrapods", "German Bight", "Barent Sea", "Northern European",
  "European Crop Pests", "Eurasian Rodent Flea", "Brazilian Ant-Tree", "South Korea Freshwater",
  "Cedar Creek Arthropods", "Jena experiment Arthropods", "New Zealand Plant-herbivore-host",
  "Marine Antarctic", "Adirondacks Lake", "Florida Keys", "Australia Pollen-insect", "Global Plant-frugivore",
  "Southwest China", "Brazilian Forest", "California Rocky Intertidal", "Baja California Estuaries",
  "Global Marine", "Argentinian Pampas"
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
rank_spacing <- 1.0

left_limit <- min(df$mean_degree)

df <- df |>
  mutate(
    rank = row_number() * rank_spacing,
    stem_start = left_limit,
    label_x = mean_degree * 1.08
  )

average_y <- 3 * rank_spacing

arrow_x <- (cartoon_xmin + cartoon_xmax) / 2

arrow_y_start <- bottom_ymax + gap * 0.3
arrow_y_end   <- top_ymin - gap * 0.1

p1 = ggplot(df, aes(x = mean_degree, y = rank)) +
  
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
  # Average label (black box with white text)
  annotate("label",
           x = m * 1.05,
           y = average_y*0.85,
           label = "Global average",
           fontface = "bold",
           hjust = 0,
           size = 4.0,
           fill = "black",
           label.colour = "black",   # box border
           colour = "white",         # text color
           label.size = 0.25,        # border thickness
           label.r = unit(0.15, "lines"),  # corner roundness
           label.padding = unit(0.25, "lines")
  ) +
  
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
         y = max(df$rank) - rank_spacing*1.7,
         shape = 18,
         size = 6,
         color = "red") +
  
  annotate("text",
           x = x_ref -1.75,
           y =20,
           label = "Simulation reference (6.25)",
           hjust = 0,
           vjust = -0.5,
           size = 3.2,
           fontface = "bold") +
  
# ---- Sparse network ----
geom_segment(data = net_sparse$edges,
           aes(x = x1, y = y1, xend = x2, yend = y2),
           inherit.aes = FALSE,
           linewidth = 0.4) +

geom_point(data = net_sparse$nodes,
           aes(x = x, y = y),
           inherit.aes = FALSE,
           size = 2) +

# ---- Dense network ----
geom_segment(data = net_dense$edges,
           aes(x = x1, y = y1, xend = x2, yend = y2),
           inherit.aes = FALSE,
           linewidth = 0.4) +

geom_point(data = net_dense$nodes,
           aes(x = x, y = y),
           inherit.aes = FALSE,
           size = 2) +

# ---- Double arrow ----
annotate("segment",
         x = arrow_x,
         xend = arrow_x,
         y = arrow_y_start,
         yend = arrow_y_end,
         arrow = arrow(type = "closed",
                       ends = "both",
                       length = unit(0.35, "cm")),
         linewidth = 1.2) +
  
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
      xlim = c(left_limit, xmax * 1.30),
      ylim = c(0, max(df$rank)),
      clip = "off"
    ) +
  
  labs(
    x = "Mean degree (log scale)",
    y = "Rank",
    fill = "Interaction type"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13, face = "bold"),
    
    axis.line = element_line(linewidth = 0.6),
    
    plot.margin = margin(10, 20, 10, 10)
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

ggsave(
  OUTPUT_FILE,
  p1,
  width = 15,
  height = 8.5,
  dpi = 600
)


