library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)

OUTROOT <- "RPlots/Plots"
TRAITS  <- c("ctmin", "lt50", "ctmax", "ltmax")
OUTPUT_FILE <- "RPlots/Plots/MetawebCompilation_PlusCartoon.png"

# ============================================================
# USER CONTROLS
# ============================================================

POS <- list(
  # Reference lines / labels
  x_ref = 10.0,
  average_y_rank_mult = 3.0,
  global_avg_label_x_mult = 1.35,
  global_avg_label_y_mult = 0.85,
  
  # Simulation reference label
  sim_label_x = 10.0,
  sim_label_y = 20.0,
  sim_label_text = "Simulation reference (250 species)\nConnectance = 0.04\nMean Degree = 10.0",
  
  # Cartoon horizontal placement (relative to data xmax)
  cartoon_xmin_mult = 0.62,
  cartoon_xmax_mult = 1.26,
  
  # Cartoon box padding (horizontal only)
  cartoon_box_x_left_div  = 1.03,
  cartoon_box_x_right_mult = 1.03,
  
  # Cartoon vertical frame (absolute in rank coordinates)
  cartoon_box_ymin = 0.5,
  cartoon_box_ymax = 20.0,
  
  # Network layout inside cartoon box
  net_height_frac = 0.30,
  gap_frac = 0.12,
  bottom_base_nudge = 0.20,
  top_base_nudge = 0.20,   # applied before top shift
  bottom_shift = -0.35,    # easy tweak
  top_shift    = -1.2,    # easy tweak
  
  # Cartoon label y offsets (relative to network tops)
  top_cartoon_label_y_nudge    = 0.2,
  bottom_cartoon_label_y_nudge =  0.50,
  
  # Arrow (center x will be computed from cartoon box)
  arrow_y_start = 8.55,
  arrow_y_end   = 12.65,
  
  # Plot limits
  xlim_right_mult = 1.20
)

STY <- list(
  # Main plot
  point_stroke = 0.4,
  lollipop_stem_alpha = 0.25,
  lollipop_stem_lwd = 0.4,
  
  # Cartoon box
  cartoon_box_fill_alpha = 0.65,
  cartoon_box_border_alpha = 0.45,
  cartoon_box_lwd = 0.0,
  
  # Cartoon networks
  cartoon_edge_lwd = 0.58,
  cartoon_edge_alpha = 0.50,
  cartoon_node_size = 4.0,
  cartoon_node_fill = "#f7f7f7",
  cartoon_node_stroke = 0.60,
  
  # Cartoon labels
  cartoon_label_size = 3.3,
  cartoon_label_r = 0.18,
  cartoon_label_pad = 0.5,
  
  # Arrow
  arrow_lwd = 1.5,
  arrow_type = "open",  # "open" or "closed"
  arrow_len_cm = 0.28
)

SHOW <- list(
  ref_line = TRUE,
  sim_reference_label = TRUE
  # If you want markers too, add a helper layer below or toggle them here
)

# ============================================================
# Data
# ============================================================

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

rank_spacing <- 1.0

df <- tibble(
  mean_degree = mean_degree,
  type = trimws(type_of_network),
  n_nodes = n_nodes,
  name = net_names
) |>
  mutate(type = factor(type, levels = legend_order)) |>
  arrange(mean_degree) |>
  mutate(rank = row_number() * rank_spacing)

# ============================================================
# Summary values / axes
# ============================================================

m <- median(df$mean_degree)

xmin <- min(df$mean_degree)
xmax <- max(df$mean_degree)
left_limit <- xmin

df <- df |>
  mutate(
    stem_start = left_limit,
    label_x = mean_degree * 1.08
  )

tickvals <- c(1, 2, 3.3191, 5, 10, 20, 50, 100, 200, 300)
tickvals <- tickvals[tickvals >= xmin & tickvals <= xmax]

average_y <- POS$average_y_rank_mult * rank_spacing

# ============================================================
# Helpers
# ============================================================

make_cartoon_network <- function(connectance = 0.10,
                                 xmin, xmax,
                                 ymin, ymax,
                                 n = 40,
                                 seed = 1) {
  set.seed(seed)
  
  # Clustered layout (looks more cartoon-like than uniform scatter)
  k <- 4
  centers <- tibble(
    cluster = 1:k,
    cx = c(0.20, 0.42, 0.65, 0.82),
    cy = c(0.70, 0.35, 0.65, 0.35)
  )
  
  nodes <- tibble(
    id = 1:n,
    cluster = sample(1:k, size = n, replace = TRUE)
  ) |>
    left_join(centers, by = "cluster") |>
    mutate(
      x = pmin(pmax(rnorm(n(), cx, 0.08), 0.05), 0.95),
      y = pmin(pmax(rnorm(n(), cy, 0.10), 0.05), 0.95),
      x = xmin + x * (xmax - xmin),
      y = ymin + y * (ymax - ymin)
    ) |>
    select(id, cluster, x, y)
  
  pairs <- as.data.frame(t(combn(nodes$id, 2)))
  names(pairs) <- c("i", "j")
  
  edges <- pairs |>
    left_join(nodes |> select(id, cluster, x, y), by = c("i" = "id")) |>
    rename(cluster_i = cluster, x1 = x, y1 = y) |>
    left_join(nodes |> select(id, cluster, x, y), by = c("j" = "id")) |>
    rename(cluster_j = cluster, x2 = x, y2 = y) |>
    mutate(
      p_edge = if_else(cluster_i == cluster_j,
                       pmin(connectance * 2.2, 0.75),
                       pmin(connectance * 0.65, 0.40)),
      keep = runif(n()) < p_edge
    ) |>
    filter(keep) |>
    select(x1, y1, x2, y2)
  
  list(nodes = nodes, edges = edges)
}

label_layers_from_df <- function(lbl_df, size = 3.3, r = 0.18, pad = 0.22) {
  lapply(seq_len(nrow(lbl_df)), function(i) {
    annotate(
      "label",
      x = lbl_df$x[i],
      y = lbl_df$y[i],
      label = lbl_df$label[i],
      hjust = lbl_df$hjust[i],
      fontface = "bold",
      size = size,
      fill = lbl_df$fill[i],
      colour = lbl_df$colour[i],
      label.r = grid::unit(r, "lines"),
      label.padding = grid::unit(pad, "lines")
    )
  })
}

maybe_layer <- function(cond, layer) {
  if (isTRUE(cond)) layer else NULL
}

# ============================================================
# Cartoon placement (easy-to-edit layout)
# ============================================================

# Horizontal extents (in x data units)
cartoon_xmin <- xmax * POS$cartoon_xmin_mult
cartoon_xmax <- xmax * POS$cartoon_xmax_mult

# Vertical frame for the cartoon box (rank coordinates)
cartoon_box_ymin <- POS$cartoon_box_ymin
cartoon_box_ymax <- POS$cartoon_box_ymax

total_height <- cartoon_box_ymax - cartoon_box_ymin
net_height   <- total_height * POS$net_height_frac
gap          <- total_height * POS$gap_frac

# Base network positions inside the cartoon frame
bottom_ymin <- cartoon_box_ymin + POS$bottom_base_nudge
bottom_ymax <- bottom_ymin + net_height

top_ymax <- cartoon_box_ymax - POS$top_base_nudge
top_ymin <- top_ymax - net_height

# User shifts (the main knobs you'll change)
bottom_ymin <- bottom_ymin + POS$bottom_shift
bottom_ymax <- bottom_ymax + POS$bottom_shift

top_ymin <- top_ymin + POS$top_shift
top_ymax <- top_ymax + POS$top_shift

# Cartoon box extents (single box around both networks + arrow)
cartoon_box <- list(
  xmin = cartoon_xmin / POS$cartoon_box_x_left_div,
  xmax = cartoon_xmax * POS$cartoon_box_x_right_mult,
  ymin = cartoon_box_ymin,
  ymax = cartoon_box_ymax
)

# Use geometric center because x-axis is log-scaled (visually centered on log scale)
x_mid <- sqrt(cartoon_box$xmin * cartoon_box$xmax)

# Arrow positions (single center line)
arrow_x <- x_mid
arrow_y_start <- POS$arrow_y_start
arrow_y_end   <- POS$arrow_y_end

# Cartoon label colors (using viridis endpoints)
vir <- viridisLite::viridis(100, option = "viridis")
low_col  <- vir[20]    # dark purple
high_col <- vir[70]  # bright yellow

# higher mean degree = purple (top), lower mean degree = yellow (bottom)
cartoon_labels <- tibble(
  x = c(x_mid, x_mid),
  y = c(top_ymax + POS$top_cartoon_label_y_nudge,
        bottom_ymax + POS$bottom_cartoon_label_y_nudge),
  label = c("HIGHER MEAN DEGREE", "LOWER MEAN DEGREE"),
  fill = c(low_col, high_col),
  colour = c("white", "white"),
  hjust = c(0.5, 0.5)
)

# Build cartoon networks (easy to tweak here)
net_specs <- tibble(
  connectance = c(0.05, 0.15),
  ymin = c(bottom_ymin, top_ymin),
  ymax = c(bottom_ymax, top_ymax),
  n = c(28, 40),
  seed = c(11, 22)
)

cartoon_nets <- lapply(seq_len(nrow(net_specs)), function(i) {
  make_cartoon_network(
    connectance = net_specs$connectance[i],
    xmin = cartoon_xmin,
    xmax = cartoon_xmax,
    ymin = net_specs$ymin[i],
    ymax = net_specs$ymax[i],
    n = net_specs$n[i],
    seed = net_specs$seed[i]
  )
})

cartoon_edges <- bind_rows(lapply(seq_along(cartoon_nets), function(i) {
  cartoon_nets[[i]]$edges |> mutate(group = i)
}))

cartoon_nodes <- bind_rows(lapply(seq_along(cartoon_nets), function(i) {
  cartoon_nets[[i]]$nodes |> mutate(group = i)
}))

# Pre-build repeatable layer lists (no repetition in the main ggplot call)
cartoon_label_layers <- label_layers_from_df(
  cartoon_labels,
  size = STY$cartoon_label_size,
  r = STY$cartoon_label_r,
  pad = STY$cartoon_label_pad
)

optional_ref_line_layer <- maybe_layer(
  SHOW$ref_line,
  geom_vline(
    xintercept = POS$x_ref,
    color = "red",
    linetype = "dashed",
    linewidth = 0.9
  )
)

optional_sim_label_layer <- maybe_layer(
  SHOW$sim_reference_label,
  annotate(
    "label",
    x = POS$sim_label_x,
    y = POS$sim_label_y,
    label = POS$sim_label_text,
    hjust = 0.5,
    size = 3.2,
    fontface = "bold",
    fill = "white",
    colour = "black",
    linewidth = 0.5,
    label.padding = unit(0.45, "lines")
  )
)

# ============================================================
# Plot (single final plotting approach)
# ============================================================
p1 <- ggplot(df, aes(x = mean_degree, y = rank)) +
  
  # ---- Optional reference content ----
optional_ref_line_layer +
  
  # ---- Median line ----
geom_vline(
  xintercept = m,
  linetype = "dashed",
  linewidth = 0.6
) +
  
  # ---- Global average label ----
annotate(
  "label",
  x = m * POS$global_avg_label_x_mult,
  y = average_y * POS$global_avg_label_y_mult,
  label = "GLOBAL AVERAGE",
  fontface = "bold",
  hjust = 0.5,
  size = 4.0,
  fill = "black",
  label.colour = "black",
  colour = "white",
  label.r = grid::unit(0.15, "lines"),
  label.padding = grid::unit(0.25, "lines")
) +
  
  # ---- Lollipop stems ----
geom_segment(
  aes(x = stem_start, xend = mean_degree, y = rank, yend = rank),
  linewidth = STY$lollipop_stem_lwd,
  alpha = STY$lollipop_stem_alpha
) +
  
  # ---- Main points ----
geom_point(
  aes(fill = type, size = n_nodes),
  shape = 21,
  color = "black",
  stroke = STY$point_stroke
) +
  
  # ---- Labels ----
geom_text(
  aes(x = label_x, label = name),
  hjust = 0,
  size = 3.1,
  color = "black",
  show.legend = FALSE
) +
  
  # ---- ONE enclosing cartoon box ----
annotate(
  "rect",
  xmin = cartoon_box$xmin,
  xmax = cartoon_box$xmax,
  ymin = cartoon_box$ymin,
  ymax = cartoon_box$ymax,
  fill = alpha("white", STY$cartoon_box_fill_alpha),
  color = alpha("black", STY$cartoon_box_border_alpha),
  linewidth = STY$cartoon_box_lwd
) +
  
  # ---- Cartoon labels (from data frame -> no repeated annotate blocks) ----
cartoon_label_layers +
  
optional_sim_label_layer +
  
  # ---- Cartoon edges + nodes (combined data -> no repeated geoms) ----
geom_segment(
  data = cartoon_edges,
  aes(x = x1, y = y1, xend = x2, yend = y2),
  inherit.aes = FALSE,
  linewidth = STY$cartoon_edge_lwd,
  alpha = STY$cartoon_edge_alpha,
  lineend = "round",
  colour = "black"
) +
  
  geom_point(
    data = cartoon_nodes,
    aes(x = x, y = y),
    inherit.aes = FALSE,
    shape = 21,
    size = STY$cartoon_node_size,
    fill = STY$cartoon_node_fill,
    colour = "black",
    stroke = STY$cartoon_node_stroke,
    show.legend = FALSE
  ) +
  
  # ---- Double arrow only ----
annotate(
  "segment",
  x = arrow_x,
  xend = arrow_x,
  y = arrow_y_start,
  yend = arrow_y_end,
  linewidth = STY$arrow_lwd,
  lineend = "round",
  colour = "black",
  arrow = arrow(
    type = STY$arrow_type,
    ends = "both",
    length = grid::unit(STY$arrow_len_cm, "cm")
  )
) +
  
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
  
  coord_cartesian(
    xlim = c(left_limit, xmax * POS$xlim_right_mult),
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
    axis.text = element_text(color = "black"),
    plot.margin = margin(10, 20, 10, 10)
  ) +
  
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        size = 4,
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