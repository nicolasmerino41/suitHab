using CSV, DataFrames

path = "Metawebs/Compilation_extendedManually.csv"  # <- make sure you're in the right folder

using CairoMakie
using KernelDensity
using Statistics

conns = [
    0.002, 0.11, 0.255, 0.027, 0.054, 0.084, 0.078, 0.116, 0.116,
    0.244, 0.289, 0.032, 0.01, 0.189, 0.213, 0.163, 0.005, 0.15, 0.21, 0.0039, 0.0453, 0.055, 0.07085
]

x = Float64.(conns)
xmin, xmax = extrema(x)

# --- choose smaller bins (explicit, predictable) ---
nbins = 15  # try 12–20; 15 is a good start for n=18
edges = range(xmin, xmax; length = nbins + 1)

# --- KDE (evaluate only on [xmin, xmax] so it doesn’t suggest negative connectance) ---
k = kde(x)
grid = range(xmin, xmax; length = 400)
dens = [pdf(k, t) for t in grid]

begin
    fig = Figure(size = (850, 450))
    ax = Axis(fig[1, 1];
        xlabel = "Connectance",
        ylabel = "Density",
        title = "Connectance distribution across systems"
    )

    # Histogram as density
    hist!(ax, x; bins = edges, normalization = :pdf)

    # KDE overlay
    lines!(ax, grid, dens; linewidth = 3)

    # Rug plot (raw data) at y=0 so you see every system
    # vlines!(ax, x; ymin = 0, ymax = 0.03 * maximum(dens), linewidth = 2)

    xlims!(ax, max(-0.01, xmin - 0.02*(xmax-xmin)), xmax + 0.02*(xmax-xmin))

    display(fig)
end
using CairoMakie, Statistics

x = Float64.([0.002, 0.11, 0.255, 0.027, 0.054, 0.084, 0.078, 0.116, 0.116,
              0.244, 0.289, 0.032, 0.01, 0.189, 0.213, 0.163, 0.005, 0.15, 0.21, 0.0039, 0.0453, 0.055, 0.07085])

xs = sort(x)
n = length(xs)

begin
    fig = Figure(size=(800, 450))
    ax = Axis(fig[1,1], xlabel="Connectance", ylabel="System (sorted rank)",
            title="Connectance across systems (sorted)")

    y = 1:n
    hlines!(ax, y; color=(:black, 0.08))             # faint guides
    lines!(ax, xs, y; linewidth=2)                   # optional connecting line
    scatter!(ax, xs, y; markersize=10)

    # Optional: show median as a vertical line
    vlines!(ax, median(xs); linestyle=:dash, linewidth=2)

    display(fig)
end

using CairoMakie
using Statistics
using Printf

mean_degree = Float64.([
    96.1992, 253.1312, 21.3818, 19.0386, 27.4189, 17.5256, 23.7921, 14.7914,
    103.4709, 11.2174, 9.25, 6.4, 16.305, 19.2512, 65.9871, 3.3191, 7.9309, 4.0606,
    4.65909, 14.3603, 22.5812, 34.1791, 8.7313
])

type_of_network = [
    "Trophic, parasitic, mutualistic",
    "Trophic", "Trophic", "Trophic", "Parasitic", "Trophic", "Parasitic",
    "Trophic", "Trophic", "Trophic", "Trophic", "Trophic, parasitic",
    "Trophic", "Trophic", "Trophic", "Mutualistic", "Mutualistic", "Mutualistic",
    "Mutualistic", "Trophic, parasitic", "Trophic, parasitic", "Trophic",
    "Mutualistic"
]

n_nodes = [
    23020, 1151, 48, 244, 148, 215, 327, 278, 446, 23, 16, 100, 846, 211, 155, 47,
    4838, 66, 88, 1809, 499, 11365, 268
]

net_names = [
    "Switzerland", "TETRA-EU", "German Bight", "Barent Sea", "Northern European",
    "European Crop Pests", "Eurasian Rodent Flea", "Brazilian Ant-Tree", "South Korea Freshwater",
    "Cedar Creek Arthropods", "Jena experiment Arthropods", "New Zealand Plant-herbivore-host",
    "Marine Antarctic", "Adirondacks Lake", "Florida Keys", "Pollen-insect Australia", "Global Plant-frugivore",
    "Southwest China", "Brazilian Forest", "California Rocky Intertidal", "Baja California Estuaries",
    "Marine Global", "Argentinian Pampas"
]

@assert length(mean_degree) == length(type_of_network) == length(n_nodes) == length(net_names)

# -----------------------------
# SETTINGS
# -----------------------------
x_ref = 6.25

# Legend order you want (this ALSO defines the color mapping)
legend_order = [
    "Trophic",
    "Parasitic",
    "Mutualistic",
    "Trophic, parasitic",
    "Trophic, parasitic, mutualistic"
]

# Normalize type strings (keeps your current capitalization as-is, just trims whitespace)
normalize_type(s::AbstractString) = strip(s)

# -----------------------------
# SORT, keeping metadata aligned
# -----------------------------
perm = sortperm(mean_degree)
xs = mean_degree[perm]
ts = normalize_type.(type_of_network[perm])
ns = Float64.(n_nodes[perm])
names_sorted = net_names[perm]

n = length(xs)
y = 1:n

m = median(xs)
q25, q75 = quantile(xs, (0.25, 0.75))

tickvals = [1, 2, 3.3191, 5, 10, 20, 50, 100, 200, 300]
tickvals = filter(t -> minimum(xs) <= t <= maximum(xs), tickvals)
ticklabs = [@sprintf("%.0f", t) for t in tickvals]

# -----------------------------
# COLOR MAPPING (fixed to legend_order)
# -----------------------------
missing_types = setdiff(unique(ts), legend_order)
@assert isempty(missing_types) "Types not in legend_order: $(collect(missing_types))"

type_to_idx = Dict(t => i for (i, t) in enumerate(legend_order))
cidx = [type_to_idx[t] for t in ts]  # THIS is used for points and legend

# -----------------------------
# SIZE MAPPING (log-scaled richness)
# -----------------------------
min_px, max_px = 7.0, 22.0
s_raw = log10.(ns)
smin, smax = minimum(s_raw), maximum(s_raw)
markersizes = (s_raw .- smin) ./ (smax - smin + eps()) .* (max_px - min_px) .+ min_px

# -----------------------------
# PLOT
# -----------------------------
begin
    fig = Figure(size = (1200, 700))
    ax = Axis(fig[1, 1];
        xlabel = "Mean number of interactions per species (log-scale)",
        ylabel = "Community Rank",
        # title  = "Mean degree across systems",
        xscale = log10,
        xticks = (tickvals, ticklabs),
        xgridvisible = false,
        ygridvisible = false
    )

    # median + IQR (draw early so it sits behind points)
    # vspan!(ax, q25, q75; color = (:black, 0.06))
    vlines!(ax, m; linestyle = :dash, linewidth = 2)

    # stems
    for (xi, yi) in zip(xs, y)
        lines!(ax, [minimum(xs), xi], [yi, yi]; color = (:black, 0.10))
    end

    # points (colored + sized) — uses the SAME cidx mapping as the legend
    scatter!(ax, xs, y;
        markersize = markersizes,
        color = cidx,
        colormap = :tab10,
        colorrange = (1, length(legend_order))
    )

    # labels next to points
    for i in 1:n
        text!(ax, xs[i], y[i];
            text = names_sorted[i],
            align = (:left, :center),
            offset = (10, 0),
            fontsize = 10
        )
    end

    # simulation reference marker + barrier line
    vlines!(ax, x_ref; color = :red, linewidth = 3, linestyle = :dash)
    scatter!(ax, [x_ref], [0.0]; marker = :star5, markersize = 18, color = :red)
    text!(ax, x_ref, 0.0;
        text = " simulation reference (6.25)",
        align = (:left, :center),
        offset = (12, 0)
    )

    # room for labels and ref point
    xlims!(ax, minimum(xs) * 0.9, maximum(xs) * 1.35)
    ylims!(ax, -1, n + 1)

    # -----------------------------
    # LEGEND UNDERNEATH (colors match because mapping is fixed)
    # We create one dummy handle per legend_order entry.
    # -----------------------------
    legend_handles = Any[]
    for (i, _) in enumerate(legend_order)
        h = scatter!(ax, [NaN], [NaN];
            color = i,
            colormap = :tab10,
            colorrange = (1, length(legend_order)),
            markersize = 10
        )
        push!(legend_handles, h)
    end

    Legend(fig[2, 1], legend_handles, legend_order, "Network type";
        orientation = :horizontal,
        tellwidth = false,
        framevisible = false,
        nbanks = 1
    )

    rowgap!(fig.layout, 10)

    display(fig)

end