using CSV, DataFrames

path = "Metawebs/Compilation.csv"  # <- make sure you're in the right folder

using CairoMakie
using KernelDensity
using Statistics

conns = [0.002, 0.11, 0.255, 0.027, 0.054, 0.084, 0.078, 0.116, 0.116,
         0.244, 0.289, 0.032, 0.01, 0.189, 0.213, 0.163, 0.005, 0.15]

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
              0.244, 0.289, 0.032, 0.01, 0.189, 0.213, 0.163, 0.005, 0.15])

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
    103.4709, 11.2174, 9.25, 6.4, 16.305, 19.2512, 65.9871, 3.3191, 7.9309, 4.0606
])

xs = sort(mean_degree)
n = length(xs)
y = 1:n

m = median(xs)
q25, q75 = quantile(xs, (0.25, 0.75))

# Choose nice log ticks in *data units*
tickvals = [1, 2, 3.3191, 5, 10, 20, 50, 100, 200, 300]
tickvals = filter(t -> minimum(xs) <= t <= maximum(xs), tickvals)
ticklabs = [@sprintf("%.0f", t) for t in tickvals]  # plain numbers, no 10^k

begin
    fig = Figure(size = (850, 480))
    ax = Axis(fig[1, 1];
        xlabel = "Mean number of interactions per species",
        ylabel = "Sorted Rank",
        title  = "Mean prey degree across systems",
        xscale = log10,
        xticks = (tickvals, ticklabs)   # <-- key line
    )

    # Lollipop stems + points
    for (xi, yi) in zip(xs, y)
        lines!(ax, [minimum(xs), xi], [yi, yi]; color = (:black, 0.10))
    end
    scatter!(ax, xs, y; markersize = 10)

    # Median + IQR markers
    vlines!(ax, m; linestyle = :dash, linewidth = 2)
    vspan!(ax, q25, q75; color = (:black, 0.06))

    # --- Simulation reference marker and barrier line ---
    # Barrier line at x_ref
    vlines!(ax, x_ref; color = :red, linewidth = 3, linestyle = :dash)

    # Put the reference point slightly below the data ranks so it stands out
    y_ref = 0.0
    scatter!(ax, [x_ref], [y_ref];
        marker = :star5, markersize = 18, color = :red
    )

    text!(ax, x_ref, y_ref;
        text = " simulation reference (6.25)",
        align = (:left, :center),
        offset = (12, 0)
    )

    xlims!(ax, minimum(xs) * 0.9, maximum(xs) * 1.05)

    display(fig)
end
