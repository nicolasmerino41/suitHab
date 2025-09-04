begin
    # assume `pool` exists; else make a quick pool
    if !@isdefined(pool)
        # pool = build_pool(140; basal_frac=0.45, seed=2, Rmin=3.0, Rmax=300.0, b0=0.12, bspread=0.08)
        pool = build_pool(14; basal_frac = 0.5, seed = 10)
    end

    S      = pool.S
    ord    = sortperm(pool.masses)              # increasing mass
    nB     = count(pool.basal)

    # build full adjacency (pred -> prey)
    A = zeros(Int, S, S)
    for s in 1:S
        for q in pool.E[s]
            A[s, q] = 1
        end
    end
    Aord = A[ord, ord]                          # sort by mass

    fig = Figure(; size = (640, 450))
    ax  = Axis(fig[1,1], title = "Metaweb adjacency (sorted by mass)",
               xlabel = "Prey (light → heavy)", ylabel = "Predators (light → heavy)")
    hm  = heatmap!(ax, Aord; colormap = :grayC)
    # mark the basal/consumer boundary (lightest nB are basal)
    vlines!(ax, [nB + 0.5], color=(:red, 0.7), linestyle=:dash)
    hlines!(ax, [nB + 0.5], color=(:red, 0.7), linestyle=:dash)
    display(fig)
end

begin
    # (re)build with the probabilistic kernel
    pool = build_pool(
        140; basal_frac=0.5, seed=2,
        R0_mean=16.0, R0_sd=0.6, sigma=0.55,
        density=0.3, pmax=0.9
    )

    S = pool.S
    B = count(pool.basal)

    # out-degree (prey per species) and consumers-only slice
    prey_counts = [length(pool.E[s]) for s in 1:S]
    prey_cons   = prey_counts[.!pool.basal]

    # in-degree (predators per species)
    pred_counts = zeros(Int, S)
    for s in 1:S, q in pool.E[s]
        pred_counts[q] += 1
    end

    Etotal = sum(prey_counts)
    conn   = Etotal / (S*(S-1))

    println("---- Metaweb summary ----")
    println("Species: ", S, "  (basal=", B, ", consumers=", S-B, ")")
    println("Edges: ", Etotal, "   Connectance: ", round(conn, digits=4))
    println("Prey/consumer mean=", round(mean(prey_cons), digits=2),
            "  q10/50/90=", quantile(prey_cons, [0.1,0.5,0.9]))
    println("Predators/species q10/50/90=",
            quantile(pred_counts, [0.1,0.5,0.9]))

    # proper histograms (binned)
    fig = Figure(; size=(800, 360))
    ax1 = Axis(fig[1,1], title="Prey per consumer", xlabel="# prey", ylabel="count")
    ax2 = Axis(fig[1,2], title="Predators per species", xlabel="# predators", ylabel="count")
    hist!(ax1, prey_cons; bins=15)
    hist!(ax2, pred_counts; bins=15)
    display(fig)
end

begin
    using CairoMakie, Statistics

    # If you don't already have `pool`, make one (using your current build_pool)
    if !@isdefined(pool)
        pool = build_pool(140; basal_frac=0.5, seed=2,
                          R0_mean=16.0, R0_sd=0.6, sigma=0.55,
                          density=0.3, pmax=0.9)
    end

    S = pool.S

    # Out-degree: number of prey per species (then take only consumers)
    prey_counts = [length(pool.E[s]) for s in 1:S]
    prey_cons   = prey_counts[.!pool.basal]

    # In-degree: number of predators per species
    pred_counts = zeros(Int, S)
    for s in 1:S, q in pool.E[s]
        pred_counts[q] += 1
    end

    # Figure: histograms of degrees
    fig = Figure(; size = (900, 380))

    ax1 = Axis(fig[1,1], title = "Out-degree (prey per consumer)",
               xlabel = "# prey", ylabel = "count")
    hist!(ax1, prey_cons; bins = max(10, ceil(Int, sqrt(length(prey_cons)))))

    ax2 = Axis(fig[1,2], title = "In-degree (predators per species)",
               xlabel = "# predators", ylabel = "count")
    hist!(ax2, pred_counts; bins = max(10, ceil(Int, sqrt(length(pred_counts)))))

    display(fig)
end

begin
    using StatsBase, CairoMakie

    S  = pool.S
    B  = count(pool.basal)
    Cn = S - B

    # out-degree (number of prey) per species
    prey_counts = [length(pool.E[s]) for s in 1:S]

    # in-degree (number of predators) per species
    pred_counts = zeros(Int, S)
    for s in 1:S, q in pool.E[s]
        pred_counts[q] += 1
    end

    # total degree = in + out
    degree = prey_counts .+ pred_counts

    # frequency of each degree value
    deg_pairs = collect(pairs(countmap(degree)))
    deg_vals  = first.(deg_pairs)
    deg_freqs = last.(deg_pairs)

    # sort by frequency (descending)
    didx = sortperm(deg_freqs; rev=true)
    xdeg = 1:length(didx)

    fig = Figure(size=(500, 360))
    ax = Axis(fig[1,1], title="Degree distribution", xlabel="degree (#links)", ylabel="count",
              xticks=(xdeg, string.(deg_vals[didx])))

    barplot!(ax, xdeg, deg_freqs[didx]; color=:seagreen)

    display(fig)
end

pool = build_pool(140; basal_frac=0.35, seed=2,
                  R0_mean=16.0, R0_sd=0.6,    # preferred ratio a bit lower
                  sigma=0.35,                  # moderate diet spread
                  density=0.18, pmax=0.7,     # thin links a bit
                  b0=0.12, bspread=0.08)

Z = climate_pass(pool, grid; τ=0.55)          # slightly softer climate filter
P = assemble(Z, pool)
th = trophic_height_per_cell(P, pool)
@show minimum(th), maximum(th), count(>(1), th)/length(th)
