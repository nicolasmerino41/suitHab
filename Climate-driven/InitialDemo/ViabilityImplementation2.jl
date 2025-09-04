# ---- species-wise deltas for a given keep mask ----
function deltas_for_mask(pool::SpeciesPool, grid::Grid; τ::Float64, keepmask)
    Z0 = climate_pass(pool, grid; τ=τ)
    P0 = assemble(Z0, pool)

    A0 = climate_area_per_species(Z0)                       # fraction of kept cells that are clim-suitable
    p0 = bsh1_cond_per_species(P0, Z0, pool)                # prey support within suitable cells

    Z1 = apply_mask(Z0, keepmask)
    P1 = assemble(Z1, pool)

    A1 = climate_area_per_species(Z1)
    p1 = bsh1_cond_per_species(P1, Z1, pool)

    cons = .!pool.basal
    ΔA = A1[cons] .- A0[cons]
    Δp = p1[cons] .- p0[cons]

    keep = .!(isnan.(Δp))                                   # drop consumers with undefined p
    return ΔA[keep], Δp[keep]
end

# utility: choose a keep mask for a given loss fraction and scheme
function keepmask_for(kind::Symbol, grid::Grid, keep::Float64; Z0=nothing, nseeds_cluster::Int=1, seed::Int=123)
    if kind === :random
        return random_mask(grid.C, keep; seed=seed)
    elseif kind === :clustered
        return clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=seed)
    elseif kind === :hotspot
        @assert Z0 !== nothing "pass Z0 (climate_pass) to keepmask_for when kind=:hotspot"
        return hotspot_mask(Z0, current_pool[], keep)        # see note below if you don't keep pool globally
    else
        error("kind must be :random, :clustered, or :hotspot")
    end
end

# simple bin counter for a barplot-based histogram
function hist_counts(x::AbstractVector; nbins::Int=30)
    xmin, xmax = extrema(x)
    if xmin == xmax
        edges = range(xmin - 0.5, xmax + 0.5; length=nbins+1)
    else
        edges = range(xmin, xmax; length=nbins+1)
    end
    counts = zeros(Int, nbins)
    for v in x
        i = clamp(searchsortedlast(edges, v), 1, nbins)
        counts[i] += 1
    end
    mids = (edges[1:end-1] .+ edges[2:end]) ./ 2
    return mids, counts ./ max(1, length(x))   # normalize to frequency
end

# main diagnostic plotter
function diagnostic_plots(kind::Symbol; grid::Grid, pool::SpeciesPool, τ::Float64=0.55,
                          loss_frac::Float64=0.5, nseeds_cluster::Int=1, seed::Int=123)

    # baseline climate pass (for hotspot scoring)
    Z0 = climate_pass(pool, grid; τ=τ)

    keep = 1.0 - loss_frac
    keepmask = if kind === :random
        random_mask(grid.C, keep; seed=seed)
    elseif kind === :clustered
        clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=seed)
    elseif kind === :hotspot
        hotspot_mask(Z0, pool, keep)
    else
        error("kind must be :random, :clustered, or :hotspot")
    end

    ΔA, Δp = deltas_for_mask(pool, grid; τ=τ, keepmask=keepmask)
    K = length(ΔA)
    β = var(ΔA) > 0 ? cov(ΔA, Δp) / var(ΔA) : 0.0
    α = mean(Δp) - β * mean(ΔA)
    ρ = (std(ΔA) > 0 && std(Δp) > 0) ? cor(ΔA, Δp) : 0.0

    # histogram of ΔA
    mids, freq = hist_counts(ΔA; nbins=30)

    begin
        fig = Figure(; size=(920,380))

        # left: ΔA histogram
        ax1 = Axis(fig[1,1], title="ΔA across consumers — $(String(kind)) (loss=$(round(loss_frac,digits=2)))",
                   xlabel="ΔA (change in climate-suitable fraction)", ylabel="frequency")
        barplot!(ax1, mids, freq)
        vlines!(ax1, [0.0], color=:gray, linestyle=:dash)

        # right: ΔA vs Δp with OLS line
        ax2 = Axis(fig[1,2], title="Δp vs ΔA — slope=$(round(β,digits=3)), ρ=$(round(ρ,digits=3))",
                   xlabel="ΔA", ylabel="Δp")
        scatter!(ax2, ΔA, Δp, markersize=5)
        if var(ΔA) > 0
            xs = range(minimum(ΔA), maximum(ΔA); length=100)
            ys = α .+ β .* xs
            lines!(ax2, xs, ys)
        end
        hlines!(ax2, [0.0], color=:gray, linestyle=:dash)
        vlines!(ax2, [0.0], color=:gray, linestyle=:dash)

        display(fig)
    end

    return (; ΔA, Δp, slope=β, intercept=α, corr=ρ, K)
end

# assume you already have grid, pool
diag_rand    = diagnostic_plots(:random;    grid=grid, pool=pool, τ=0.55, loss_frac=0.5, seed=101)
diag_cluster = diagnostic_plots(:clustered; grid=grid, pool=pool, τ=0.55, loss_frac=0.5, nseeds_cluster=1, seed=202)
diag_hotspot = diagnostic_plots(:hotspot;   grid=grid, pool=pool, τ=0.55, loss_frac=0.5)
