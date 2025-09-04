"Per-cell climate-suitable share across consumers."
"Share of consumers that have climate suitability at cell c."
function cons_climate_share(Z::BitMatrix, pool::SpeciesPool)
    cons = .!pool.basal
    C = size(Z,2)
    denom = max(1, count(cons))
    share = zeros(Float64, C)
    @inbounds for c in 1:C
        share[c] = count(s -> cons[s] && Z[s,c], 1:size(Z,1)) / denom
    end
    return share
end

"Random mask matched by climate-share quantiles."
function climate_matched_random_mask(grid::Grid, pool::SpeciesPool; τ, keep_frac, nbins::Int=20, seed::Int=11)
    Z = climate_pass(pool, grid; τ=τ)
    share = cons_climate_share(Z, pool)
    C = grid.C
    edges = quantile(share, range(0,1; length=nbins+1))
    bins  = [findall(edges[i] .<= share .< edges[i+1]) for i in 1:nbins]
    rng = MersenneTwister(seed)
    keep = falses(C)
    for idx in bins
        if !isempty(idx)
            k = round(Int, keep_frac * length(idx))
            keep[rand(rng, idx, k)] .= true
        end
    end
    return keep
end

"Within each bin, cluster selections near random seeds using index-distance."
function climate_matched_clustered_mask(grid::Grid, pool::SpeciesPool; τ, keep_frac, nbins::Int=20, nseeds::Int=2, seed::Int=11)
    Z = climate_pass(pool, grid; τ=τ)
    share = cons_climate_share(Z, pool)
    C = grid.C
    edges = quantile(share, range(0,1; length=nbins+1))
    bins  = [findall(edges[i] .<= share .< edges[i+1]) for i in 1:nbins]
    rng = MersenneTwister(seed)
    keep = falses(C)
    for idx in bins
        if isempty(idx); continue; end
        k = round(Int, keep_frac * length(idx))
        # seeds inside the bin
        seeds = rand(rng, idx, min(nseeds, length(idx)))
        # index-distance to nearest seed (cheap 1D surrogate for spatial distance)
        d = [minimum(abs(i - s) for s in seeds) for i in idx]
        order = sortperm(d)
        keep[idx[order[1:k]]] .= true
    end
    return keep
end

# >>> RUN: compare to your earlier masks
Z0 = climate_pass(pool, grid; τ=τ)
kmR = climate_matched_random_mask(grid, pool; τ=τ, keep_frac=keep)
kmC = climate_matched_clustered_mask(grid, pool; τ=τ, keep_frac=keep, nseeds=2)
kmH = climate_matched_hotspot_mask

dR = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmR)
dC = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmC)
@show dR dC  # any dC-dR is now interaction (+ synergy) by design

function control2_curves(pool, grid; τ, losses=0.05:0.05:0.8, nbins=20, nseeds=2)
    ex_clim = Float64[]; ex_inter = Float64[]; ex_syn = Float64[]
    for f in losses
        keep = 1.0 - f
        kmR = climate_matched_random_mask(grid, pool; τ=τ, keep_frac=keep, nbins=nbins, seed=11)
        kmC = climate_matched_clustered_mask(grid, pool; τ=τ, keep_frac=keep, nbins=nbins, nseeds=nseeds, seed=11)
        dR = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmR)
        dC = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmC)
        push!(ex_clim,  dC.clim  - dR.clim)
        push!(ex_inter, dC.inter - dR.inter)
        push!(ex_syn,   dC.syn   - dR.syn)
    end
    return (; losses=collect(losses), ex_clim, ex_inter, ex_syn)
end

# --- plot ---
begin
    grid = make_grid(60,60; seed=11)
    pool = build_pool(200; basal_frac=0.35, seed=1,
        sigma=0.22, density=0.12, pmax=0.70,
        niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.04,
        b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04)

    τ = 0.55
    cc = control2_curves(pool, grid; τ=τ, losses=0.05:0.05:0.8)

    fig = Figure(; size=(640,330))
    ax = Axis(fig[1,1], title="Climate-matched Clustered − Random",
              xlabel="Area lost (fraction)", ylabel="Excess component")
    lines!(ax, cc.losses, cc.ex_clim,  label="Climate")
    lines!(ax, cc.losses, cc.ex_inter, label="Interaction")
    lines!(ax, cc.losses, cc.ex_syn,   label="Synergy")
    hlines!(ax, [0.0], color=(:black,0.3), linestyle=:dot)
    axislegend(ax, position=:lt)
    display(fig)
end

# === climate-matched masks (already have these two) =========================
# cons_climate_share, climate_matched_random_mask, climate_matched_clustered_mask
# ... keep your existing definitions ...

"Climate-matched hotspot mask: same per-bin counts as random, but pick top-share cells in each bin."
function climate_matched_hotspot_mask(grid::Grid, pool::SpeciesPool;
                                      τ, keep_frac, nbins::Int=20, power::Float64=2.0, seed::Int=11)
    Z = climate_pass(pool, grid; τ=τ)
    share = cons_climate_share(Z, pool)        # length = grid.C
    C = grid.C
    edges = quantile(share, range(0,1; length=nbins+1))
    bins  = [findall(edges[i] .<= share .< edges[i+1]) for i in 1:nbins]
    keep  = falses(C)
    for idx in bins
        if isempty(idx); continue; end
        k = round(Int, keep_frac * length(idx))
        if k > 0
            # rank within bin by (share^power), highest first
            sc   = share[idx].^power
            ord  = sortperm(sc; rev=true)
            sel  = idx[ord[1:k]]
            keep[sel] .= true
        end
    end
    return keep
end

"Component-wise excess for climate-matched HOTSPOT vs matched-RANDOM."
function control2_hotspot_curves(pool, grid; τ, losses=0.05:0.05:0.8, nbins=20, power=2.0)
    ex_clim = Float64[]; ex_inter = Float64[]; ex_syn = Float64[]
    for f in losses
        keep = 1.0 - f
        kmR = climate_matched_random_mask(grid, pool; τ=τ, keep_frac=keep, nbins=nbins, seed=11)
        kmH = climate_matched_hotspot_mask(grid, pool; τ=τ, keep_frac=keep, nbins=nbins, power=power, seed=11)
        dR = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmR)
        dH = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmH)
        push!(ex_clim,  dH.clim  - dR.clim)
        push!(ex_inter, dH.inter - dR.inter)
        push!(ex_syn,   dH.syn   - dR.syn)
    end
    return (; losses=collect(losses), ex_clim, ex_inter, ex_syn)
end

# ---- Plot (same style as your clustered control-2) -------------------------
begin
    grid = make_grid(60,60; seed=11)
    pool = build_pool(200; basal_frac=0.35, seed=1,
        sigma=0.22, density=0.12, pmax=0.70,
        niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.04,
        b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04)

    τ = 0.55
    ccH = control2_hotspot_curves(pool, grid; τ=τ, losses=0.05:0.05:0.8, nbins=20, power=2.5)

    fig = Figure(; size=(640,330))
    ax = Axis(fig[1,1], title="Climate-matched Hotspot − Random",
              xlabel="Area lost (fraction)", ylabel="Excess component")
    lines!(ax, ccH.losses, ccH.ex_clim,  label="Climate")
    lines!(ax, ccH.losses, ccH.ex_inter, label="Interaction")
    lines!(ax, ccH.losses, ccH.ex_syn,   label="Synergy")
    hlines!(ax, [0.0], color=(:black,0.3), linestyle=:dot)
    axislegend(ax, position=:lb)
    display(fig)
end
