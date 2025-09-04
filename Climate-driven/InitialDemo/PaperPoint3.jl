# --- helper: mean ΔBSH over consumers for a given keep mask -------------
function mean_dBSH_for_mask(pool::SpeciesPool, grid::Grid; τ::Float64, keepmask)
    Z0 = climate_pass(pool, grid; τ=τ)
    P0 = assemble(Z0, pool)
    cons = .!pool.basal
    B0 = mean(bsh1_per_species(P0, Z0, pool)[cons])

    Z1 = apply_mask(Z0, keepmask)
    P1 = assemble(Z1, pool)
    B1 = mean(bsh1_per_species(P1, Z1, pool)[cons])

    return B1 - B0
end

# --- Point 3: excess-damage curves (paired differences) -------------------
"""
excess_damage_curves(pool_seeds, mask_seeds;
    grid, τ, S, basal_frac,
    losses=0.0:0.05:0.8,
    nseeds_cluster=1,
    score_type::Symbol=:vulnerable,
    build_kwargs...,
)

Returns NamedTuple with mean, 10% and 90% quantiles for:
 - ΔB_random, ΔB_clustered, ΔB_hotspot
 - Excess_clustered = clustered - random
 - Excess_hotspot   = hotspot  - random
"""
function excess_damage_curves(pool_seeds, mask_seeds;
    grid::Grid, τ::Float64, S::Int, basal_frac::Float64,
    losses = collect(0.0:0.05:0.8),
    nseeds_cluster::Int = 1,
    score_type::Symbol = :vulnerable,
    build_kwargs...
)
    C = grid.C

    # storage
    μR   = Float64[];  loR   = Float64[];  hiR   = Float64[]
    μC   = Float64[];  loC   = Float64[];  hiC   = Float64[]
    μH   = Float64[];  loH   = Float64[];  hiH   = Float64[]
    μExcC= Float64[];  loExcC= Float64[];  hiExcC= Float64[]
    μExcH= Float64[];  loExcH= Float64[];  hiExcH= Float64[]

    for f in losses
        keep = 1.0 - f
        valsR = Float64[]; valsC = Float64[]; valsH = Float64[]
        exC   = Float64[]; exH   = Float64[]

        for ps in pool_seeds, ms in mask_seeds
            # pool
            pool = build_pool(S; basal_frac=basal_frac, seed=ps, build_kwargs...)

            # hotspot score (depends on the pool)
            score = score_type === :support ?
                prey_support_score(pool, grid; τ=τ) :
                prey_vulnerability_score(pool, grid; τ=τ)

            # paired masks with the same seed
            kmR = random_mask(C, keep; seed=ms)
            kmC = clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=ms)
            kmH = hotspot_clustered_mask_bestfirst(grid, keep; score=score,
                                                   nseeds=nseeds_cluster, seed=ms)

            dR = mean_dBSH_for_mask(pool, grid; τ=τ, keepmask=kmR)
            dC = mean_dBSH_for_mask(pool, grid; τ=τ, keepmask=kmC)
            dH = mean_dBSH_for_mask(pool, grid; τ=τ, keepmask=kmH)

            push!(valsR, dR); push!(valsC, dC); push!(valsH, dH)
            push!(exC, dC - dR); push!(exH, dH - dR)
        end

        append!(μR,   mean(valsR));  append!(loR,   quantile(valsR, 0.10)); append!(hiR,   quantile(valsR, 0.90))
        append!(μC,   mean(valsC));  append!(loC,   quantile(valsC, 0.10)); append!(hiC,   quantile(valsC, 0.90))
        append!(μH,   mean(valsH));  append!(loH,   quantile(valsH, 0.10)); append!(hiH,   quantile(valsH, 0.90))
        append!(μExcC,mean(exC));    append!(loExcC,quantile(exC, 0.10));   append!(hiExcC,quantile(exC, 0.90))
        append!(μExcH,mean(exH));    append!(loExcH,quantile(exH, 0.10));   append!(hiExcH,quantile(exH, 0.90))
    end

    return (
        loss = losses,
        ΔB = (rand=(μ=μR, lo=loR, hi=hiR),
              clust=(μ=μC, lo=loC, hi=hiC),
              hot=(μ=μH, lo=loH, hi=hiH)),
        Excess = (clust=(μ=μExcC, lo=loExcC, hi=hiExcC),
                  hot=(μ=μExcH,  lo=loExcH,  hi=hiExcH))
    )
end

# --- quick plot for Excess curves ------------------------------------------
function plot_excess(ex; title_str="Excess damage (clustered/hotspot − random)")
    begin
        fig = Figure(; size=(760,360))
        ax  = Axis(fig[1,1], title=title_str, xlabel="Area lost (fraction)", ylabel="Excess ΔBSH")

        # ribbons (10–90%)
        band!(ax, ex.loss, ex.Excess.clust.lo, ex.Excess.clust.hi; transparency=true, label="Clustered − Random")
        band!(ax, ex.loss, ex.Excess.hot.lo,   ex.Excess.hot.hi;   transparency=true, label="Hotspot − Random")

        lines!(ax, ex.loss, ex.Excess.clust.μ)
        lines!(ax, ex.loss, ex.Excess.hot.μ)

        hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
        axislegend(ax, position=:lb)
        display(fig)
    end
end

pool_seeds = 1:5
mask_seeds = 1:30
τ = 0.66
losses = collect(0.0:0.05:0.8)

ex = excess_damage_curves(
    pool_seeds, mask_seeds;
    grid=grid, τ=τ, S=220, basal_frac=0.35, losses=losses,
    nseeds_cluster=1, score_type=:vulnerable,
    # low redundancy + high synchrony:
    sigma=0.12, density=0.05, pmax=0.50,
    niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.01,
    b0_basal=0.05, bspread_basal=0.01, b0_cons=0.14, bspread_cons=0.05,
    R0_mean=3.0, R0_sd=0.10
)

plot_excess(ex; title_str="Excess damage curves — LowR/HighS")
