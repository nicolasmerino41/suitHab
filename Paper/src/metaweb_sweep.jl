# src/metaweb_sweep.jl
# Sweep metaweb structure by thinning/capping prey lists and compare AM vs BAM.
# Uses only your package APIs: Grids, Metawebs, HL, BSH.

module MetawebSweep

using Random, Statistics
using CairoMakie
using ..Grids
using ..Metawebs
using ..BSH
using ..HL

export run_metaweb_sweep
export run_realities_Aonly

# ------------------------ utilities on SpeciesPool ----------------------------

"Deep copy a metaweb (fields you use)."
function _copy_pool(pool::Metawebs.SpeciesPool)
    Metawebs.SpeciesPool(
        pool.S,
        copy(pool.masses),
        copy(pool.basal),
        [copy(v) for v in pool.prey],
    )
end

"""
    thin_cap_metaweb!(pool; keep_prob, diet_cap, rng)

Bernoulli-thin each existing link with probability `keep_prob` (0–1),
then cap each consumer diet at `diet_cap`. Ensures ≥1 prey per consumer by
back-filling the nearest lower-mass prey.
"""
function thin_cap_metaweb!(pool::Metawebs.SpeciesPool; keep_prob::Float64, diet_cap::Int, rng::AbstractRNG)
    S = pool.S
    order = sortperm(pool.masses)             # prey→pred order
    rank  = zeros(Int, S); rank[order] = 1:S

    @inbounds for s in 1:S
        pool.basal[s] && continue
        # Bernoulli thinning
        if keep_prob < 0.999
            kept = Int[]
            for q in pool.prey[s]
                (rand(rng) < keep_prob) && push!(kept, q)
            end
            pool.prey[s] = kept
        end
        # hard cap
        if isfinite(diet_cap) && length(pool.prey[s]) > diet_cap
            shuffle!(rng, pool.prey[s])
            pool.prey[s] = pool.prey[s][1:diet_cap]
        end
        # ensure ≥1 prey (nearest lower-mass fallback)
        if isempty(pool.prey[s])
            r = rank[s]
            if r > 1
                push!(pool.prey[s], order[r-1])
            else
                bas = findall(pool.basal)
                isempty(bas) || push!(pool.prey[s], bas[1])
            end
        end
    end
    return pool
end

# ----------------------------- sweep driver -----------------------------------

"""
    run_metaweb_sweep(grids; S, basal_frac, archetype=:mid,
                      loss_fracs=0.2:0.1:0.8, fstar=0.6,
                      keep_probs=[1.0,0.7,0.5,0.3], caps=[9999,6,4,2],
                      pars=BSH.BAMParams(...), outdir="figs/sweep")

- `grids` :: Vector of (name::String, grid::Grids.Grid) you already built in run_all.
- Builds a base pool per grid with `Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype)`.
- For each (keep_prob, cap): copies the pool, applies `thin_cap_metaweb!`,
  computes relative loss curves via `BSH.relative_loss_curves`, and at `f*`
  records DiD = (relBAM - relAM) and whether the worst HL geometry flips.

Produces a heatmap PNG per grid and returns a Dict of results.
"""
function run_metaweb_sweep(; 
    grids::Vector{Tuple{String,Grids.Grid}},
    S::Int, basal_frac::Float64, archetype::Symbol=:mid,
    loss_fracs=0.2:0.1:0.8, fstar::Float64=0.6,
    keep_probs=[1.0,0.7,0.5,0.3],
    caps=[9999,6,4,2],
    pars::BSH.BAMParams=BSH.BAMParams(; τA=0.5, τB=0.35, τocc=0.2, γ=3.0, movement=:off, T=8),
    outdir::String="Paper/figs/sweep"
)

    isdir(outdir) || mkpath(outdir)
    xs = collect(loss_fracs)
    # pick nearest index to f*
    at_index = argmin(abs.(xs .- fstar))

    out = Dict{String,Any}()

    for (gname, grid) in grids
        rng_pool = MersenneTwister(hash((:pool, gname, S, basal_frac, archetype)))
        base_pool = Metawebs.build_metaweb_archetype(rng_pool; S=S, basal_frac=basal_frac, archetype=archetype)

        # store as (cols=caps, rows=keep_probs) to match heatmap(x=caps, y=keep_probs, Z)
        DiD   = zeros(Float64, length(caps), length(keep_probs))
        Rflip = falses(       length(caps), length(keep_probs))

        for (iy, kp) in enumerate(keep_probs), (ix, cap) in enumerate(caps)
            pool = _copy_pool(base_pool)
            thin_cap_metaweb!(
                pool; keep_prob=kp, diet_cap=cap,
                rng=MersenneTwister(hash((:thin, gname, kp, cap)))
            )

            rel = BSH.relative_loss_curves(
                MersenneTwister(hash((:curves, gname, kp, cap))),
                pool, grid, pars;
                loss_fracs=loss_fracs, seed_A=1
            )

            # DiD at f*: average across geometries
            dids = Float64[]
            # track which geometry is worst under AM and under BAM
            worstAM, worstBAM = "", ""
            bestAM, bestBAM   = +Inf, +Inf # more negative is "worse", so compare raw values

            for geom in (:random, :clustered, :front)
                relAM  = rel[geom].relAM[at_index]
                relBAM = rel[geom].relBAM[at_index]
                push!(dids, relBAM - relAM)

                if relAM < bestAM;   bestAM  = relAM;  worstAM  = String(geom); end
                if relBAM < bestBAM; bestBAM = relBAM; worstBAM = String(geom); end
            end

            DiD[ix, iy]   = mean(dids)
            Rflip[ix, iy] = (worstAM != worstBAM)
        end

        out[gname] = (; DiD, Rflip, keep_probs, caps, fstar, loss_fracs)

        # ---- plot for this grid ----
        fig = Figure(; size=(1050,380))

        ax1 = Axis(fig[1,1], title="DiD (BAM−AM) at f*=$(round(fstar, digits=2))",
                   xlabel="diet cap", ylabel="keep prob")
        heatmap!(ax1, 1:length(caps), 1:length(keep_probs), DiD; colormap=:viridis)
        ax1.xticks = (1:length(caps), string.(caps))
        ax1.yticks = (1:length(keep_probs), string.(keep_probs))
        Colorbar(fig[1,2], label="DiD")

        ax2 = Axis(fig[1,3], title="Rank flip (worst HL under BAM ≠ AM)",
                   xlabel="diet cap", ylabel="keep prob")
        heatmap!(ax2, 1:length(caps), 1:length(keep_probs), Float64.(Rflip); colormap=[:white, :orange])
        ax2.xticks = (1:length(caps), string.(caps))
        ax2.yticks = (1:length(keep_probs), string.(keep_probs))

        save(joinpath(outdir, "metaweb_sweep_$(gname)_f$(Int(round(100*fstar))).png"), fig)
    end

    out
end

"""
run_realities_Aonly(; rng, grid, realities=:default, loss_fracs=0.2:0.1:0.8,
                    S=200, basal_frac=0.25, seed_A=1, Npools=12, pool_seed0=1,
                    geoms=(:random,:clustered,:front))

Returns Dict{Symbol,NamedTuple} mapping each geometry to (loss, ABM, MAB, BAM),
each a vector of A-only means (mean over consumers / original area).
"""
function run_realities_Aonly(; rng, grid, realities=:default,
    loss_fracs=0.2:0.1:0.8, S::Int=200, basal_frac::Float64=0.25,
    seed_A::Int=1, Npools::Int=12, pool_seed0::Int=1,
    geoms::Tuple=(:random,:clustered,:front), τA::Float64=0.5)

    # --- define three pool generators using your Metawebs API
    build_ABM = () -> Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:low)   # tight niches / lower R
    build_MAB = () -> Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:mid)   # broader niches
    build_BAM = () -> Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:high)  # high R (biotic structure)

    out = Dict{Symbol,NamedTuple}()

    for g in geoms
        A_ABM = Float64[]; A_MAB = Float64[]; A_BAM = Float64[]
        for f in loss_fracs
            keepfrac = 1 - f
            keep = g === :random    ? HL.random_mask(rng, grid.C, keepfrac) :
                   g === :clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                      HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)

            vals_ABM = Float64[]; vals_MAB = Float64[]; vals_BAM = Float64[]

            for k in 0:Npools-1
                # new pool per replicate
                pool_ABM = build_ABM()
                pool_MAB = build_MAB()
                pool_BAM = build_BAM()

                # same A-building path as BSH uses (no assembly here)
                A_ABM_full = BSH.abiotic_matrix(pool_ABM, grid; seed=seed_A + k)
                A_MAB_full = BSH.abiotic_matrix(pool_MAB, grid; seed=seed_A + k)
                A_BAM_full = BSH.abiotic_matrix(pool_BAM, grid; seed=seed_A + k)

                A_ABM_mask = @view A_ABM_full[:, keep]
                A_MAB_mask = @view A_MAB_full[:, keep]
                A_BAM_mask = @view A_BAM_full[:, keep]

                # mean climate-suitable area over consumers, normalized by ORIGINAL area
                push!(vals_ABM, BSH.mean_Aonly_over_consumers_area0(A_ABM_mask, pool_ABM, grid.C; τA=τA))
                push!(vals_MAB, BSH.mean_Aonly_over_consumers_area0(A_MAB_mask, pool_MAB, grid.C; τA=τA))
                push!(vals_BAM, BSH.mean_Aonly_over_consumers_area0(A_BAM_mask, pool_BAM, grid.C; τA=τA))

            end

            push!(A_ABM, mean(vals_ABM))
            push!(A_MAB, mean(vals_MAB))
            push!(A_BAM, mean(vals_BAM))
        end
        out[g] = (loss=collect(loss_fracs), ABM=A_ABM, MAB=A_MAB, BAM=A_BAM)
    end
    return out
end

end # module
