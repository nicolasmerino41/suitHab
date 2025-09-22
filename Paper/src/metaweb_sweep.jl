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
export run_realities_suitable_area
export run_blend_realities_suitable_area, simplex_grid
export blend_geometry_delta_at_fstar, blend_elasticity_at_fstar

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

# -- helper: mean (over consumers) of occupied cells in keep-mask / ORIGINAL area
_mean_area0_from_P(P::BitMatrix, pool, grid, keep::BitVector) = begin
    cons = findall(!, pool.basal)
    isempty(cons) && return 0.0
    vals = Float64[]
    @inbounds for s in cons
        push!(vals, sum(@view P[s, keep]) / grid.C)
    end
    return mean(vals)
end

# -- build a "reality" = baseline occupancy map P with dominance tilted to A, M, or B
function _build_reality_P(; rng::AbstractRNG, grid,
    S::Int, basal_frac::Float64, seed_pool::Int, seed_A::Int,
    which::Symbol,  # :ABM | :MAB | :BAM
    pars_ABM::BSH.BAMParams, pars_MAB::BSH.BAMParams, pars_BAM::BSH.BAMParams,
    agg::Symbol=:mean, kreq::Int=1
)
    # pool choice can be neutral or you can vary redundancy; I keep it neutral (= :mid)
    pool = Metawebs.build_metaweb_archetype(rng; S=S, basal_frac=basal_frac, archetype=:mid)
    A    = BSH.abiotic_matrix(pool, grid; seed=seed_A)
    keep_all = trues(grid.C)

    if which === :ABM
        # abiotic-dominated: movement OFF or weak, stricter τA, no prey gate
        P, _ = BSH.assemble_AM(pool, grid, A, keep_all; pars=pars_ABM)
        return (pool=pool, P=P)
    elseif which === :MAB
        # movement-dominated: strong component gate (large T), otherwise AM
        P, _ = BSH.assemble_AM(pool, grid, A, keep_all; pars=pars_MAB)
        return (pool=pool, P=P)
    elseif which === :BAM
        # biotic-dominated: strong prey sufficiency; movement optional
        out = BSH.assemble_BAM(pool, grid, A, keep_all; pars=pars_BAM, agg=agg, kreq=kreq)
        return (pool=pool, P=out.P)
    else
        error("which must be :ABM, :MAB, or :BAM")
    end
end

"""
run_realities_suitable_area(; rng, grid, loss_fracs, S, basal_frac,
    Npools=12, seed_pool0=1, seed_A=1,
    geoms=(:random,:clustered,:front),
    pars_ABM=BSH.BAMParams(τA=0.55, movement=:off),
    pars_MAB=BSH.BAMParams(τA=0.50, movement=:component, T=12),
    pars_BAM=BSH.BAMParams(τA=0.50, τB=0.50, movement=:off),
    agg=:kofn, kreq=2)

For each geometry and loss fraction, builds Npools independent *baseline* occupancy
maps that represent ABM / MAB / BAM “worlds”, then **clips** them with the mask
(no re-assembly), and returns the mean (over consumers) suitable area / original area.

Returns Dict{Symbol,NamedTuple} mapping geometry to:
  (loss, ABM::Vector, MAB::Vector, BAM::Vector)
"""
function run_realities_suitable_area(; rng::AbstractRNG, grid,
    loss_fracs=0.2:0.1:0.8, S::Int=200, basal_frac::Float64=0.25,
    Npools::Int=12, seed_pool0::Int=1, seed_A::Int=1,
    geoms::Tuple=(:random,:clustered,:front),
    pars_ABM=BAMParams(; τA=0.55, movement=:off),
    pars_MAB=BAMParams(; τA=0.50, movement=:component, T=12),
    pars_BAM=BAMParams(; τA=0.50, τB=0.50, movement=:off),
    agg::Symbol=:kofn, kreq::Int=2
)
    out = Dict{Symbol,NamedTuple}()

    for g in geoms
        yABM = Float64[]; yMAB = Float64[]; yBAM = Float64[]
        for f in loss_fracs
            keepfrac = 1 - f
            keep = g === :random    ? HL.random_mask(rng, grid.C, keepfrac) :
                   g === :clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                      HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)

            valsA = Float64[]; valsM = Float64[]; valsB = Float64[]

            for k in 0:Npools-1
                rseed = seed_pool0 + k
                # ABM world
                AB = _build_reality_P(; rng=MersenneTwister(hash((:ABM, rseed))),
                    grid, S, basal_frac, seed_pool=rseed, seed_A=seed_A+k,
                    which=:ABM, pars_ABM, pars_MAB, pars_BAM, agg, kreq)
                push!(valsA, _mean_area0_from_P(AB.P, AB.pool, grid, keep))

                # MAB world
                MB = _build_reality_P(; rng=MersenneTwister(hash((:MAB, rseed))),
                    grid, S, basal_frac, seed_pool=rseed, seed_A=seed_A+k,
                    which=:MAB, pars_ABM, pars_MAB, pars_BAM, agg, kreq)
                push!(valsM, _mean_area0_from_P(MB.P, MB.pool, grid, keep))

                # BAM world
                BB = _build_reality_P(; rng=MersenneTwister(hash((:BAM, rseed))),
                    grid, S, basal_frac, seed_pool=rseed, seed_A=seed_A+k,
                    which=:BAM, pars_ABM, pars_MAB, pars_BAM, agg, kreq)
                push!(valsB, _mean_area0_from_P(BB.P, BB.pool, grid, keep))
            end

            push!(yABM, mean(valsA)); push!(yMAB, mean(valsM)); push!(yBAM, mean(valsB))
        end
        out[g] = (loss=collect(loss_fracs), ABM=yABM, MAB=yMAB, BAM=yBAM)
    end
    return out
end

"Quick grid on the A/M/B simplex (step in {0,0.1,0.2,...,1}, sum≈1)."
function simplex_grid(step::Float64=0.1)
    W = Tuple{Float64,Float64,Float64}[]
    vals = collect(0.0:step:1.0)
    for a in vals, m in vals
        b = 1.0 - a - m
        if b ≥ -1e-9 && b ≤ 1.0 + 1e-9
            push!(W, (round(a;digits=3), round(m;digits=3), round(max(0.0,min(1.0,b));digits=3)))
        end
    end
    unique(W)
end

"Build the three baseline maps once."
function _baseline_maps(; rng, grid, S::Int, basal_frac::Float64, seed_pool::Int, seed_A::Int,
    pars_ABM::BAMParams, pars_MAB::BAMParams, pars_BAM::BAMParams,
    agg::Symbol=:kofn, kreq::Int=2)

    pool = Metawebs.build_metaweb_archetype(rng; S=S, basal_frac=basal_frac, archetype=:mid)
    A    = BSH.abiotic_matrix(pool, grid; seed=seed_A)
    keep_all = trues(grid.C)

    P_A, _ = BSH.assemble_AM(pool, grid, A, keep_all; pars=pars_ABM)
    P_M, _ = BSH.assemble_AM(pool, grid, A, keep_all; pars=pars_MAB)
    outB   = BSH.assemble_BAM(pool, grid, A, keep_all; pars=pars_BAM, agg=agg, kreq=kreq)
    P_B    = outB.P
    (; pool, P_A, P_M, P_B)
end

"Blend species rows from P_A, P_M, P_B according to weights (wA,wM,wB)."
function _blend_P(pool, P_A::BitMatrix, P_M::BitMatrix, P_B::BitMatrix,
                  wA::Float64, wM::Float64, rng::AbstractRNG)
    S, C = size(P_A)
    P = falses(S, C)
    cons = findall(!, pool.basal)
    n = length(cons)
    nA = round(Int, wA * n); nM = round(Int, wM * n)
    nB = max(0, n - nA - nM)
    idx = copy(cons); shuffle!(rng, idx)
    SA, SM, SB = idx[1:nA], idx[nA+1:nA+nM], idx[nA+nM+1:end]
    P[SA, :] .= P_A[SA, :]
    P[SM, :] .= P_M[SM, :]
    P[SB, :] .= P_B[SB, :]
    # keep basal identical across blends (use P_A rows, which equal P_M in AM)
    bas = findall(pool.basal)
    P[bas, :] .= P_A[bas, :]
    return P
end

"""
run_blend_realities_suitable_area(; rng, grid, loss_fracs, S, basal_frac,
    weights::Vector{Tuple{Float64,Float64,Float64}},  # list of (wA,wM,wB)
    geoms=(:random,:clustered,:front), Nrep=8, seed0=1,
    pars_ABM=BSH.BAMParams(τA=0.55, movement=:off),
    pars_MAB=BSH.BAMParams(τA=0.50, movement=:component, T=14),
    pars_BAM=BSH.BAMParams(τA=0.50, τB=0.50, movement=:off),
    agg=:kofn, kreq=2)

Returns Dict{Symbol,Dict{Tuple=>NamedTuple}} so that
out[:random][(wA,wM,wB)] = (loss, y)
where y is the curve of suitable area (mean over consumers / original area).
"""
function run_blend_realities_suitable_area(; rng::AbstractRNG, grid,
    loss_fracs=0.2:0.1:0.8, S::Int=200, basal_frac::Float64=0.25,
    weights::Vector{Tuple{Float64,Float64,Float64}} = simplex_grid(0.25),
    geoms::Tuple=(:random,:clustered,:front), Nrep::Int=8, seed0::Int=1,
    pars_ABM::BSH.BAMParams=BSH.BAMParams(; τA=0.55, movement=:off),
    pars_MAB::BSH.BAMParams=BSH.BAMParams(; τA=0.50, movement=:component, T=14),
    pars_BAM::BSH.BAMParams=BSH.BAMParams(; τA=0.50, τB=0.50, movement=:off),
    agg::Symbol=:kofn, kreq::Int=2
)
    out = Dict{Symbol,Dict{Tuple{Float64,Float64,Float64},NamedTuple}}(
        g => Dict{Tuple{Float64,Float64,Float64},NamedTuple}() for g in geoms)

    for g in geoms
        for w in weights
            wA, wM, wB = w
            y = Float64[]
            for f in loss_fracs
                keepfrac = 1 - f
                keep = g === :random    ? HL.random_mask(rng, grid.C, keepfrac) :
                       g === :clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                          HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)
                repvals = Float64[]
                for r in 0:Nrep-1
                    base = _baseline_maps(
                        ; rng=MersenneTwister(hash((:base, g, w, r, seed0))),
                          grid, S, basal_frac, seed_pool=seed0+r, seed_A=seed0+r,
                          pars_ABM, pars_MAB, pars_BAM, agg, kreq)
                    Pmix = _blend_P(base.pool, base.P_A, base.P_M, base.P_B, wA, wM,
                                    MersenneTwister(hash((:assign, g, w, r, seed0))))
                    push!(repvals, _mean_area0_from_P(Pmix, base.pool, grid, keep))
                end
                push!(y, mean(repvals))
            end
            out[g][w] = (loss=collect(loss_fracs), y=y)
        end
    end
    return out
end

"""
blend_geometry_delta_at_fstar(blend, fstar; geomA=:random, geomB=:front)
Given the `blend` dict from run_blend_realities_suitable_area, return
a Vector of (wA,wM,wB, delta) at f* where delta = y_geomA - y_geomB.
"""
function blend_geometry_delta_at_fstar(blend::Dict{Symbol,<:Any}, fstar::Float64;
    geomA::Symbol=:random, geomB::Symbol=:front)

    loss = first(values(blend[geomA])).loss
    k = argmin(abs.(loss .- fstar))
    out = Tuple{Float64,Float64,Float64,Float64}[]
    for w in keys(blend[geomA])
        ya = blend[geomA][w].y[k]
        yb = blend[geomB][w].y[k]
        push!(out, (w[1], w[2], 1.0 - w[1] - w[2], ya - yb))
    end
    out
end

"""
Numerical derivative dy/df at f* for each weight triple and geometry.
Returns Dict{Symbol,Vector{(wA,wM,wB,slope)}}.
"""
function blend_elasticity_at_fstar(blend::Dict{Symbol,<:Any}, fstar::Float64)
    res = Dict{Symbol, Vector{Tuple{Float64,Float64,Float64,Float64}}}()
    loss = first(values(first(values(blend)))).loss
    k = argmin(abs.(loss .- fstar))
    k = clamp(k, 2, length(loss)-1)
    for g in keys(blend)
        acc = Tuple{Float64,Float64,Float64,Float64}[]
        for w in keys(blend[g])
            y = blend[g][w].y
            slope = (y[k+1]-y[k-1]) / (loss[k+1]-loss[k-1])
            push!(acc, (w[1], w[2], 1.0 - w[1] - w[2], slope))
        end
        res[g] = acc
    end
    res
end

end # module
