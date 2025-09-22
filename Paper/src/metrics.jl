module Metrics

using Statistics
using ..BSH
using ..Grids
using ..Metawebs
using ..HL
import ..BSH: BAMParams

export worst_geometry, rank_flip, did_curves, regime_summary
export dist_metrics
export ks_best_fstar, ks_distance
export abs_bsh_vs_loss
export mixture_response_at_f, mixture_elasticity_at_f, mixture_geom_delta
export bootstrap_mixture_maps

"Pick geometry with largest relative loss (most negative)."
function worst_geometry(rel_by_geom::Dict{Symbol,Any}; at_index::Int)
    # expect keys :random, :clustered, :front
    losses = Dict(g => rel_by_geom[g].relBAM[at_index] for g in keys(rel_by_geom))
    # more negative = worse
    gmin = first(keys(losses))
    for (g,v) in losses
        if v < losses[gmin]; gmin = g; end
    end
    gmin
end

"Rank flip between AM and BAM at index."
function rank_flip(rel_by_geom::Dict{Symbol,Any}; at_index::Int)
    lossesAM = Dict(g => rel_by_geom[g].relAM[at_index] for g in keys(rel_by_geom))
    lossesB  = Dict(g => rel_by_geom[g].relBAM[at_index] for g in keys(rel_by_geom))
    argminAM = first(keys(lossesAM))
    for (g,v) in lossesAM; if v < lossesAM[argminAM]; argminAM=g; end; end
    argminB = first(keys(lossesB))
    for (g,v) in lossesB; if v < lossesB[argminB]; argminB=g; end; end
    (argminAM, argminB, argminAM != argminB)
end

"Compute DiD curves (BAM - AM) per geometry."
function did_curves(rel_by_geom::Dict{Symbol,Any})
    out = Dict{Symbol,Any}()
    for (g,r) in rel_by_geom
        out[g] = (; x=r.x, did = r.relBAM .- r.relAM)
    end
    out
end

"Small regime summary at f*: DiD mean over geometries, and flip?"
function regime_summary(rel_by_geom::Dict{Symbol,Any}; at_index::Int)
    di = [rel_by_geom[g].relBAM[at_index] - rel_by_geom[g].relAM[at_index] for g in keys(rel_by_geom)]
    _,_,flip = rank_flip(rel_by_geom; at_index)
    (; meanDiD=mean(di), flip=flip)
end

"KS distance, tail exceedance, and Gini for per-species relative losses."
function dist_metrics(relAM::Vector{Float64}, relBAM::Vector{Float64};
        tail_cut::Float64=0.8)
    # convert to positive losses x = -Δ/BSH0 in [0,1]
    xA = -relAM; xB = -relBAM
    # KS distance
    sA, sB = sort(xA), sort(xB)
    nA, nB = length(sA), length(sB)
    i=j=1; d=0.0
    while i≤nA || j≤nB
        vA = i≤nA ? sA[i] : Inf
        vB = j≤nB ? sB[j] : Inf
        v = min(vA,vB)
        FA = i>nA ? 1.0 : (i-1)/nA
        FB = j>nB ? 1.0 : (j-1)/nB
        d = max(d, abs(FA-FB))
        (vA==v) && (i+=1)
        (vB==v) && (j+=1)
    end
    # tail exceedance
    pA = mean(xA .> tail_cut)
    pB = mean(xB .> tail_cut)
    # Gini
    function gini(x)
        y = sort(x); n=length(y)
        n==0 && return 0.0
        2*sum((1:n) .* y)/(n*sum(y) + eps()) - (n+1)/n
    end
    gA, gB = gini(xA), gini(xB)
    (; KS=d, tailA=pA, tailB=pB, giniA=gA, giniB=gB, tailDiff=pB-pA, giniDiff=gB-gA)
end

"Two-sample KS distance between empirical CDFs of x and y on [0,1]."
function ks_distance(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    xs = sort!(copy(x)); ys = sort!(copy(y))
    grid = sort!(unique!(vcat(xs, ys)))
    Fx = [searchsortedlast(xs, t)/length(xs) for t in grid]
    Fy = [searchsortedlast(ys, t)/length(ys) for t in grid]
    maximum(abs.(Fx .- Fy))
end

"Pick f* that maximizes KS between AM and BAM for a chosen geometry."
function ks_best_fstar(; rng, pool, grid, pars, loss_fracs, geometry::Symbol=:front, seed_A::Int=1)
    best_f, best_KS = first(loss_fracs), -Inf
    for f in loss_fracs
        rAM, rBAM = per_species_relative_loss(rng, pool, grid, pars; fstar=f, geometry=geometry, seed_A=seed_A)
        # use positive-loss convention L = -Δ/BSH0 ∈ [0,1]
        LAM  = clamp.(-rAM, 0.0, 1.0)
        LBAM = clamp.(-rBAM, 0.0, 1.0)
        KS = ks_distance(LAM, LBAM)
        if KS > best_KS
            best_KS = KS; best_f = f
        end
    end
    (; f=best_f, KS=best_KS)
end

"""
abs_bsh_vs_loss(; rng, pool, grid, pars, loss_fracs, seed_A=1, A_fn=nothing, geoms=(:random,:clustered,:front))

Returns Dict{Symbol,NamedTuple} with fields:
  :loss  -> Vector{Float64}
  :AM    -> Vector{Float64}  (mean over consumers, normalized by ORIGINAL area)
  :BAM   -> Vector{Float64}  (same)
Uses the SAME A-construction path you use elsewhere (seed_A and optional A_fn).
"""
function abs_bsh_vs_loss(; rng, pool, grid, pars::BAMParams,
    loss_fracs::AbstractVector{<:Real}, seed_A::Int=1, A_fn=nothing,
    geoms::Tuple=(:random,:clustered,:front))

    # build A once, exactly as your BSH routines do
    A = isnothing(A_fn) ? abiotic_matrix(pool, grid; niche_width=0.12, seed=seed_A) :
                          A_fn(pool, grid; seed=seed_A)

    out = Dict{Symbol,NamedTuple}()

    for g in geoms
        am = Float64[]; bm = Float64[]
        for f in loss_fracs
            keepfrac = 1 - f
            keep = g === :random    ? HL.random_mask(rng, grid.C, keepfrac) :
                   g === :clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                      HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)

            bAM = mean_BSH_over_consumers_area0(rng, pool, grid, pars;
                    A=A, keepmask=keep, mode=:AM)
            bBM = mean_BSH_over_consumers_area0(rng, pool, grid, pars;
                    A=A, keepmask=keep, mode=:BAM)
            push!(am, bAM); push!(bm, bBM)
        end
        out[g] = (loss=collect(loss_fracs), AM=am, BAM=bm)
    end
    return out
end

# --- src/metrics.jl additions ---
export mixture_response_at_f, mixture_elasticity_at_f, mixture_geom_delta

"Retained suitable area at f* for a given A/B/M mixture (weights sum to 1)."
function mixture_response_at_f(; rng, pool, grid, pars::BSH.BAMParams,
    wA::Float64, wB::Float64, wM::Float64, fstar::Float64, geometry::Symbol,
    seed_A::Int=1, A_fn=abiotic_matrix, agg::Symbol=:mean, kreq::Int=1)

    @assert isapprox(wA + wB + wM, 1.0; atol=1e-8)
    keepfrac = 1 - fstar
    keep = geometry === :random    ? HL.random_mask(rng, grid.C, keepfrac) :
           geometry === :clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                     HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)

    # Build A once; reuse for A-only, AM, BAM
    A = A_fn(pool, grid; seed=seed_A)

    # A-only (climate pass only; average over consumers vs ORIGINAL area)
    A_mask = @view A[:, keep]
    yA = BSH.mean_Aonly_over_consumers_area0(A_mask, pool, grid.C)

    # AM (abiotic + movement)
    yM = BSH.mean_BSH_over_consumers_area0(rng, pool, grid, pars; A=A, keepmask=keep, mode=:AM)

    # BAM (abiotic + movement + biotic; pass agg/kreq to control stringency)
    yB = BSH.mean_BSH_over_consumers_area0(rng, pool, grid, pars; A=A, keepmask=keep, mode=:BAM)

    return wA*yA + wM*yM + wB*yB
end

"Finite-difference elasticity at f* for a given mixture and geometry."
function mixture_elasticity_at_f(; rng, pool, grid, pars, wA, wB, wM, fstar, geometry,
    δ::Float64=0.02, kwargs...)

    f1 = clamp(fstar - δ/2, 0.0, 0.95)
    f2 = clamp(fstar + δ/2, 0.0, 0.95)
    y1 = mixture_response_at_f(; rng, pool, grid, pars, wA, wB, wM, fstar=f1, geometry, kwargs...)
    y2 = mixture_response_at_f(; rng, pool, grid, pars, wA, wB, wM, fstar=f2, geometry, kwargs...)
    return (y2 - y1) / (f2 - f1)
end

"Geometry contrast at f*: Δ = y_random − y_other for a mixture."
function mixture_geom_delta(; rng, pool, grid, pars, wA, wB, wM, fstar, other::Symbol, kwargs...)
    yR = mixture_response_at_f(; rng, pool, grid, pars, wA, wB, wM, fstar, geometry=:random, kwargs...)
    yO = mixture_response_at_f(; rng, pool, grid, pars, wA, wB, wM, fstar, geometry=other,  kwargs...)
    return yR - yO
end

"""
bootstrap_mixture_maps(; rng, pool, grid, pars, weights, fstar, A_fn, agg, kreq,
                       geometries=(:random,:clustered,:front),
                       A_seeds=1:24, pool_seeds=nothing)

Returns:
  elastic_mean[g], elastic_lo[g], elastic_hi[g] :: Dict{Symbol, Vector{Float64}}
  dRF_mean, dRF_lo, dRF_hi :: Vector{Float64}
  dRC_mean, dRC_lo, dRC_hi :: Vector{Float64}
(‘lo/hi’ are p10/p90 across seeds.)
"""
function bootstrap_mixture_maps(; rng, pool, grid, pars,
        weights::Vector{NTuple{3,Float64}},
        fstar::Float64, A_fn, agg::Symbol, kreq::Int,
        geometries::Tuple=(:random,:clustered,:front),
        A_seeds = 1:24, pool_seeds = nothing)

    # Helper to get pool per replicate (fixed pool if pool_seeds = nothing)
    function pool_for(rep)
        if pool_seeds === nothing
            return pool
        else
            seed = pool_seeds[ ((rep-1) % length(pool_seeds)) + 1 ]
            return Metawebs.build_metaweb_archetype(MersenneTwister(seed);
                                                    S=pool.S, basal_frac=mean(pool.basal), archetype=:mid)
        end
    end

    # Containers
    elast_vals = Dict(g => [Float64[] for _ in eachindex(weights)] for g in geometries)
    dRF_vals   = [Float64[] for _ in eachindex(weights)]
    dRC_vals   = [Float64[] for _ in eachindex(weights)]

    nrep = length(A_seeds)

    for (rep, seedA) in enumerate(A_seeds)
        p = pool_for(rep)
        for (k, (a,b,m)) in enumerate(weights)
            # elasticities per geometry
            for g in geometries
                e = mixture_elasticity_at_f(; rng, pool=p, grid, pars,
                        wA=a, wB=b, wM=m, fstar, geometry=g,
                        A_fn=A_fn, agg=agg, kreq=kreq, seed_A=seedA)
                push!(elast_vals[g][k], e)
            end
            # deltas
            rf = mixture_geom_delta(; rng, pool=p, grid, pars, wA=a, wB=b, wM=m,
                                    fstar, other=:front,     A_fn=A_fn, agg=agg, kreq=kreq, seed_A=seedA)
            rc = mixture_geom_delta(; rng, pool=p, grid, pars, wA=a, wB=b, wM=m,
                                    fstar, other=:clustered, A_fn=A_fn, agg=agg, kreq=kreq, seed_A=seedA)
            push!(dRF_vals[k], rf);  push!(dRC_vals[k], rc)
        end
    end

    # Summarize
    function summarize(vs)
        μ  = [mean(v) for v in vs]
        lo = [quantile(v,0.10) for v in vs]
        hi = [quantile(v,0.90) for v in vs]
        μ, lo, hi
    end

    elastic_mean = Dict{Symbol,Vector{Float64}}()
    elastic_lo   = Dict{Symbol,Vector{Float64}}()
    elastic_hi   = Dict{Symbol,Vector{Float64}}()
    for g in geometries
        μ, lo, hi = summarize(elast_vals[g])
        elastic_mean[g] = μ; elastic_lo[g] = lo; elastic_hi[g] = hi
    end
    dRF_mean, dRF_lo, dRF_hi = summarize(dRF_vals)
    dRC_mean, dRC_lo, dRC_hi = summarize(dRC_vals)

    return (; elastic_mean, elastic_lo, elastic_hi,
            dRF_mean, dRF_lo, dRF_hi,
            dRC_mean, dRC_lo, dRC_hi)
end

end # module
