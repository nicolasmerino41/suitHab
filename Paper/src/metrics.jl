module Metrics

using Statistics

export worst_geometry, rank_flip, did_curves, regime_summary
export dist_metrics
export ks_best_fstar, ks_distance

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
        rAM, rBAM = BSH.per_species_relative_loss(rng, pool, grid, pars; fstar=f, geometry=geometry, seed_A=seed_A)
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

end # module
