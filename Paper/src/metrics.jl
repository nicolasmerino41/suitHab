module Metrics

using Statistics

export worst_geometry, rank_flip, did_curves, regime_summary

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

end # module
