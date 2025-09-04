"Return integer trophic level per species on the metaweb (basal=1)."
function trophic_level_metaweb(pool::SpeciesPool)
    S = pool.S
    ord = sortperm(pool.masses)      # increasing mass (acyclic)
    tl  = ones(Int, S)               # basal default = 1
    for s in ord                      # from small → large
        prey = pool.E[s]
        isempty(prey) && continue
        tl[s] = 1 + maximum(tl[prey]) # longest-path style level
    end
    return tl
end

"""
decomp_by_TL(pool, grid; τ, keepmask, include_basal=false)

Returns:
  tl_vals      :: Vector{Int}          unique trophic levels present
  n_per_tl     :: Vector{Int}          #species per TL used
  mean_dAcl    :: Vector{Float64}
  mean_dInt    :: Vector{Float64}
  mean_dSyn    :: Vector{Float64}
  mean_dB      :: Vector{Float64}
(All on area basis; averages are over species in each TL.)
"""
function decomp_by_TL(pool::SpeciesPool, grid::Grid;
                      τ::Float64=0.5, keepmask::BitVector,
                      include_basal::Bool=true)
    Cfull = grid.C

    # before/after components
    Z0 = climate_pass(pool, grid; τ=τ);  P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask);       P1 = assemble(Z1, pool)
    A0, p0, B0 = components_area(P0, Z0, pool, Cfull)
    A1, p1, B1 = components_area(P1, Z1, pool, Cfull)

    # species trophic level (metaweb)
    tl = trophic_level_metaweb(pool)

    # which species to include
    use = include_basal ? trues(pool.S) : .!pool.basal

    # per-species pieces
    dA_only   = similar(B0); dInt_only = similar(B0); dSyn = similar(B0); dB = similar(B0)
    @inbounds for s in 1:pool.S
        if use[s]
            d = decompose_midpoint(A0[s], p0[s], B0[s], A1[s], p1[s], B1[s])
            dB[s]       = d.dB
            dA_only[s]  = d.dA_only
            dInt_only[s]= d.dInt_only
            dSyn[s]     = d.synergy
        else
            dB[s] = dA_only[s] = dInt_only[s] = dSyn[s] = NaN
        end
    end

    # group by integer TL
    tls = sort(unique(tl[use]))
    mean_or_nan(v,i) = isempty(i) ? NaN : mean(@view v[i])

    mean_dAcl = Float64[]; mean_dInt = Float64[]; mean_dSyn = Float64[]; mean_dB = Float64[]; n_per = Int[]
    for L in tls
        idx = findall((tl .== L) .& use)
        push!(n_per, length(idx))
        push!(mean_dAcl, mean_or_nan(dA_only, idx))
        push!(mean_dInt, mean_or_nan(dInt_only, idx))
        push!(mean_dSyn, mean_or_nan(dSyn, idx))
        push!(mean_dB,   mean_or_nan(dB, idx))
    end
    return (tl_vals = tls, n_per_tl = n_per,
            mean_dAcl = mean_dAcl, mean_dInt = mean_dInt,
            mean_dSyn = mean_dSyn, mean_dB = mean_dB)
end

"Simple TL bar plot for one mask."
function plot_TL_partition(pool::SpeciesPool, grid::Grid; τ=0.55, keepmask,
                           title_str::String="ΔBSH by trophic level (area basis)")
    r = decomp_by_TL(pool, grid; τ=τ, keepmask=keepmask, include_basal=false)
    begin
        fig = Figure(; size=(860,360))
        ax  = Axis(fig[1,1], title=title_str, xlabel="Trophic level", ylabel="Contribution")
        xs  = r.tl_vals
        barplot!(ax, xs .- 0.30, r.mean_dAcl; width=0.25, color=:steelblue, label="Climate-only")
        barplot!(ax, xs,          r.mean_dInt; width=0.25, color=:orange,    label="Interaction-only")
        barplot!(ax, xs .+ 0.30,  r.mean_dSyn; width=0.25, color=:seagreen,  label="Synergy")
        hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
        axislegend(ax, position=:rb)
        display(fig)
    end
    return r
end
