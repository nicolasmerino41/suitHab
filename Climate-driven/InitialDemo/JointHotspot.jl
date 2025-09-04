# Joint (climate × prey-support) hotspot removal
# Z0: climate pass at τ (BitMatrix S×C), pool: SpeciesPool, keep∈(0,1]
function joint_hotspot_mask(Z0::BitMatrix, pool::SpeciesPool, keep::Float64)
    S, C = size(Z0)
    score = zeros(Float64, C)
    @inbounds for c in 1:C
        for s in 1:S
            Z0[s,c] || continue
            pool.basal[s] && continue
            prey = pool.E[s]
            isempty(prey) && continue
            # fraction of s's prey with climate-suitability in the same cell
            supp = count(q -> Z0[q,c], prey) / length(prey)
            score[c] += supp
        end
    end
    nkeep = clamp(round(Int, keep*C), 0, C)
    idx = sortperm(score; rev=true)   # highest score first (worst to remove)
    keepmask = falses(C)
    keepmask[idx[1:nkeep]] .= true
    return keepmask
end

"""
apply_viability(p; mode=:threshold, ρ=0.6, h=0.2, k=6)

Map prey-support p (fraction of suitable cells that have ≥1 prey present) to a
viability score φ(p) in [0,1]. NaNs in p are propagated as NaN.

modes:
  :threshold  -> hard cutoff:  φ = 1(p ≥ ρ)
  :saturating -> logistic around ρ with width h
  :hill       -> Hill/sigmoid: φ = p^k / (p^k + ρ^k)

Arguments
  p    :: AbstractVector{<:Real}   # typically from bsh1_cond_per_species
  ρ    :: Real    # location of the cutoff/half-saturation (0..1)
  h    :: Real    # logistic width (only for :saturating), >0
  k    :: Real    # Hill exponent (only for :hill), >0
"""
function apply_viability(p; mode::Symbol=:threshold, ρ::Real=0.6, h::Real=0.2, k::Real=6)
    φ = similar(p, Float64)
    if mode === :threshold
        @inbounds for i in eachindex(p)
            φ[i] = isnan(p[i]) ? NaN : (p[i] ≥ ρ ? 1.0 : 0.0)
        end
    elseif mode === :saturating
        hw = max(float(h), 1e-6)
        @inbounds for i in eachindex(p)
            φ[i] = isnan(p[i]) ? NaN : 1.0 / (1.0 + exp(-(p[i] - ρ)/hw))
        end
    elseif mode === :hill
        kk = max(float(k), 1e-6)
        ρk = float(ρ)^kk + eps()
        @inbounds for i in eachindex(p)
            if isnan(p[i])
                φ[i] = NaN
            else
                pk = float(p[i])^kk
                φ[i] = pk / (pk + ρk)
            end
        end
    else
        error("apply_viability: unknown mode $(mode). Use :threshold, :saturating, or :hill.")
    end
    return φ
end

# mean B across consumers for a given Z and viability
function mean_B_from_Z(Z::BitMatrix, pool::SpeciesPool; viability=(;mode=:threshold, ρ=0.6, h=0.2))
    P = assemble(Z, pool)
    A = climate_area_per_species(Z)
    p = bsh1_cond_per_species(P, Z, pool)
    φ = apply_viability(p; viability...)
    B = A .* φ
    cons = .!pool.basal
    return mean(B[cons][.!isnan.(B[cons])])
end

# One-shot ΔB (relative to intact) for a given keep mask
function deltaB_for_mask(pool::SpeciesPool, grid::Grid; τ=0.55, keepmask, viability=(;mode=:threshold, ρ=0.6, h=0.2))
    Z0 = climate_pass(pool, grid; τ=τ)
    B0 = mean_B_from_Z(Z0, pool; viability=viability)
    Z1 = apply_mask(Z0, keepmask)
    B1 = mean_B_from_Z(Z1, pool; viability=viability)
    return B1 - B0
end

# Excess damage curves vs loss
function excess_damage_curves(grid::Grid, pool::SpeciesPool;
        τ::Float64=0.55, losses=0.0:0.05:0.8, nseeds_cluster::Int=1, seed::Int=123,
        viability=(;mode=:threshold, ρ=0.6, h=0.2))

    Z0 = climate_pass(pool, grid; τ=τ)

    ΔBr = Float64[]
    ΔBc = Float64[]
    ΔBh = Float64[]  # joint climate×prey hotspot

    for f in losses
        keep = 1 - f
        km_r = random_mask(grid.C, keep; seed=seed)
        km_c = clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=seed+1)
        km_h = joint_hotspot_mask(Z0, pool, keep)

        push!(ΔBr, deltaB_for_mask(pool, grid; τ=τ, keepmask=km_r, viability=viability))
        push!(ΔBc, deltaB_for_mask(pool, grid; τ=τ, keepmask=km_c, viability=viability))
        push!(ΔBh, deltaB_for_mask(pool, grid; τ=τ, keepmask=km_h, viability=viability))
    end

    excess_clustered = ΔBc .- ΔBr
    excess_hotspot   = ΔBh .- ΔBr

    return (; losses=collect(losses), ΔBr, ΔBc, ΔBh, excess_clustered, excess_hotspot)
end

# --- Run & plot ---
# (Use the same pool you used for the diagnostic; τ and viability match your threshold=0.6.)
curves = excess_damage_curves(grid, pool; τ=0.55, losses=0.0:0.05:0.8,
                              nseeds_cluster=1,
                              viability=(;mode=:threshold, ρ=0.6, h=0.2))

begin
    fig = Figure(; size=(820,360))
    ax  = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Excess ΔBSH (vs random)",
               title="Excess damage curve")
    lines!(ax, curves.losses, curves.excess_clustered, label="Clustered − Random")
    lines!(ax, curves.losses, curves.excess_hotspot,   label="Joint hotspot − Random")
    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax, position=:lt)
    display(fig)
end
