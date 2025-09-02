function assemble_cascade(Z::BitMatrix, pool::SpeciesPool; eta::Union{Int,Vector{Int}}=1)
    P = assemble(Z, pool)
    S, C = size(P)
    changed = true
    while changed
        changed = false
        occ = vec(sum(P; dims=2))                # regional occupancy per sp.
        extinct = if eta isa Int
            findall(<(eta), occ)
        else
            [s for s in 1:S if occ[s] < eta[s]]
        end
        if !isempty(extinct)
            for s in extinct
                if any(@view P[s,:])
                    P[s,:] .= false
                    changed = true
                end
            end
        end
        # local prey support
        @inbounds for c in 1:C, s in 1:S
            if P[s,c] && !pool.basal[s]
                has_pre = false
                for q in pool.E[s]
                    if P[q,c]; has_pre = true; break; end
                end
                if !has_pre
                    P[s,c] = false
                    changed = true
                end
            end
        end
    end
    return P
end

# mean BSH over consumers on the FULL grid denominator
mean_BSH_consumers(P::BitMatrix, Z::BitMatrix, pool::SpeciesPool) =
    mean(filter(!isnan, bsh1_per_species(P, Z, pool)[.!pool.basal]))

# Build per-species regional thresholds from baseline
function baseline_thresholds(Zfull::BitMatrix, pool::SpeciesPool; α::Float64=0.6, ηmin::Int=1)
    P0   = assemble(Zfull, pool)
    occ0 = vec(sum(P0; dims=2))                     # baseline occupancy
    η    = similar(occ0, Int)
    @inbounds for s in eachindex(occ0)
        η[s] = occ0[s] == 0 ? 0 : max(ηmin, ceil(Int, α * occ0[s]))
    end
    return η, P0
end

function decomp_clim_int_syn_relative(pool::SpeciesPool, grid::Grid;
                                      τ::Float64, keepmask::BitVector,
                                      α::Float64=0.6, ηmin::Int=1)
    # full grid climate and species-specific thresholds
    Zfull = climate_pass(pool, grid; τ=τ)
    η, P00 = baseline_thresholds(Zfull, pool; α=α, ηmin=ηmin)

    # states: (mask M, cascade R)
    # B00: no loss, no cascade
    Z00 = Zfull
    P00 = P00
    B00 = mean_BSH_consumers(P00, Z00, pool)

    # helper to apply mask
    applyM(Z, M) = (Z .& reshape(M,1,:))

    # B10: loss only (mask on, cascade off)
    Z10 = applyM(Zfull, keepmask)
    P10 = assemble(Z10, pool)
    B10 = mean_BSH_consumers(P10, Z10, pool)

    # B01: cascade only (no loss, cascade on) — should ≈ B00
    P01 = assemble_cascade(Z00, pool; eta=η)
    B01 = mean_BSH_consumers(P01, Z00, pool)

    # B11: loss + cascade
    Z11 = Z10
    P11 = assemble_cascade(Z11, pool; eta=η)
    B11 = mean_BSH_consumers(P11, Z11, pool)

    clim = B10 - B00
    inter = B01 - B00
    syn  = (B11 - B10) - (B01 - B00)   # non-additive remainder
    total = B11 - B00
    return (; total, clim, inter, syn)
end

grid = make_grid(60,60; seed=11)
S = 200; basal_frac = 0.35

# Make prey clumpy & diets narrow (low redundancy, high synchrony)
pool = build_pool(S; basal_frac=basal_frac, seed=1,
                  sigma=0.18, density=0.10, pmax=0.65,
                  niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.03,
                  b0_basal=0.08, bspread_basal=0.02, b0_cons=0.10, bspread_cons=0.03)

τ = 0.62                    # stricter climate suitability
keep = 0.5                  # 50% loss
power = 3.0
kmR = random_mask(grid.C, keep; seed=101)
kmC = clustered_mask(grid, keep; nseeds=1, seed=101)
# kmH = consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep, power=power, what=:cons, seed=101)
kmH = consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep,
                            power=3.0, seed=101)


# tougher persistence: must keep ≥ α of baseline cells (try 0.6–0.8)
α = 0.7; ηmin = 2

dR = decomp_clim_int_syn_relative(pool, grid; τ=τ, keepmask=kmR, α=α, ηmin=ηmin)
dC = decomp_clim_int_syn_relative(pool, grid; τ=τ, keepmask=kmC, α=α, ηmin=ηmin)
dH = decomp_clim_int_syn_relative(pool, grid; τ=τ, keepmask=kmH, α=α, ηmin=ηmin)

# stack & draw (your polygon bars routine)
Y = vcat([dR.clim dR.inter dR.syn],
         [dC.clim dC.inter dC.syn],
         [dH.clim dH.inter dH.syn])

begin
    fig = Figure(; size=(900,340))
    ax = Axis(fig[1,1], ylabel="ΔBSH (50% loss)", xlabel="Scenario")
    x = 1:3
    stackedbars_poly!(ax, collect(x), Y; width=0.8, colors=[:dodgerblue,:orange,:forestgreen])
    ax.xticks = (x, ["Random","Clustered","Hotspot"])
    ax.xticklabelrotation[] = π/12
    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax, [PolyElement(color=:dodgerblue),PolyElement(color=:orange),PolyElement(color=:forestgreen)],
                ["Climate","Interaction","Synergy"]; position=:lt)
    display(fig)
end

# how many species fail the baseline-relative threshold after loss?
function failed_counts(pool, grid; τ, keepmask, α=0.7, ηmin=2)
    Zfull = climate_pass(pool, grid; τ=τ)
    η, _ = baseline_thresholds(Zfull, pool; α=α, ηmin=ηmin)
    Z = Zfull .& reshape(keepmask,1,:)
    P = assemble(Z, pool)
    occ = vec(sum(P; dims=2))
    fails = sum(occ .< η)
    fails_consumers = sum((occ .< η) .& .!pool.basal)
    return (; fails_total=fails, fails_consumers=fails_consumers)
end

failed_counts(pool, grid; τ=τ, keepmask=kmH, α=α, ηmin=ηmin)
