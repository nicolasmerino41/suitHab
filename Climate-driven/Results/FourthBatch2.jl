# Local assembly you already have:
# assemble(Z::BitMatrix, pool::SpeciesPool) -> P::BitMatrix

"""
assemble_cascade(Z, pool; eta::Int)
Local assembly + regional persistence:
- Start from assemble(Z, pool)
- Repeatedly remove any species with regional occupancy < eta
  and prune consumers that lose all prey locally, until fixed point.
"""
function assemble_cascade(Z::BitMatrix, pool::SpeciesPool; eta::Int=1)
    P = assemble(Z, pool)
    S, C = size(P)
    changed = true
    while changed
        changed = false
        # regional occupancy (cells with presence) per species
        occ = vec(sum(P; dims=2))             # length S
        extinct = findall(<(eta), occ)        # occ < eta
        if !isempty(extinct)
            for s in extinct
                if any(@view P[s, :])
                    P[s, :] .= false
                    changed = true
                end
            end
        end
        # after removing regionally extinct spp, enforce local prey support
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

# mean BSH over consumers; Z is sized S×C over the FULL grid
mean_BSH_consumers(P::BitMatrix, Z::BitMatrix, pool::SpeciesPool) = begin
    b = bsh1_per_species(P, Z, pool)[.!pool.basal]   # your existing function
    mean(filter(!isnan, b))
end

"""
Evaluate mean consumer BSH under a given mask & cascade switch.
- keepmask :: BitVector length C (or `nothing` for no loss)
- eta_frac :: fraction of *kept* cells required for regional persistence (≥0)
"""
function eval_state(pool::SpeciesPool, grid::Grid; τ::Float64,
                    keepmask::Union{Nothing,BitVector}, cascade::Bool,
                    eta_frac::Float64 = 0.1)
    Zfull = climate_pass(pool, grid; τ=τ)                # S×C over full grid
    M = keepmask === nothing ? trues(grid.C) : keepmask  # length C
    # zero climate in removed cells but keep full width (fraction basis)
    Z = falses(size(Zfull));  Z .= Zfull .& reshape(M, 1, :)
    if cascade
        nkept = count(M)
        eta   = max(1, ceil(Int, eta_frac * nkept))
        P = assemble_cascade(Z, pool; eta=eta)
    else
        P = assemble(Z, pool)
    end
    return mean_BSH_consumers(P, Z, pool)
end

"""
Decompose ΔBSH into Climate, Interaction, and Synergy (non-additive remainder)
between baseline (no loss, no cascade) and post-loss (loss, cascade).
Returns (total, climate, interaction, synergy).
"""
function decomp_clim_int_syn(pool::SpeciesPool, grid::Grid; τ::Float64,
                             keepmask::BitVector, eta_frac::Float64=0.02)
    B00 = eval_state(pool, grid; τ=τ, keepmask=nothing,   cascade=false, eta_frac=eta_frac)
    B10 = eval_state(pool, grid; τ=τ, keepmask=keepmask,  cascade=false, eta_frac=eta_frac)
    B01 = eval_state(pool, grid; τ=τ, keepmask=nothing,   cascade=true,  eta_frac=eta_frac)
    B11 = eval_state(pool, grid; τ=τ, keepmask=keepmask,  cascade=true,  eta_frac=eta_frac)

    clim = B10 - B00
    inter = B01 - B00
    syn = (B11 - B10) - (B01 - B00)        # inclusion–exclusion remainder
    total = B11 - B00
    return (; total, clim, inter, syn)
end

# --- polygon-based stacked bars you already used ---
function stackedbars_poly!(ax::Axis, x::AbstractVector, Y::AbstractMatrix;
                           width=0.8, colors=[:dodgerblue, :orange, :forestgreen])
    n, k = size(Y)
    halfw = width/2
    for i in 1:n
        pos = 0.0; neg = 0.0
        for j in 1:k
            h = Y[i,j]
            x0 = x[i]-halfw; x1 = x[i]+halfw
            if h ≥ 0
                y0, y1 = pos, pos+h
                poly!(ax, Point2f[(x0,y0),(x1,y0),(x1,y1),(x0,y1)]; color=colors[j])
                pos = y1
            else
                y0, y1 = neg+h, neg
                poly!(ax, Point2f[(x0,y0),(x1,y0),(x1,y1),(x0,y1)]; color=colors[j])
                neg = y0
            end
        end
    end
    ax
end

# ---- run for one representative loss (e.g., 50%) ----
begin
    grid = make_grid(60,60; seed=11)
    S = 200; basal_frac = 0.35
    pool = build_pool(S; basal_frac=basal_frac, seed=1,
                      sigma=0.1, density=0.07, pmax=0.70,
                      niche_mode=:bimodal, mu_basal_centers=(0.25,0.75),
                      mu_basal_sd=0.04, b0_basal=0.08, bspread_basal=0.02,
                      b0_cons=0.12,  bspread_cons=0.04)

    keep = 0.5
    τ = 0.55
    hotspot_power = 2.5
    nseeds_cluster = 1
    eta_frac = 0.03   # 3% of kept cells required for regional persistence

    kmR = random_mask(grid.C, keep; seed=101)
    kmC = clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=101)
    kmH = consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep, power=hotspot_power, what=:cons, seed=101)

    dR = decomp_clim_int_syn(pool, grid; τ=τ, keepmask=kmR, eta_frac=eta_frac)
    dC = decomp_clim_int_syn(pool, grid; τ=τ, keepmask=kmC, eta_frac=eta_frac)
    dH = decomp_clim_int_syn(pool, grid; τ=τ, keepmask=kmH, eta_frac=eta_frac)

    # stack Climate, Interaction, Synergy
    Y = vcat(
        [dR.clim  dR.inter  dR.syn],
        [dC.clim  dC.inter  dC.syn],
        [dH.clim  dH.inter  dH.syn],
    )

    fig = Figure(; size=(900, 340))
    ax  = Axis(fig[1,1], ylabel="ΔBSH (50% loss)", xlabel="Scenario")
    x   = 1:3
    stackedbars_poly!(ax, collect(x), Y; width=0.8, colors=[:dodgerblue, :orange, :forestgreen])
    ax.xticks = (x, ["Random","Clustered","Hotspot"])
    ax.xticklabelrotation[] = π/12
    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax,
        [PolyElement(color=:dodgerblue), PolyElement(color=:orange), PolyElement(color=:forestgreen)],
        ["Climate","Interaction","Synergy"]; position=:lt)
    display(fig)
end

