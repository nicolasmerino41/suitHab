# -----------------------------
# (A) Nonlinear viability metric
# -----------------------------
# We still compute A_s (climate area fraction) and p_s (prey support among suitable),
# but the outcome is now: U(A,p) = A * φ(p), with φ nonlinear.
# Choose either a saturating φ (default) or a hard threshold.

# Saturating φ(p) = p / (p + h)   (h is "half-saturation"; larger h => stronger nonlinearity)
φ_saturating(p; h=0.20) = p <= 0 ? 0.0 : p / (p + h)

# Threshold φ(p) = 1[p ≥ ρ]        (ρ is the required prey support level)
φ_threshold(p; ρ=0.5) = (isnan(p) || p < ρ) ? 0.0 : 1.0

"""
U(A, p; mode=:saturating, h=0.20, ρ=0.5)

Return the viability outcome for a species given:
- A: fraction of climatically suitable cells (0-1)
- p: fraction of those suitable cells that are prey-supported (0-1), with NaN allowed when A=0
Nonlinear φ creates genuine curvature → synergy becomes non-zero.
"""
function U(A::Float64, p::Float64; mode::Symbol=:saturating, h::Float64=0.20, ρ::Float64=0.5)
    if !isfinite(A) || A <= 0.0
        return 0.0
    end
    if mode === :saturating
        return A * φ_saturating(isnan(p) ? 0.0 : p; h=h)
    elseif mode === :threshold
        return A * φ_threshold(isnan(p) ? 0.0 : p; ρ=ρ)
    else
        error("Unknown mode = $mode. Use :saturating or :threshold.")
    end
end

# ---------------------------------------------------------
# (B) Symmetric (midpoint) decomposition for GENERAL U(A,p)
# ---------------------------------------------------------
# For any outcome U(A,p), define:
#   ΔU        = U(A1,p1) - U(A0,p0)
#   Climate   = U(A1, p̄) - U(A0, p̄)      (vary A at midpoint p̄)
#   Interaction = U(Ā, p1) - U(Ā, p0)     (vary p at midpoint Ā)
#   Synergy   = ΔU - Climate - Interaction
# This equals 0 only when U is bilinear (our old B = A*p). With nonlinear φ, synergy ≠ 0.

struct Decomp{T}
    dB::T; dA_only::T; dInt_only::T; synergy::T
end

function decompose_viability(A0::Float64, p0::Float64, A1::Float64, p1::Float64;
                             mode::Symbol=:saturating, h::Float64=0.20, ρ::Float64=0.5)
    Ā   = 0.5*(A0 + A1)
    p̄   = 0.5*(p0 + p1)
    U0  = U(A0, p0; mode=mode, h=h, ρ=ρ)
    U1  = U(A1, p1; mode=mode, h=h, ρ=ρ)
    ΔU  = U1 - U0
    dCl = U(A1, p̄; mode=mode, h=h, ρ=ρ) - U(A0, p̄; mode=mode, h=h, ρ=ρ)
    dIn = U(Ā,  p1; mode=mode, h=h, ρ=ρ) - U(Ā,  p0; mode=mode, h=h, ρ=ρ)
    dSy = ΔU - dCl - dIn
    return Decomp(ΔU, dCl, dIn, dSy)
end

# ---------------------------------------------------
# (C) Cell-removal masks (random / clustered / hotspot)
# ---------------------------------------------------

# (You already have random_mask and clustered_mask)
# Hotspot: rank cells by "prey-supported climate potential" using climate pass only.
function hotspot_mask(Z0::BitMatrix, pool::SpeciesPool, keep_frac::Float64)
    S, C = size(Z0)
    cons = .!pool.basal
    score = zeros(Float64, C)
    @inbounds for c in 1:C
        sc = 0.0
        for s in findall(cons)
            if Z0[s,c]
                # prey potentially present in same cell under climate pass?
                has = false
                for q in pool.E[s]
                    if Z0[q,c]; has = true; break; end
                end
                sc += has ? 1.0 : 0.0
            end
        end
        score[c] = sc
    end
    nkeep = max(1, round(Int, keep_frac * C))
    idx   = sortperm(score; rev=true)       # high score kept first
    keep  = falses(C); keep[idx[1:nkeep]] .= true
    return keep
end

# ---------------------------------------------------
# (D) One-step decomposition at a given keep mask
# ---------------------------------------------------
"""
decomp_at_viability(pool, grid; τ, keepmask, mode=:saturating, h=0.20, ρ=0.5)

Returns mean (over consumers) of ΔU, climate-only, interaction-only, and synergy.
"""
function decomp_at_viability(pool::SpeciesPool, grid::Grid;
                             τ::Float64=0.5, keepmask::BitVector,
                             mode::Symbol=:saturating, h::Float64=0.20, ρ::Float64=0.5)
    Z0 = climate_pass(pool, grid; τ=τ)
    P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask)
    P1 = assemble(Z1, pool)

    cons = .!pool.basal

    # A_s and p_s under "before" and "after"
    A0 = climate_area_per_species(Z0)[cons]
    p0 = bsh1_cond_per_species(P0, Z0, pool)[cons]
    A1 = climate_area_per_species(Z1)[cons]
    p1 = bsh1_cond_per_species(P1, Z1, pool)[cons]

    # species-wise decomposition, then average
    d = map((a0,p0v,a1,p1v)->decompose_viability(a0,p0v,a1,p1v; mode=mode, h=h, ρ=ρ), A0,p0,A1,p1)
    dB   = mean(getfield.(d, :dB))
    dCl  = mean(getfield.(d, :dA_only))
    dIn  = mean(getfield.(d, :dInt_only))
    dSy  = mean(getfield.(d, :synergy))
    return (dB=dB, dA_only=dCl, dInt_only=dIn, synergy=dSy)
end

# ---------------------------------------------------
# (E) Curves vs loss for different spatial strategies
# ---------------------------------------------------
"""
decomp_curve_viability(kind; grid, pool, τ, losses, mode, h, ρ, nseeds_cluster)

kind ∈ (:random, :clustered, :hotspot)
Returns NamedTuple of vectors over `losses`.
"""
function decomp_curve_viability(kind::Symbol; grid::Grid, pool::SpeciesPool,
                                τ::Float64=0.5, losses=0.0:0.05:0.8,
                                mode::Symbol=:saturating, h::Float64=0.20, ρ::Float64=0.5,
                                nseeds_cluster::Int=6)
    Z0 = climate_pass(pool, grid; τ=τ)
    out = (; loss = Float64[], dA_only=Float64[], dInt_only=Float64[], synergy=Float64[], dTotal=Float64[])
    for f in losses
        keep = 1.0 - f
        km = if kind === :random
            random_mask(grid.C, keep; seed=round(Int, 10_000*f)+11)
        elseif kind === :clustered
            clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=round(Int, 10_000*f)+37)
        elseif kind === :hotspot
            hotspot_mask(Z0, pool, keep)
        else
            error("Unknown kind = $kind")
        end
        r = decomp_at_viability(pool, grid; τ=τ, keepmask=km, mode=mode, h=h, ρ=ρ)
        push!(out.loss, f)
        push!(out.dA_only,     r.dA_only)
        push!(out.dInt_only, r.dInt_only)
        push!(out.synergy,     r.synergy)
        push!(out.dTotal,       r.dB)
    end
    return out
end

# ---------------------------------------------------
# (F) Minimal demo / plotting
# ---------------------------------------------------
# Example params (you can swap in your RS regimes). Choose a setup that
# creates low redBndancy + high synchrony so that hotspot targeting bites.

grid = make_grid(60, 60; seed=11)

pool = build_pool(200;
    basal_frac=0.35, seed=1,
    # diets: few prey per consumer
    sigma=0.22, density=0.12, pmax=0.70, R0_mean=10.0, R0_sd=0.25,
    # prey synchrony via basal niches (bimodal)
    niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.04,
    b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04
)

losses = 0.0:0.05:0.8
τ = 0.55

R = decomp_curve_viability(:random;   grid=grid, pool=pool, τ=τ, losses=losses, mode=:threshold, h=0.1, ρ=0.6)
C = decomp_curve_viability(:clustered;grid=grid, pool=pool, τ=τ, losses=losses, mode=:threshold, h=0.1, ρ=0.6, nseeds_cluster=1)
H = decomp_curve_viability(:hotspot;  grid=grid, pool=pool, τ=τ, losses=losses, mode=:threshold, h=0.1, ρ=0.6)

begin
    fig = Figure(; size=(1440,420))
    lab = ["Climate-only","Interaction-only","Synergy"]
    for (j,(name,dat)) in enumerate((("Random",R),("Clustered",C),("Hotspot",H)))
        ax = Axis(fig[1,j], title=name, xlabel="Area lost (fraction)", ylabel=j==1 ? "ΔBSH decomposition" : "")
        lines!(ax, dat.loss, dat.dA_only,     label=lab[1])
        lines!(ax, dat.loss, dat.dInt_only, label=lab[2])
        lines!(ax, dat.loss, dat.synergy,     label=lab[3])
        hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
        if j==3; axislegend(ax, position=:lt); end
    end
    display(fig)
end

# --- Excess damage (strategy minus random) --------------------
function totals_curve(kind::Symbol; grid::Grid, pool::SpeciesPool, τ=0.55,
                      losses=0.0:0.05:0.8, mode::Symbol=:saturating, h=0.25,
                      nseeds_cluster::Int=1)
    Z0 = climate_pass(pool, grid; τ=τ)
    tot = Float64[]
    for f in losses
        keep = 1.0 - f
        km = if kind === :random
            random_mask(grid.C, keep; seed=round(Int, 10_000*f)+11)
        elseif kind === :clustered
            clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=round(Int, 10_000*f)+37)
        elseif kind === :hotspot
            hotspot_mask(Z0, pool, keep)
        else
            error("kind must be :random, :clustered or :hotspot")
        end
        r = decomp_at_viability(pool, grid; τ=τ, keepmask=km, mode=mode, h=h)
        push!(tot, r.dB)  # mean ΔU over consumers
    end
    return collect(losses), tot
end

# Example run (same pool & grid you used above)
losses = 0.0:0.05:0.8
Lr, Ur = totals_curve(:random;   grid=grid, pool=pool, τ=τ, losses=losses, mode=:saturating, h=0.25)
Lc, Uc = totals_curve(:clustered;grid=grid, pool=pool, τ=τ, losses=losses, mode=:saturating, h=0.25, nseeds_cluster=1)
Lh, Uh = totals_curve(:hotspot;  grid=grid, pool=pool, τ=τ, losses=losses, mode=:saturating, h=0.25)

ExC = Uc .- Ur   # clustered - random
ExH = Uh .- Ur   # hotspot   - random

begin
    fig = Figure(; size=(780,360))
    ax  = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Excess ΔU vs random")
    lines!(ax, losses, ExC, label="Clustered - Random")
    lines!(ax, losses, ExH, label="Hotspot - Random")
    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax, position=:lb)
    display(fig)
end
