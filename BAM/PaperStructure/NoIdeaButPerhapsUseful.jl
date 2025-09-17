# --- AM branch: never builds a biotic gate (B≡1) ---------------------------
"AM-only assembler (no biotic gate)."
function assemble_AM(pool::SpeciesPool, grid::Grid, A::Matrix{Float64}, keep::BitVector;
                     α::Float64=1.0, μ::Float64=0.0, τA::Float64=0.5, τocc::Float64=0.35,
                     mp::MovementParams=MovementParams())
    S, C = pool.S, grid.C
    M = movement_matrix(pool, grid, A, keep; mp=mp, τA=τA)
    P = falses(S, C)
    for s in 1:S, i in 1:C
        keep[i] || continue
        score = (A[s,i]^α) * (M[s,i]^μ)
        P[s,i] = score ≥ τocc
    end
    return P, ones(Float64, S, C)  # B is identically 1 for AM
end

# --- BAM branch with a HARD fraction gate on diets -------------------------
"""
assemble_BAM_hard(...; θB, kB) uses a fraction-of-prey gate:
  frac = (#prey present)/(#prey total)
  B = frac^kB, and the cell is eligible only if frac ≥ θB.

Set β>0 to make B matter in Score = A^α * B^β * M^μ.
"""
function assemble_BAM_hard(pool::SpeciesPool, grid::Grid, A::Matrix{Float64}, keep::BitVector;
                           α::Float64=1.0, β::Float64=4.0, μ::Float64=0.0,
                           τA::Float64=0.5, τocc::Float64=0.35,
                           θB::Float64=0.75, kB::Int=2,
                           mp::MovementParams=MovementParams())
    S, C = pool.S, grid.C
    M = movement_matrix(pool, grid, A, keep; mp=mp, τA=τA)

    # Build basal first, then consumers, exactly like your assemble_BAM does
    order = sortperm(pool.masses)  # prey -> predators
    P = falses(S, C)
    Bsup = ones(Float64, S, C)

    for s in order
        prey = pool.E[s]
        for i in 1:C
            keep[i] || continue
            if isempty(prey) || pool.basal[s]
                frac = 1.0
            else
                # count present prey at cell i
                present = 0
                for q in prey
                    present += P[q,i] ? 1 : 0
                end
                frac = length(prey) == 0 ? 1.0 : present / length(prey)
            end

            B = frac^kB
            Bsup[s,i] = B

            eligible = (isempty(prey) || pool.basal[s]) ? true : (frac ≥ θB)
            if eligible
                score = (A[s,i]^α) * (B^β) * (M[s,i]^μ)
                P[s,i] = score ≥ τocc
            else
                P[s,i] = false
            end
        end
    end
    return P, Bsup
end

# --- tiny diagnostic: is the biotic gate ever binding? ---------------------
"Share of (s,i) with A≥τA & kept where the biotic gate fails (B below θB)."
function biotic_gate_activation(B::Matrix{Float64}, A::Matrix{Float64}, keep::BitVector;
                                τA::Float64, θB::Float64)
    S, C = size(B)
    num = 0; den = 0
    @inbounds for s in 1:S, i in 1:C
        if keep[i] && A[s,i] ≥ τA
            den += 1
            if B[s,i] < θB; num += 1; end
        end
    end
    return den == 0 ? 0.0 : num/den
end

"Front-on-climate: keep the coolest or warmest tail by climate quantile."
function front_climate_tail_mask(grid::Grid, keep_frac::Float64; tail::Symbol=:warm)
    C = grid.C
    clim = grid.climate
    q = tail === :warm ? quantile(clim, 1 - keep_frac) : quantile(clim, keep_frac)
    keep = BitVector(undef, C)
    if tail === :warm
        @inbounds for i in 1:C; keep[i] = clim[i] ≥ q; end
    else
        @inbounds for i in 1:C; keep[i] = clim[i] ≤ q; end
    end
    keep
end

"""
Plot AM (A×M) vs BAM (A×B×M) for three HL geometries.
Uses assemble_AM() and assemble_BAM_hard() so the maths truly differs.

Key knobs to make BAM differ:
- low link density or diet_cap=2 when building the pool
- θB≈0.75, kB≥2, β≥6
- front = climate-tail (warm by default)
"""
function plot_AM_vs_BAM(; nx=120, ny=120, S=150, basal_frac=0.45,
                         # axes & pool
                         A_level::Symbol=:divergent,
                         B_level::Symbol=:strong,
                         density::Float64=0.06, diet_cap::Int=2,
                         # gates
                         α=1.0, μ=0.0, τA=0.50, τocc=0.50,
                         β=6.0, θB=0.75, kB=2,
                         # loss grid
                         loss_fracs=0.2:0.1:0.8)

    grid = make_grid(nx, ny; seed=42)

    # --- build pool with lean diets (forces non-redundancy) ----------------
    rng = MersenneTwister(123)
    pool = build_pool_from_axes(rng; S, basal_frac, A_level, B_level)
    # enforce sparse metaweb
    for s in 1:pool.S; if !pool.basal[s]; pool.E[s] = pool.E[s][1:min(diet_cap, length(pool.E[s]))]; end; end
    # also random-drop links to reach target density
    if density < 0.3
        for s in 1:pool.S
            if !pool.basal[s] && !isempty(pool.E[s])
                keepn = max(1, round(Int, density * length(pool.E[s]) / 0.3))
                shuffle!(rng, pool.E[s]); pool.E[s] = pool.E[s][1:keepn]
            end
        end
    end

    A = abiotic_matrix(pool, grid)
    mp = MovementParams(; mode=:none, T=8)        # movement OFF for clarity

    # --- helper to compute ΔBSH curve for a given mask maker --------------
    function curve_AM_BAM(make_mask::Function)
        x = collect(loss_fracs)
        dAM = similar(x, Float64); dBAM = similar(x, Float64)
        gate0 = Float64[]; gate1 = Float64[]  # diagnostics at each f

        base = trues(grid.C)
        P0A,_ = assemble_AM(pool, grid, A, base; α, μ, τA, τocc, mp)
        P0B,B0 = assemble_BAM_hard(pool, grid, A, base; α, β, μ, τA, τocc, θB, kB, mp)
        BSH0A = mean(sum(P0A, dims=2)) / grid.C
        BSH0B = mean(sum(P0B, dims=2)) / grid.C
        push!(gate0, biotic_gate_activation(B0, A, base; τA, θB))

        for (k,f) in enumerate(x)
            keep = make_mask(1 - f)
            P1A,_ = assemble_AM(pool, grid, A, keep; α, μ, τA, τocc, mp)
            P1B,B1 = assemble_BAM_hard(pool, grid, A, keep; α, β, μ, τA, τocc, θB, kB, mp)
            BSH1A = mean(sum(P1A, dims=2)) / grid.C
            BSH1B = mean(sum(P1B, dims=2)) / grid.C
            dAM[k]  = BSH1A - BSH0A
            dBAM[k] = BSH1B - BSH0B
            push!(gate1, biotic_gate_activation(B1, A, keep; τA, θB))
        end
        return (; x, dAM, dBAM, gate0=gate0[1], gate1_last=gate1[end])
    end

    # random / clustered / front-on-climate
    rR = curve_AM_BAM(f -> random_mask(MersenneTwister(100), grid.C, f))
    rC = curve_AM_BAM(f -> clustered_mask(MersenneTwister(200), grid, f; nseeds=6))
    rF = curve_AM_BAM(f -> front_climate_tail_mask(grid, f; tail=:warm))

    println("Gate activation baseline→@f=0.6  (R/C/F):")
    println("  Random   : ", round(rR.gate0,digits=2), " → ", round(rR.gate1_last,digits=2))
    println("  Clustered: ", round(rC.gate0,digits=2), " → ", round(rC.gate1_last,digits=2))
    println("  Front    : ", round(rF.gate0,digits=2), " → ", round(rF.gate1_last,digits=2))

    # -------------------------- figure ------------------------------------
    fig = Figure(; size=(1400, 520))
    Label(fig[0,1:2], "HL effect with vs without biotic — A=$(A_level), M=off"; fontsize=18)

    ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="ΔBSH (mean consumers)",
               title="Climate-only (B OFF)")
    lines!(ax1, rR.x, rR.dAM, label="Random")
    lines!(ax1, rC.x, rC.dAM, label="Clustered")
    lines!(ax1, rF.x, rF.dAM, label="Front (climate-tail)")
    axislegend(ax1, position=:lb, framevisible=false)

    ax2 = Axis(fig[1,2], xlabel="Area lost (fraction)", ylabel="ΔBSH",
               title="Full ABM (B ON) — hard diet gate")
    lines!(ax2, rR.x, rR.dBAM, label="Random")
    lines!(ax2, rC.x, rC.dBAM, label="Clustered")
    lines!(ax2, rF.x, rF.dBAM, label="Front (climate-tail)")
    axislegend(ax2, position=:lb, framevisible=false)

    display(fig)
    # return fig
end

# 1) Make the figure and print gate activation (expect front >> random/clustered)
fig = plot_AM_vs_BAM(; A_level=:divergent, B_level=:strong,
                     density=0.06, diet_cap=2,    # sparse diets
                     β=6.0, θB=0.75, kB=2,       # hard biotic gate
                     τA=0.50, τocc=0.50,         # moderate abiotic gate
                     loss_fracs=0.2:0.1:0.8)

# 2) If BAM ≈ AM, tighten once more (or flip the cut to the tail your consumers use)
fig = plot_AM_vs_BAM(; β=8.0, θB=0.80, kB=3, density=0.05, diet_cap=2)

# --- baseline equalization for a fair AM vs BAM comparison -----------------
# "Find τocc for BAM so the baseline mean BSH matches a target."
function calibrate_tauocc_BAM(target::Float64, pool::SpeciesPool, grid::Grid, A::Matrix{Float64}, keep::BitVector;
                              α::Float64, β::Float64, μ::Float64, τA::Float64, θB::Float64, kB::Int,
                              mp::MovementParams, lo::Float64=0.20, hi::Float64=0.80, iters::Int=26)
    # monotone: higher τocc -> fewer presences -> smaller BSH
    for _ in 1:iters
        mid = (lo+hi)/2
        Pmid, _ = assemble_BAM_hard(pool, grid, A, keep; α, β, μ, τA, τocc=mid, θB, kB, mp)
        val = mean(sum(Pmid, dims=2)) / grid.C
        if val > target
            lo = mid   # occupancy still too high -> tighten threshold
        else
            hi = mid
        end
    end
    return (lo+hi)/2
end

"Share of climate-suitable kept cells where B < θB (biotic gate fails)."
function pfail_B(B::Matrix{Float64}, A::Matrix{Float64}, keep::BitVector; τA::Float64, θB::Float64)
    S, C = size(B)
    num = 0; den = 0
    @inbounds for s in 1:S, i in 1:C
        if keep[i] && A[s,i] ≥ τA
            den += 1
            if B[s,i] < θB; num += 1; end
        end
    end
    return den == 0 ? 0.0 : num/den
end

"""
AM vs BAM with fair baseline, relative losses, DiD and biotic gate diagnostics.

- Baseline equalization: τocc for BAM is tuned so BSH0(BAM) ≈ BSH0(AM).
- Optional species restriction to 'biotic-ready at baseline' (ready≥ϕ_ready).
- Outputs three panels: relative losses, DiD(f) = ΔBAM - ΔAM, and P_fail(f).
"""
function plot_AM_vs_BAM_normalized(; nx=120, ny=120, S=300, basal_frac=0.35,
                                   A_level::Symbol=:divergent,
                                   density::Float64=0.06, diet_cap::Int=2,
                                   α::Float64=1.0, μ::Float64=0.0, τA::Float64=0.50,
                                   τocc_AM::Float64=0.50, β::Float64=6.0, θB::Float64=0.75, kB::Int=2,
                                   ϕ_ready::Float64=0.8, restrict_ready::Bool=true,
                                   loss_fracs=0.2:0.1:0.8)

    # grid & pool (lean diets to expose trophic sensitivity)
    grid = make_grid(nx, ny; seed=42)
    rng  = MersenneTwister(123)
    pool = build_pool_from_axes(rng; S, basal_frac, A_level, B_level=:strong)
    for s in 1:pool.S
        if !pool.basal[s]
            pool.E[s] = pool.E[s][1:min(diet_cap, length(pool.E[s]))]
            if density < 0.3 && length(pool.E[s]) > 1
                keepn = max(1, round(Int, density * length(pool.E[s]) / 0.3))
                shuffle!(rng, pool.E[s]); pool.E[s] = pool.E[s][1:keepn]
            end
        end
    end

    A  = abiotic_matrix(pool, grid)
    mp = MovementParams(; mode=:none, T=8)   # keep movement off for a pure trophic test
    cons = consumer_mask(pool)

    base = trues(grid.C)

    # --- Baseline AM
    P0A,_ = assemble_AM(pool, grid, A, base; α, μ, τA, τocc=τocc_AM, mp)
    BSH0A = vec(sum(P0A, dims=2)) ./ grid.C

    # --- Baseline BAM with τocc calibrated to match AM baseline mean
    target_mean = mean(BSH0A[cons])
    τocc_BAM = calibrate_tauocc_BAM(target_mean, pool, grid, A, base;
                                    α, β, μ, τA, θB, kB, mp)
    P0B, B0 = assemble_BAM_hard(pool, grid, A, base; α, β, μ, τA, τocc=τocc_BAM, θB, kB, mp)
    BSH0B = vec(sum(P0B, dims=2)) ./ grid.C

    # --- choose species set (restrict to 'biotic-ready' if asked)
    ready = trues(pool.S)
    if restrict_ready
        # fraction of climate-suitable baseline cells that pass the B gate
        fracA = zeros(Float64, pool.S); fracB = zeros(Float64, pool.S)
        @inbounds for s in 1:pool.S, i in 1:grid.C
            if A[s,i] ≥ τA
                fracA[s] += 1
                if B0[s,i] ≥ θB; fracB[s] += 1; end
            end
        end
        for s in 1:pool.S
            if fracA[s] > 0
                ready[s] = (fracB[s]/fracA[s]) ≥ ϕ_ready
            else
                ready[s] = false
            end
        end
    end
    keep_spp = cons .& ready

    # --- per-geometry curves
    function curves(make_mask)
        x    = collect(loss_fracs)
        relA = similar(x); relB = similar(x)
        did  = similar(x); pf  = similar(x)
        BSH0A_mean = mean(BSH0A[keep_spp])
        BSH0B_mean = mean(BSH0B[keep_spp])

        for (k,f) in enumerate(x)
            keep = make_mask(1 - f)

            P1A,_ = assemble_AM(pool, grid, A, keep; α, μ, τA, τocc=τocc_AM, mp)
            P1B,B1 = assemble_BAM_hard(pool, grid, A, keep; α, β, μ, τA, τocc=τocc_BAM, θB, kB, mp)

            BSH1A = vec(sum(P1A, dims=2)) ./ grid.C
            BSH1B = vec(sum(P1B, dims=2)) ./ grid.C

            ΔA = mean(BSH1A[keep_spp]) - BSH0A_mean
            ΔB = mean(BSH1B[keep_spp]) - BSH0B_mean
            relA[k] = ΔA / max(BSH0A_mean, eps())
            relB[k] = ΔB / max(BSH0B_mean, eps())
            did[k]  = ΔB - ΔA
            pf[k]   = pfail_B(B1, A, keep; τA, θB)
        end
        return (; x, relA, relB, did, pf, τocc_BAM)
    end

    makeR = f -> random_mask(MersenneTwister(100), grid.C, f)
    makeC = f -> clustered_mask(MersenneTwister(200), grid, f; nseeds=6)
    makeF = f -> front_climate_tail_mask(grid, f; tail=:warm)

    rR = curves(makeR); rC = curves(makeC); rF = curves(makeF)

    # ----------------- figure -----------------
    fig = Figure(; size=(1500, 520))
    Label(fig[0,1:3], "AM vs BAM (baseline equalized) — A=$(A_level), M=off  |  τocc_BAM=$(round(rR.τocc_BAM,digits=2))"; fontsize=18)

    ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Relative loss (Δ/BSH₀)",
               title="Relative ΔBSH: AM (dashed) vs BAM (solid)")
    for (col, rr) in enumerate((rR,rC,rF))
        lines!(ax1, rr.x, rr.relB; label=col==1 ? "BAM - Random" : nothing)
        lines!(ax1, rr.x, rr.relA; linestyle=:dash, label=col==1 ? "AM - Random" : nothing)
    end
    lines!(ax1, rC.x, rC.relB; color=:orange, label="BAM - Clustered")
    lines!(ax1, rC.x, rC.relA; color=:orange, linestyle=:dash, label="AM - Clustered")
    lines!(ax1, rF.x, rF.relB; color=:seagreen, label="BAM - Front")
    lines!(ax1, rF.x, rF.relA; color=:seagreen, linestyle=:dash, label="AM - Front")
    axislegend(ax1, position=:lb, framevisible=false)

    ax2 = Axis(fig[1,2], xlabel="Area lost (fraction)", ylabel="DiD(f) = ΔBAM − ΔAM",
               title="Biotic amplification (difference-in-differences)")
    lines!(ax2, rR.x, rR.did; label="Random")
    lines!(ax2, rC.x, rC.did; label="Clustered")
    lines!(ax2, rF.x, rF.did; label="Front")
    hlines!(ax2, [0.0]; linestyle=:dot, color=:gray)
    axislegend(ax2, position=:lb, framevisible=false)

    ax3 = Axis(fig[1,3], xlabel="Area lost (fraction)", ylabel="P_fail(f)",
               title="Biotic gate failure among A-suitable kept cells")
    lines!(ax3, rR.x, rR.pf; label="Random")
    lines!(ax3, rC.x, rC.pf; label="Clustered")
    lines!(ax3, rF.x, rF.pf; label="Front")
    axislegend(ax3, position=:lt, framevisible=false)

    display(fig)
    # return fig
end

fig = plot_AM_vs_BAM_normalized(; A_level=:divergent,
                                density=0.04, diet_cap=2,  # lean metaweb
                                β=6.0, θB=0.75, kB=2,     # meaningful biotic gate
                                τA=0.50, τocc_AM=0.50,    # same AM gate as before
                                restrict_ready=true, ϕ_ready=0.8)
