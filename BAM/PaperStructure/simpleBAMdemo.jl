# colors for geometries
const COL = Dict(:random => RGBf(0.24,0.49,0.86),
                 :clustered => RGBf(0.90,0.55,0.18),
                 :front => RGBf(0.18,0.70,0.48))

# ---------- helpers ----------

# bisection to match BAM baseline (mean BSH at full keep) to AM baseline
function _equalize_tauocc!(pool, grid, A; bam::BAMParams, mp::MovementParams, target_mean::Float64)
    keep_all = trues(grid.C)
    lo, hi = 0.05, 0.90
    for _ in 1:25
        mid = (lo+hi)/2
        bam2 = BAMParams(; α=bam.α, β=bam.β, μ=bam.μ, γ=bam.γ, τA=bam.τA, τocc=mid)
        P, Bsup, Score, M = assemble_BAM(pool, grid, A, keep_all; bam=bam2, mp=mp)
        S = species_stats(pool, grid, A, keep_all, P, Bsup, M, bam2)
        m = mean(S.BSH[.!pool.basal])
        if m > target_mean
            lo = mid
        else
            hi = mid
        end
    end
    return (lo+hi)/2
end

# build masks
function _mask_for(grid::Grid, rng::AbstractRNG, hl_kind::Symbol, keep_frac; nseeds_cluster::Int, front_axis::Symbol, front_noise::Float64)
    return hl_kind === :random    ? random_mask(rng, grid.C, keep_frac) :
           hl_kind === :clustered ? clustered_mask(rng, grid, keep_frac; nseeds=nseeds_cluster) :
           hl_kind === :front     ? frontlike_mask(rng, grid, keep_frac; axis=front_axis, noise=front_noise) :
           error("Unknown HL $(hl_kind)")
end

# per-species relative loss (BSH1 vs BSH0) for AM or BAM
function _rel_losses(pool, grid, A, keep; bam::BAMParams, mp::MovementParams)
    keep0 = trues(grid.C)
    P0, B0, _, M0 = assemble_BAM(pool, grid, A, keep0; bam=bam, mp=mp)
    S0 = species_stats(pool, grid, A, keep0, P0, B0, M0, bam)
    P1, B1, _, M1 = assemble_BAM(pool, grid, A, keep;  bam=bam, mp=mp)
    S1 = species_stats(pool, grid, A, keep,  P1, B1, M1, bam)
    cons = .!pool.basal
    rel = similar(S1.BSH)
    rel .= 0.0
    @inbounds for s in 1:pool.S
        if cons[s]
            b0 = S0.BSH[s]
            b1 = S1.BSH[s]
            rel[s] = (b1 - b0) / (b0 + 1e-12)
        end
    end
    return rel[cons]
end

# P_fail among A-suitable kept cells (mean across consumers)
function _pfail(pool, grid, A, keep; bam::BAMParams, mp::MovementParams, b_fail::Float64=0.5)
    P, Bsup, _, M = assemble_BAM(pool, grid, A, keep; bam=bam, mp=mp)
    cons = .!pool.basal
    out = Float64[]
    @inbounds for s in 1:pool.S
        cons[s] || continue
        idx = [i for i in 1:grid.C if keep[i] && (A[s,i] ≥ bam.τA)]
        if isempty(idx)
            push!(out, 0.0)
        else
            push!(out, mean(b -> b < b_fail ? b : b_fail, view(Bsup, s, idx)))
        end
    end
    return mean(out)
end

# ---------- main computation ----------

"""
run_simple_bam_demo(; nx,ny,S, A_level, B_level, specialist_frac, diet_cap, β, γ,
                     fstar, fracs, M_level, rng_seeds...)

Builds a low-redundancy web, equalizes BAM baseline to AM, and returns:
- rel_AM[g], rel_BAM[g]  (relative per-species losses at f*=keep loss)
- pfail[g][k]            (curve of P_fail over fracs, BAM only)
- worst_AM, worst_BAM    (argmin ΔBSH at f*)
- plus the figure.
"""
function run_simple_bam_demo(; nx::Int=60, ny::Int=60, S::Int=150, basal_frac::Float64=0.45,
    A_level::Symbol=:divergent, B_level::Symbol=:strong, M_level::Symbol=:off,
    specialist_frac::Float64=0.30, diet_cap::Int=2, β::Float64=7.0, γ::Float64=8.0,
    τA::Float64=0.5, τocc_AM::Float64=0.35, τocc_BAM_init::Float64=0.35,
    fracs = collect(0.2:0.1:0.8), fstar::Float64=0.60,
    seeds_mask = 1:4, nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
    sim_seed::Int=2025)

    grid = make_grid(nx, ny; seed=42)

    # --- build pool (low redundancy with a slice of specialists)
    rng_pool = MersenneTwister(sim_seed)
    pool = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level, specialist_frac)
    # also cap general diets gently
    _trim_diets!(rng_pool, pool, diet_cap)
    A = abiotic_matrix(pool, grid)

    # --- set AM vs BAM parameters
    pars_AM  = bam_from_axes(; B_level=:none, M_level, τA=τA, τocc=τocc_AM, β_override=0.0, γ_override=2.0)
    bamAM, mpAM = pars_AM.bam, pars_AM.mp

    pars_BAM = bam_from_axes(; B_level, M_level, τA=τA, τocc=τocc_BAM_init, β_override=β, γ_override=γ)
    bamBAM, mpBAM = pars_BAM.bam, pars_BAM.mp

    # --- equalize BAM baseline to AM (mean BSH at keep=all)
    keep_all = trues(grid.C)
    P0A,B0A,_,M0A = assemble_BAM(pool, grid, A, keep_all; bam=bamAM, mp=mpAM)
    SA0 = species_stats(pool, grid, A, keep_all, P0A, B0A, M0A, bamAM)
    target = mean(SA0.BSH[.!pool.basal])
    τeq = _equalize_tauocc!(pool, grid, A; bam=bamBAM, mp=mpBAM, target_mean=target)
    bamBAM = BAMParams(; α=bamBAM.α, β=bamBAM.β, μ=bamBAM.μ, γ=bamBAM.γ, τA=bamBAM.τA, τocc=τeq)

    # --- curves and per-species at f*
    geoms = (:random, :clustered, :front)
    rel_AM = Dict{Symbol,Vector{Float64}}()
    rel_BAM = Dict{Symbol,Vector{Float64}}()
    pfail = Dict{Symbol,Vector{Float64}}()
    dBSH_AM = Dict{Symbol,Vector{Float64}}()
    dBSH_BAM = Dict{Symbol,Vector{Float64}}()

    for g in geoms
        rng_m = MersenneTwister(hash((sim_seed,:mask,g)))
        # curves
        dA = Float64[]; dB = Float64[]; pF = Float64[]
        for f in fracs
            keep = _mask_for(grid, rng_m, g, 1-f; nseeds_cluster, front_axis, front_noise)
            # AM
            P0,B0,_,M0 = assemble_BAM(pool, grid, A, keep_all; bam=bamAM, mp=mpAM)
            S0 = species_stats(pool, grid, A, keep_all, P0, B0, M0, bamAM)
            P1,B1,_,M1 = assemble_BAM(pool, grid, A, keep;    bam=bamAM, mp=mpAM)
            S1 = species_stats(pool, grid, A, keep,    P1, B1, M1, bamAM)
            push!(dA, mean(S1.BSH[.!pool.basal]) - mean(S0.BSH[.!pool.basal]))
            # BAM
            P0b,B0b,_,M0b = assemble_BAM(pool, grid, A, keep_all; bam=bamBAM, mp=mpBAM)
            S0b = species_stats(pool, grid, A, keep_all, P0b, B0b, M0b, bamBAM)
            P1b,B1b,_,M1b = assemble_BAM(pool, grid, A, keep;    bam=bamBAM, mp=mpBAM)
            S1b = species_stats(pool, grid, A, keep,    P1b, B1b, M1b, bamBAM)
            push!(dB, mean(S1b.BSH[.!pool.basal]) - mean(S0b.BSH[.!pool.basal]))
            push!(pF, _pfail(pool, grid, A, keep; bam=bamBAM, mp=mpBAM, b_fail=0.5))
        end
        dBSH_AM[g] = dA
        dBSH_BAM[g] = dB
        pfail[g] = pF

        # per-species relative losses at f*
        k = findall(x->isapprox(x,fstar; atol=1e-6), fracs)
        if isempty(k); error("f*=$(fstar) not in fracs"); end
        keepf = _mask_for(grid, rng_m, g, 1-fracs[k[1]]; nseeds_cluster, front_axis, front_noise)
        rel_AM[g]  = _rel_losses(pool, grid, A, keepf; bam=bamAM,  mp=mpAM)
        rel_BAM[g] = _rel_losses(pool, grid, A, keepf; bam=bamBAM, mp=mpBAM)
    end

    # worst geometry at f*
    kstar = findall(x->isapprox(x,fstar; atol=1e-6), fracs)[1]
    worst_AM  = argmin([dBSH_AM[g][kstar] for g in geoms]) |> i->geoms[i]
    worst_BAM = argmin([dBSH_BAM[g][kstar] for g in geoms]) |> i->geoms[i]

    # ---------- plot ----------
    fig = Figure(; size=(1400, 520))
    Label(fig[0,:], "Simple test — A=$(A_level), B=$(B_level), M=$(M_level)  |  baseline equalized  |  f*=$(fstar)", fontsize=18)

    # (1) CDF panel
    ax1 = Axis(fig[1,1], xlabel="relative loss (−ΔBSH/BSH₀)", ylabel="Fraction of consumers", title="Per-species loss CDF at f*")
    for g in geoms
        # plot CDF of positive losses
        L_AM  = -rel_AM[g]
        L_BAM = -rel_BAM[g]
        for (L,ls) in ((L_AM,:dash), (L_BAM,:solid))
            v = sort(L)
            y = range(0,1; length=length(v))
            lines!(ax1, v, y; color=COL[g], linestyle=ls, label=(ls==:solid ? String(g)*"  BAM" : String(g)*"  AM"))
        end
    end
    axislegend(ax1; position=:rb, framevisible=false)

    # (2) Rank flip mini-panel
    ax2 = Axis(fig[1,2], xlabel="", ylabel="", title="Worst geometry at f*", xticks=(1:2, ["AM","BAM"]))
    hidedecorations!(ax2, grid=false)
    xs = [1.0, 2.0]; ys = [0.5, 0.5]
    # draw the three geometry symbols at each x; highlight the worst
    function _mark!(x, worst)
        scatter!(ax2, [x],[0.6]; marker=:circle, markersize=10, color=(worst==:random ? COL[:random] : RGBAf(0,0,0,0.1)))
        scatter!(ax2, [x],[0.5]; marker=:rect,   markersize=10, color=(worst==:clustered ? COL[:clustered] : RGBAf(0,0,0,0.1)))
        scatter!(ax2, [x],[0.4]; marker=:utriangle, markersize=11, color=(worst==:front ? COL[:front] : RGBAf(0,0,0,0.1)))
    end
    _mark!(1.0, worst_AM); _mark!(2.0, worst_BAM)
    text!(ax2, 1.0, 0.85, text=string(worst_AM), align=(:center,:bottom), color=:gray25)
    text!(ax2, 2.0, 0.85, text=string(worst_BAM), align=(:center,:bottom), color=:gray25)

    # (3) P_fail curves (BAM only)
    ax3 = Axis(fig[1,3], xlabel="Area lost (fraction)", ylabel="P_fail (A-suitable kept cells)", title="Biotic insufficiency (BAM)")
    for g in geoms
        lines!(ax3, fracs, pfail[g]; color=COL[g], label=String(g))
    end
    axislegend(ax3; position=:rb, framevisible=false)

    display(fig)
    return (; fig, rel_AM, rel_BAM, pfail, dBSH_AM, dBSH_BAM, worst_AM, worst_BAM)
end

res = run_simple_bam_demo(
    nx=60, ny=60, S=150, basal_frac=0.45,
    A_level=:divergent,       # climate gradient so "front" is meaningful
    B_level=:strong,          # we use overrides below anyway
    M_level=:off,             # keep movement out of the story
    specialist_frac=0.30,     # ~30% true specialists
    diet_cap=2,               # low redundancy elsewhere
    β=7.0, γ=8.0,             # biotic gate bites but is smooth
    τA=0.5, τocc_AM=0.35,     # same as your usual
    fracs=collect(0.2:0.1:0.8),
    fstar=0.60,               # comparison loss fraction
    seeds_mask=1:4, nseeds_cluster=6, front_axis=:x, front_noise=0.04,
    sim_seed=2025
)