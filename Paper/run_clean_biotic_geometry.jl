# === run_clean_biotic_geometry.jl ============================================
# Purpose: Minimal, clean tests for geometry-specific biotic effects under HL
# Uses ONLY your modules (no refactors), adds a few helper functions locally.
#
# Outputs go to Paper/clean/
# ==============================================================================
# plotting
include("../SetUp.jl")
include("src/grids.jl");      using .Grids
include("src/metawebs.jl");   using .Metawebs
include("src/hl.jl");         using .HL
include("src/bsh.jl");        using .BSH
include("src/metrics.jl");    using .Metrics

# ---------------- SETTINGS (feel free to tweak) ----------------
rng         = MersenneTwister(123)
nx, ny      = 120, 120
S           = 200
basal_frac  = 0.1
outdir      = "Paper/clean"; isdir(outdir) || mkpath(outdir)

# Grids
grid_grad   = Grids.make_grid_gradient(nx, ny; seed=11)
grid_patch  = Grids.make_grid_patchy(nx, ny;   seed=12)
grid_mosa = Grids.make_grid_mosaic(nx, ny;   seed=13)
grid_ridge  = Grids.make_grid_ridge(nx, ny;    seed=14)
grid        = grid_grad                         # main grid for clean tests

# Metawebs (redundancy ladder)
pool_low  = Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:low)
pool_mid  = Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:mid)
pool_high = Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:high)

# BAM parameters — movement OFF for the cleanest biotic readout
pars = BSH.BAMParams(; τA=0.5, τB=0.5, τocc=0.2, γ=3.0, movement=:off, T=8)

# HL configs
loss_fracs = 0.2:0.05:0.8
fstar      = 0.60
geoms      = (:random, :clustered, :front)

# Abiotic matrix builder (use your default; swap to `abiotic_matrix_aligned` if needed)
A_fn = BSH.abiotic_matrix

# ============================================================================ #
#                            HELPER FUNCTIONS                                  #
# ============================================================================ #

"""
Prey sufficiency rate, conditioned on A-suitable & kept, at one f and geometry.
Uses AM presence of prey (no B gate) and requires at least `kreq` prey.
Returns a single number in [0,1].
"""
function prey_sufficiency_rate(rng, pool, grid, pars; f::Float64, geom::Symbol,
                               seed_A::Int=1, kreq::Int=1, A_fn=BSH.abiotic_matrix)
    A = A_fn(pool, grid; seed=seed_A)

    keepfrac = 1 - f
    keep = geom === :random    ? HL.random_mask(rng, grid.C, keepfrac) :
           geom === :clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                 HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)

    # A-suitable & kept cells (movement is OFF here)
    Aok = (A .>= pars.τA) .& reshape(keep, 1, :)

    # AM presence (for prey availability)
    P_AM, _ = BSH.assemble_AM(pool, grid, A, keep; pars=pars)

    cons = findall(!, pool.basal)
    prey_sufficient = 0; denom = 0
    for s in cons
        pr = pool.prey[s]
        isempty(pr) && continue
        # count present-prey in eligible cells for consumer s
        Ks = sum(P_AM[pr, :]; dims=1)          # prey present per cell
        elig = vec(Aok[s, :])                  # cells eligible for s
        # Only consider eligible cells
        if any(elig)
            Ss = vec(Ks) .>= kreq
            prey_sufficient += count(Ss[elig])
            denom += count(elig)
        end
    end
    return denom == 0 ? 0.0 : prey_sufficient / denom
end

"""
Geometry-differential BAM effect (DiD vs FRONT) at f*:
BGE_g = (Yg_BAM - Yfront_BAM) - (Yg_AM - Yfront_AM).
Y* are means over consumers (normalized by ORIGINAL area), computed with your helpers.
"""
function biotic_geometry_effect_at_fstar(rng, pool, grid, pars; fstar, A_fn)
    A = A_fn(pool, grid; seed=17)
    function Y(mode::Symbol, geom::Symbol)
        keepfrac = 1 - fstar
        keep = geom === :random    ? HL.random_mask(rng, grid.C, keepfrac) :
               geom === :clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                     HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)
        return BSH.mean_BSH_over_consumers_area0(rng, pool, grid, pars; A=A, keepmask=keep, mode=mode)
    end
    YAM = Dict(g => Y(:AM, g)  for g in geoms)
    YBM = Dict(g => Y(:BAM, g) for g in geoms)
    BGE = Dict{Symbol,Float64}()
    for g in geoms
        if g == :random; continue; end
        BGE[g] = (YBM[g] - YBM[:random]) - (YAM[g] - YAM[:random])
    end
    return (; YAM, YBM, BGE)
end

"""
Per-species penalty vs diet-size bins at f*:
Returns Dict geometry => (kbin, mean_penalty, lo, hi)
Penalty is (b_BAM - b_AM) per species (fraction of ORIGINAL area), consumers only.
Bins: k=1, k=2, k>=3.
"""
function penalty_vs_diet_bins(rng, pool, grid, pars; fstar, geom::Symbol, A_fn)
    A = A_fn(pool, grid; seed=33)
    keepfrac = 1 - fstar
    keep = geom === :random    ? HL.random_mask(rng, grid.C, keepfrac) :
           geom === :clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                 HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)

    PAM, bAM   = BSH.assemble_AM(pool,  grid, A, keep; pars=pars)
    PBM, bBM   = begin
        tmp = BSH.assemble_BAM(pool, grid, A, keep; pars=pars, agg=:kofn, kreq=2)
        (tmp.P, tmp.bsh)
    end

    cons = findall(!, pool.basal)
    ks   = [length(pool.prey[s]) for s in cons]
    pen  = bBM[cons] .- bAM[cons]

    # bins
    function binlab(k)
        k == 1 && return "k=1"
        k == 2 && return "k=2"
        return "k≥3"
    end
    bins = Dict("k=1"=>Float64[], "k=2"=>Float64[], "k≥3"=>Float64[])
    for (k, p) in zip(ks, pen)
        push!(bins[binlab(min(k, 3))], p)
    end
    labels = ["k=1","k=2","k≥3"]
    μ  = [isempty(bins[l]) ? NaN : mean(bins[l])         for l in labels]
    lo = [isempty(bins[l]) ? NaN : quantile(bins[l],0.10) for l in labels]
    hi = [isempty(bins[l]) ? NaN : quantile(bins[l],0.90) for l in labels]
    return (kbin=labels, mean=μ, lo=lo, hi=hi)
end

"""
Co-retention ECDF: for each consumer s and eligible cell c (kept & A≥τA),
count K_sc = # of that consumer's prey present under AM.
Return the ECDF arrays for each geometry at the same f.
"""
function coretention_ecdf(rng, pool, grid, pars; f::Float64, geoms, A_fn)
    A = A_fn(pool, grid; seed=44)
    function Kvals(geom)
        keepfrac = 1 - f
        keep = geom === :random    ? HL.random_mask(rng, grid.C, keepfrac) :
               geom === :clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                     HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)
        P_AM, _ = BSH.assemble_AM(pool, grid, A, keep; pars=pars)
        cons = findall(!, pool.basal)
        out = Int[]
        for s in cons
            pr = pool.prey[s]
            isempty(pr) && continue
            elig = [keep[i] & (A[s,i] ≥ pars.τA) for i in 1:grid.C]
            if any(elig)
                Ks = sum(P_AM[pr, :]; dims=1)
                append!(out, vec(Ks)[elig])
            end
        end
        return out
    end
    data = Dict(g => Kvals(g) for g in geoms)
    # turn into ECDF grids
    ecdfs = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()
    for g in geoms
        v = sort!(Float64.(data[g]))
        if isempty(v)
            ecdfs[g] = (Float64[], Float64[])
        else
            xs = v
            ys = (0:length(v)-1) ./ length(v)
            ecdfs[g] = (xs, ys)
        end
    end
    return ecdfs
end

"""
Redundancy threshold curve (single grid & geometry):
For each pool in {low, mid, high}, compute Δ = Y_BAM - Y_AM at f*.
Return redundancy proxy (0.95-quantile of diet size) and Δ.
"""
function redundancy_threshold_curve(rng, pools::Vector, grid, pars; fstar, geom::Symbol, A_fn)
    R = Float64[]; Δ = Float64[]
    for pool in pools
        diets = [length(pool.prey[s]) for s in 1:pool.S if !pool.basal[s]]
        push!(R, isempty(diets) ? 0.0 : quantile(diets, 0.95))

        A = A_fn(pool, grid; seed=97)
        keepfrac = 1 - fstar
        keep = geom === :random    ? HL.random_mask(rng, grid.C, keepfrac) :
               geom === :clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                     HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)
        yAM = BSH.mean_BSH_over_consumers_area0(rng, pool, grid, pars; A=A, keepmask=keep, mode=:AM)
        yBM = BSH.mean_BSH_over_consumers_area0(rng, pool, grid, pars; A=A, keepmask=keep, mode=:BAM)
        push!(Δ, yBM - yAM)
    end
    return (R=R, Δ=Δ)
end

# ============================================================================ #
#                                  FIGURES                                     #
# ============================================================================ #

# -- Figure A: Prey-sufficiency vs f (isolated biotic effect) -----------------
function figA_prey_sufficiency_vs_f(rng, pool, grid, pars; loss_fracs, geoms, A_fn, kreq::Int=2)
    fig = Figure(; size=(900,320))
    ax  = Axis(fig[1,1], xlabel="area lost f", ylabel="P(sufficient prey | kept ∩ A)", title="Prey sufficiency among eligible cells")
    for (col, g) in enumerate(geoms)
        ys = [prey_sufficiency_rate(rng, pool, grid, pars; f=f, geom=g, kreq=kreq, A_fn=A_fn) for f in loss_fracs]
        lines!(ax, collect(loss_fracs), ys; label=String(g))
    end
    axislegend(ax, position=:rt)
    save(joinpath(outdir, "FigA_prey_sufficiency_vs_f.png"), fig); fig
end

# -- Figure B: Biotic Geometry Effect bars at f* -------------------------------
function figB_BGE_bars(rng, pool, grid, pars; fstar, A_fn)
    res = biotic_geometry_effect_at_fstar(rng, pool, grid, pars; fstar=fstar, A_fn=A_fn)
    labels = ["front","clustered"]
    vals   = [res.BGE[:front], res.BGE[:clustered]]
    fig = Figure(; size=(520,360))
    ax  = Axis(fig[1,1], ylabel="BGE (DiD vs front) at f*=$(fstar)", title="Biotic geometry effect")
    barplot!(ax, 1:2, vals)
    ax.xticks = (1:2, labels)
    hlines!(ax, [0.0], color=:gray, linestyle=:dash)
    save(joinpath(outdir, "FigB_BGE_bars.png"), fig); fig
end

# -- Figure C: Redundancy threshold curve (single geometry) -------------------
function figC_redundancy_threshold(rng, grid, pars; fstar, geom, A_fn)
    pools = [pool_low, pool_mid, pool_high]
    RΔ = redundancy_threshold_curve(rng, pools, grid, pars; fstar=fstar, geom=geom, A_fn=A_fn)
    fig = Figure(; size=(520,360))
    ax  = Axis(fig[1,1], xlabel="redundancy proxy R (0.95-quantile diet size)", ylabel="Δ (BAM − AM) at f*=$(fstar)",
               title="Redundancy threshold — $(String(geom))")
    scatterlines!(ax, RΔ.R, RΔ.Δ)
    hlines!(ax, [0.0], color=:gray, linestyle=:dash)
    save(joinpath(outdir, "FigC_redundancy_threshold.png"), fig); fig
end

# -- Figure D: Diet-size bins × geometry (k=1,2,≥3) ---------------------------
function figD_diet_bins(pool, grid, pars; fstar, geoms, A_fn)
    fig = Figure(; size=(900,320))
    for (j,g) in enumerate(geoms)
        res = penalty_vs_diet_bins(rng, pool, grid, pars; fstar=fstar, geom=g, A_fn=A_fn)
        ax  = Axis(fig[1,j], title=String(g),
                ylabel=j==1 ? "mean (BAM−AM)" : "", xlabel="diet-size bin")
        valid = .!isnan.(res.mean)
        x = 1:sum(valid)
        barplot!(ax, x, res.mean[valid])
        ax.xticks = (x, res.kbin[valid])
        hlines!(ax, [0.0], color=:gray, linestyle=:dash)
    end
    save(joinpath(outdir, "FigD_diet_bins.png"), fig); fig
end

# -- Figure E: Co-retention ECDF (intuitive biotic driver) --------------------
function figE_coretention_ecdf(rng, pool, grid, pars; fstar, geoms, A_fn)
    EC = coretention_ecdf(rng, pool, grid, pars; f=fstar, geoms=geoms, A_fn=A_fn)
    fig = Figure(; size=(520,360))
    ax  = Axis(fig[1,1], xlabel="K = # prey present (AM) in eligible cell", ylabel="ECDF",
               title="Prey co-retention at f*=$(fstar)")
    for g in geoms
        xs, ys = EC[g]
        isempty(xs) && continue
        lines!(ax, xs, ys; label=String(g))
    end
    axislegend(ax, position=:lt)
    save(joinpath(outdir, "FigE_coretention_ecdf.png"), fig); fig
end

# ============================================================================ #
#                                  RUN IT                                      #
# ============================================================================ #

@info "Running clean geometry analyses…"
grid        = grid_ridge
# Choose your showcase pool (mid is a good default)
pool = pool_low

# A) Prey sufficiency vs f
figA = figA_prey_sufficiency_vs_f(rng, pool, grid, pars; loss_fracs, geoms, A_fn, kreq=2)

# B) Biotic geometry effect at f*
figB = figB_BGE_bars(rng, pool, grid, pars; fstar=fstar, A_fn=A_fn)

# C) Redundancy threshold curve (single geometry = random)
figC = figC_redundancy_threshold(rng, grid, pars; fstar=fstar, geom=:front, A_fn=A_fn)

# D) Diet-size bins × geometry at f*
figD = figD_diet_bins(pool, grid, pars; fstar=fstar, geoms=geoms, A_fn=A_fn)

# E) Co-retention ECDF at f*
figE = figE_coretention_ecdf(rng, pool, grid, pars; fstar=fstar, geoms=geoms, A_fn=A_fn)

@info "All clean figures written to $(outdir)"
# ============================================================================ #
