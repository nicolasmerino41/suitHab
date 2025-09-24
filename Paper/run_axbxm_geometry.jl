# ==========================================================
# run_axbxm_geometry_makie.jl
# Geometry-sensitivity across A×B×M "realities" (Makie only)
# Adds BM(noA): Biotic + Movement with NO abiotic gate.
# Figures are saved as PNGs in the working directory.
# ==========================================================
using Random, Statistics, Printf
using CairoMakie
const Mke = CairoMakie

include("src/grids.jl");     using .Grids
include("src/hl.jl");        using .HL
include("src/metawebs.jl");  using .Metawebs
include("src/bsh.jl");       using .BSH
include("src/metrics.jl");   using .Metrics

# ---------------------------
# CONFIG
# ---------------------------
rng           = MersenneTwister(7)
S             = 175
basal_frac    = 0.25
archetype     = :mid
nx, ny        = 30, 30

# two climate backdrops
grid_gradient = Grids.make_grid_gradient(nx, ny; seed=42)
grid_ridge    = Grids.make_grid_ridge(nx, ny; seed=45)

# HL loss level to compare geometries
fstar         = 0.60

# movement thresholds to test (scale-match vs geometry)
T_set         = [4, 16]    # "small components" vs "large components"

# alignment off/on for A-construction
align_set     = [0.0, 0.8]

# B-gate choice
agg           = :mean
kreq          = 1

# bootstrap seeds for A (and optionally pool)
A_seeds       = 1:24
pool_seeds    = nothing  # set to e.g., 1:10 to also vary the pool

# A×B×M mixtures for the panels/slices
weights = NTuple{3,Float64}[
    (1.0, 0.0, 0.0),  # A-led
    (0.0, 1.0, 0.0),  # B-led (BAM)
    (0.0, 0.0, 1.0),  # M-led (BAM)
    (0.5, 0.5, 0.0),  # A+B
    (0.5, 0.0, 0.5),  # A+M
    (0.0, 0.5, 0.5),  # B+M (with A)
]
weight_labels = ["A","B","M","A+B","A+M","B+M","BM(noA)"]  # <- NEW label added at the end

# ---------------------------
# BUILD POOL
# ---------------------------
pool = Metawebs.build_metaweb_archetype(rng; S=S, basal_frac=basal_frac, archetype=archetype)
@info "Pool diagnostics" Metawebs.metaweb_diagnostics(pool)

# ---------------------------
# helpers
# ---------------------------

# A-constructor closure (alignment on/off), keep same basal bias you used
function make_A_fn(align_val::Float64)
    (p::Metawebs.SpeciesPool, g::Grids.Grid; seed::Int=1) ->
        BSH.abiotic_matrix_aligned(p, g;
            niche_basal=0.10, niche_cons=0.12,
            bias_basal=0.7, align=align_val, seed=seed)
end

# ---------- NEW: BM(noA) mechanics in this script (no changes to modules) ----------

"Movement component gate on KEPT cells only (no A). Returns a BitVector over cells."
function comp_gate_keep_only(grid::Grids.Grid, keep::BitVector, T::Int)
    C = grid.C; nx, ny = grid.nx, grid.ny
    M = falses(C)
    seen = falses(C)
    for i in 1:C
        (keep[i] && !seen[i]) || continue
        comp = Int[]; q = [i]; seen[i] = true
        while !isempty(q)
            v = popfirst!(q); push!(comp, v)
            for nb in HL.neighbors4(v, nx, ny)
                if keep[nb] && !seen[nb]
                    seen[nb] = true; push!(q, nb)
                end
            end
        end
        if length(comp) ≥ T
            for v in comp; M[v] = true; end
        end
    end
    return M
end

"Assembly with Biotic + Movement only (NO abiotic suitability)."
function assemble_BM_noA(pool::Metawebs.SpeciesPool, grid::Grids.Grid,
                         keep::BitVector; τB::Float64, T::Int)
    C = grid.C; S = pool.S
    Mgate = comp_gate_keep_only(grid, keep, T)
    P = falses(S, C)
    preyfrac = zeros(Float64, S, C)

    order = sortperm(pool.masses)  # bottom-up
    for s in order
        if pool.basal[s]
            @inbounds for i in 1:C
                P[s,i] = keep[i] & Mgate[i]   # no A
            end
        else
            pr = pool.prey[s]
            if isempty(pr)
                @inbounds for i in 1:C; P[s,i] = false; end
            else
                @inbounds for i in 1:C
                    if keep[i] & Mgate[i]
                        m = mean(@view P[pr,i])
                        preyfrac[s,i] = m
                        P[s,i] = (m ≥ τB)
                    else
                        P[s,i] = false
                    end
                end
            end
        end
    end
    bsh = [sum(@view P[s,:]) / C for s in 1:S]
    return (; P, bsh, preyfrac)
end

"Mean over consumers of retained area (BM with NO A)."
function mean_BM_noA_over_consumers(pool::Metawebs.SpeciesPool, grid::Grids.Grid,
                                    keep::BitVector; τB::Float64, T::Int)
    C = grid.C
    if count(keep) == 0; return 0.0; end
    bm = assemble_BM_noA(pool, grid, keep; τB=τB, T=T)
    return mean(bm.bsh[.!pool.basal])
end

"Bootstrapped ΔRF/ΔRC for BM(noA) at f*; mirrors what bootstrap_mixture_maps returns for mixtures."
function bmnoa_contrasts_bootstrap(; rng, pool, grid::Grids.Grid, τB::Float64, T::Int,
                                   fstar::Float64, geometries::Tuple=(:random,:clustered,:front),
                                   A_seeds=1:24, pool_seeds=nothing)
    nx, ny = grid.nx, grid.ny
    keepmask(g, f) = g===:random    ? HL.random_mask(rng, grid.C, 1-f) :
                      g===:clustered ? HL.clustered_mask(rng, nx, ny, 1-f; nseeds=8) :
                                       HL.front_mask(rng, grid.xy, 1-f; axis=:x, noise=0.05)

    yR = Float64[]; yC = Float64[]; yF = Float64[]
    nrep = length(A_seeds)
    for rep in 1:nrep
        p = if pool_seeds === nothing
            pool
        else
            seed = pool_seeds[ ((rep-1) % length(pool_seeds)) + 1 ]
            Metawebs.build_metaweb_archetype(MersenneTwister(seed);
                S=pool.S, basal_frac=mean(pool.basal), archetype=:mid)
        end
        kR = keepmask(:random,    fstar)
        kC = keepmask(:clustered, fstar)
        kF = keepmask(:front,     fstar)
        push!(yR, mean_BM_noA_over_consumers(p, grid, kR; τB=τB, T=T))
        push!(yC, mean_BM_noA_over_consumers(p, grid, kC; τB=τB, T=T))
        push!(yF, mean_BM_noA_over_consumers(p, grid, kF; τB=τB, T=T))
    end
    # contrasts
    dRF = yR .- yF
    dRC = yR .- yC
    summarize(v) = (mean(v), quantile(v,0.10), quantile(v,0.90))
    μRF, loRF, hiRF = summarize(dRF)
    μRC, loRC, hiRC = summarize(dRC)
    return (; μRF, loRF, hiRF, μRC, loRC, hiRC)
end

# ---------- end of BM(noA) helpers ----------

# run one setting (grid, align, T): call your bootstrap, then append BM(noA)
function run_setting(; rng, pool, grid::Grids.Grid, align::Float64, T::Int)
    pars = BSH.BAMParams(τA=0.5, τB=0.35, movement=:component, T=T)
    A_fn = make_A_fn(align)
    base = Metrics.bootstrap_mixture_maps(; rng, pool, grid, pars,
        weights=weights, fstar=fstar, A_fn=A_fn, agg=agg, kreq=kreq,
        geometries=(:random, :clustered, :front),
        A_seeds=A_seeds, pool_seeds=pool_seeds)

    # BM(noA) bootstrap at f*
    extra = bmnoa_contrasts_bootstrap(; rng, pool, grid, τB=pars.τB, T=T,
                                      fstar=fstar, geometries=(:random,:clustered,:front),
                                      A_seeds=A_seeds, pool_seeds=pool_seeds)

    # Append to ΔRF/ΔRC vectors (kept as the 7th bar)
    dRF_mean = vcat(base.dRF_mean, extra.μRF)
    dRF_lo   = vcat(base.dRF_lo,   extra.loRF)
    dRF_hi   = vcat(base.dRF_hi,   extra.hiRF)

    dRC_mean = vcat(base.dRC_mean, extra.μRC)
    dRC_lo   = vcat(base.dRC_lo,   extra.loRC)
    dRC_hi   = vcat(base.dRC_hi,   extra.hiRC)

    return (; base.elastic_mean, base.elastic_lo, base.elastic_hi,
            dRF_mean, dRF_lo, dRF_hi, dRC_mean, dRC_lo, dRC_hi)
end

# geometry sensitivity index (optional)
gsi(ΔRF::Vector{<:Real}, ΔRC::Vector{<:Real}) = [max(abs(ΔRF[i]), abs(ΔRC[i])) for i in eachindex(ΔRF)]

# ---------------------------
# RUN: all combinations
# ---------------------------
settings = [(grid_gradient, "Gradient"), (grid_ridge, "Ridge")]
results  = Dict{String,Any}()

for (grid, gname) in settings
    for T in T_set, align in align_set
        key = @sprintf("%s | T=%d | align=%.1f", gname, T, align)
        @info "Running $key ..."
        results[key] = run_setting(; rng, pool, grid, align, T)
    end
end

# ---------------------------
# PLOTTING (Makie; required layout structure)
# ---------------------------
# Draw vertical error bars at (x, y) spanning [lo, hi]
function verrorbars!(ax, xs::AbstractVector, y::AbstractVector,
                     lo::AbstractVector, hi::AbstractVector;
                     whisker::Float64 = 0.25, color=:black, lw::Real=1.5)
    @assert length(xs)==length(y)==length(lo)==length(hi)
    for (x, yi, l, h) in zip(xs, y, lo, hi)
        Mke.lines!(ax, [x, x], [l, h]; color=color, linewidth=lw)          # main segment
        Mke.lines!(ax, [x - whisker/2, x + whisker/2], [l, l]; color=color, linewidth=lw) # whisker low
        Mke.lines!(ax, [x - whisker/2, x + whisker/2], [h, h]; color=color, linewidth=lw) # whisker high
    end
    return nothing
end

# bar + vertical error bars helper
function bars_with_err!(ax, y::AbstractVector, lo::AbstractVector, hi::AbstractVector;
                        xtlbls::Vector{String}=String[], title::AbstractString="")
    n  = length(y)
    xs = collect(1:n)
    Mke.barplot!(ax, xs, y)
    verrorbars!(ax, xs, y, lo, hi; whisker=0.35, color=:black, lw=1.2)
    Mke.hlines!(ax, 0.0; color=:gray, linestyle=:dash, linewidth=1)
    !isempty(xtlbls) && (ax.xticks = (xs, xtlbls))
    ax.title = title
    return ax
end

function plot_panels_for_grid(gname::String; saveprefix::String)
    # 2×2: rows = align {off,on}, cols = T {small,large}
    begin
        fig = Mke.Figure(; size=(1200, 700))
        # ΔRF
        for (col, T) in enumerate(T_set)
            for (row, align) in enumerate(align_set)
                key = @sprintf("%s | T=%d | align=%.1f", gname, T, align)
                res = results[key]
                ax = Mke.Axis(fig[row, col], ylabel = row==1 ? "ΔRF = Y(random)-Y(front)" : "")
                bars_with_err!(ax, res.dRF_mean, res.dRF_lo, res.dRF_hi;
                               xtlbls=weight_labels,
                               title=@sprintf("ΔRF  (T=%d, align=%.1f)", T, align))
            end
        end
        display(fig)
        Mke.save("$(saveprefix)_ΔRF_panels.png", fig)
    end

    begin
        fig = Mke.Figure(; size=(1200, 700))
        # ΔRC
        for (col, T) in enumerate(T_set)
            for (row, align) in enumerate(align_set)
                key = @sprintf("%s | T=%d | align=%.1f", gname, T, align)
                res = results[key]
                ax = Mke.Axis(fig[row, col], ylabel = row==1 ? "ΔRC = Y(random)-Y(clustered)" : "")
                bars_with_err!(ax, res.dRC_mean, res.dRC_lo, res.dRC_hi;
                               xtlbls=weight_labels,
                               title=@sprintf("ΔRC  (T=%d, align=%.1f)", T, align))
            end
        end
        display(fig)
        Mke.save("$(saveprefix)_ΔRC_panels.png", fig)
    end
end

function plot_slices_for_grid(gname::String; saveprefix::String)
    # indices on the 6-mixture set (BM(noA) is not on these axes)
    idxA  = 1
    idxAB = 4
    idxB  = 2
    idxAM = 5
    idxM  = 3
    xs = [0.0, 0.5, 1.0]
    xtlbls = ["0","0.5","1"]

    for Δname in ("ΔRF","ΔRC")
        begin
            fig = Mke.Figure(; size=(1100, 450))
            for (col, T) in enumerate(T_set)
                ax = Mke.Axis(fig[1, col], xlabel="weight along axis", title=@sprintf("%s slices (T=%d)", Δname, T))
                Mke.hlines!(ax, 0.0; color=:gray, linestyle=:dash)
                for (row, align) in enumerate(align_set)
                    key = @sprintf("%s | T=%d | align=%.1f", gname, T, align)
                    res = results[key]
                    if Δname == "ΔRF"
                        yB  = [res.dRF_mean[idxA], res.dRF_mean[idxAB], res.dRF_mean[idxB]]
                        loB = [res.dRF_lo[idxA],  res.dRF_lo[idxAB],  res.dRF_lo[idxB]]
                        hiB = [res.dRF_hi[idxA],  res.dRF_hi[idxAB],  res.dRF_hi[idxB]]
                        yM  = [res.dRF_mean[idxA], res.dRF_mean[idxAM], res.dRF_mean[idxM]]
                        loM = [res.dRF_lo[idxA],  res.dRF_lo[idxAM],  res.dRF_lo[idxM]]
                        hiM = [res.dRF_hi[idxA],  res.dRF_hi[idxAM],  res.dRF_hi[idxM]]
                    else
                        yB  = [res.dRC_mean[idxA], res.dRC_mean[idxAB], res.dRC_mean[idxB]]
                        loB = [res.dRC_lo[idxA],  res.dRC_lo[idxAB],  res.dRC_lo[idxB]]
                        hiB = [res.dRC_hi[idxA],  res.dRC_hi[idxAB],  res.dRC_hi[idxB]]
                        yM  = [res.dRC_mean[idxA], res.dRC_mean[idxAM], res.dRC_mean[idxM]]
                        loM = [res.dRC_lo[idxA],  res.dRC_lo[idxAM],  res.dRC_lo[idxM]]
                        hiM = [res.dRC_hi[idxA],  res.dRC_hi[idxAM],  res.dRC_hi[idxM]]
                    end

                    # B-axis
                    if row == 1
                        Mke.lines!(ax, xs, yB; color=:blue, label=@sprintf("B-axis, align=%.1f", align))
                        verrorbars!(ax, xs, yB, loB, hiB; color=:blue, whisker=0.08, lw=1.0)
                    else
                        Mke.lines!(ax, xs, yB; color=:blue, linestyle=:dash, label=@sprintf("B-axis, align=%.1f", align))
                        verrorbars!(ax, xs, yB, loB, hiB; color=:blue, whisker=0.08, lw=1.0)
                    end

                    # M-axis
                    if row == 1
                        Mke.lines!(ax, xs, yM; color=:orange, label=@sprintf("M-axis, align=%.1f", align))
                        verrorbars!(ax, xs, yM, loM, hiM; color=:orange, whisker=0.08, lw=1.0)
                    else
                        Mke.lines!(ax, xs, yM; color=:orange, linestyle=:dash, label=@sprintf("M-axis, align=%.1f", align))
                        verrorbars!(ax, xs, yM, loM, hiM; color=:orange, whisker=0.08, lw=1.0)
                    end
                end
                ax.xticks = (xs, xtlbls)
                Mke.axislegend(ax, position=:rt, framevisible=false)
            end
            display(fig)
            Mke.save("$(saveprefix)_$(Δname)_slices.png", fig)
        end
    end
end

# ---------------------------
# PRODUCE & SAVE FIGURES
# ---------------------------
for (_, gname) in settings
    saveprefix = replace(lowercase(gname), " "=>"_")
    plot_panels_for_grid(gname; saveprefix=saveprefix)
    plot_slices_for_grid(gname; saveprefix=saveprefix)  # unchanged (BM(noA) is off-axis)
end

println("\nSaved:")
println("  gradient_ΔRF_panels.png")
println("  gradient_ΔRC_panels.png")
println("  gradient_ΔRF_slices.png")
println("  gradient_ΔRC_slices.png")
println("  ridge_ΔRF_panels.png")
println("  ridge_ΔRC_panels.png")
println("  ridge_ΔRF_slices.png")
println("  ridge_ΔRC_slices.png")
