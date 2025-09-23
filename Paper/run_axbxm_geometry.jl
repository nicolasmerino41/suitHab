# ==========================================================
# run_axbxm_geometry_makie.jl
# Geometry-sensitivity across A×B×M "realities" (Makie only)
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
basal_frac    = 0.45
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

# A×B×M “realities” (weights sum to 1)
weights = NTuple{3,Float64}[
    (1.0, 0.0, 0.0),  # A-led
    (0.0, 1.0, 0.0),  # B-led
    (0.0, 0.0, 1.0),  # M-led
    (0.5, 0.5, 0.0),  # A+B
    (0.5, 0.0, 0.5),  # A+M
    (0.0, 0.5, 0.5),  # B+M
]
weight_labels = ["A","B","M","A+B","A+M","B+M"]

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

# run one setting (grid, align, T) via your bootstrap helper
function run_setting(; rng, pool, grid::Grids.Grid, align::Float64, T::Int)
    pars = BSH.BAMParams(τA=0.5, τB=0.35, movement=:component, T=T)
    A_fn = make_A_fn(align)
    Metrics.bootstrap_mixture_maps(; rng, pool, grid, pars,
        weights=weights, fstar=fstar, A_fn=A_fn, agg=agg, kreq=kreq,
        geometries=(:random, :clustered, :front),
        A_seeds=A_seeds, pool_seeds=pool_seeds)
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
# bar + vertical error bars helper
function bars_with_err!(ax, y::AbstractVector, lo::AbstractVector, hi::AbstractVector;
                        xtlbls::Vector{String}=String[], title::AbstractString="")
    n  = length(y)
    xs = collect(1:n)

    # bars
    Mke.barplot!(ax, xs, y)

    # vertical error bars using quantile bounds
    verrorbars!(ax, xs, y, lo, hi; whisker=0.35, color=:black, lw=1.2)

    # zero line + cosmetics
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
        # fig.suptitle = @sprintf("%s — Geometry contrast vs A×B×M realities (ΔRF)", gname)
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
        
        Label(
            fig[0, :],
            @sprintf("%s — Geometry contrast vs A×B×M realities (ΔRC)", gname);
        )
        display(fig)
        Mke.save("$(saveprefix)_ΔRC_panels.png", fig)
    end
end

function plot_slices_for_grid(gname::String; saveprefix::String)
    idxA  = findfirst(==("A"),   weight_labels)
    idxAB = findfirst(==("A+B"), weight_labels)
    idxB  = findfirst(==("B"),   weight_labels)
    idxAM = findfirst(==("A+M"), weight_labels)
    idxM  = findfirst(==("M"),   weight_labels)
    xs = [0.0, 0.5, 1.0]
    xtlbls = ["0","0.5","1"]

    for Δname in ("ΔRF","ΔRC")
        fig = Mke.Figure(; size=(1100, 450))
        for (col, T) in enumerate(T_set)
            ax = Mke.Axis(fig[1, col], xlabel="weight along axis", title=@sprintf("%s slices (T=%d)", Δname, T))
            Mke.hlines!(ax, 0.0; color=:gray, linestyle=:dash)
            for align in align_set
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
                Mke.lines!(ax, xs, yB; color=:blue, label=@sprintf("B-axis, align=%.1f", align))
                verrorbars!(ax, xs, yB, loB, hiB; color=:blue, whisker=0.08, lw=1.0)

                # M-axis
                Mke.lines!(ax, xs, yM; color=:orange, linestyle=:dash, label=@sprintf("M-axis, align=%.1f", align))
                verrorbars!(ax, xs, yM, loM, hiM; color=:orange, whisker=0.08, lw=1.0)


            end
            ax.xticks = (xs, xtlbls)
            Mke.axislegend(ax, position=:rt, framevisible=false)
        end
        Label(
            fig[0, :], @sprintf("%s — %s slices across A×B×M", gname, Δname);
            fontsize=22, tellwidth=false
        )

        display(fig)
        Mke.save("$(saveprefix)_$(Δname)_slices.png", fig)
    end
end

# Draw vertical error bars at (x, y) spanning [lo, hi]
function verrorbars!(ax, xs::AbstractVector, y::AbstractVector,
                     lo::AbstractVector, hi::AbstractVector;
                     whisker::Float64 = 0.25, color=:black, lw::Real=1.5)
    @assert length(xs)==length(y)==length(lo)==length(hi)
    for (x, yi, l, h) in zip(xs, y, lo, hi)
        # main vertical segment
        Mke.lines!(ax, [x, x], [l, h]; color=color, linewidth=lw)
        # whiskers at ends
        Mke.lines!(ax, [x - whisker/2, x + whisker/2], [l, l]; color=color, linewidth=lw)
        Mke.lines!(ax, [x - whisker/2, x + whisker/2], [h, h]; color=color, linewidth=lw)
    end
    return nothing
end

# ---------------------------
# PRODUCE & SAVE FIGURES
# ---------------------------
for (_, gname) in settings
    saveprefix = replace(lowercase(gname), " "=>"_")
    plot_panels_for_grid(gname; saveprefix=saveprefix)
    plot_slices_for_grid(gname; saveprefix=saveprefix)
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
