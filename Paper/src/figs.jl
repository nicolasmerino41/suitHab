module Figs

using CairoMakie
using Statistics, Printf

export fig1_core, fig2_cdfs, fig3_rank_smallmultiples, fig4_regimemap
export fig2_metric_numbers
export fig_abs_bsh_vs_loss
export fig_realities_Aonly

palette_geom() = Dict(:random => RGBf(0.26,0.47,0.96),
                      :clustered => RGBf(0.96,0.59,0.12),
                      :front => RGBf(0.13,0.69,0.47))

# Fig.1 — AM vs BAM curves by geometry (+ optional placebo) + P_fail
function fig1_core(relcurves::Dict{Symbol,Any},
                   pfail::Dict{Symbol,Any};
                   placebo::Union{Nothing,Dict{Symbol,Any}}=nothing,
                   with_placebo::Bool=true,
                   title::String="")
    pal = palette_geom()
    fig = Figure(; size=(1600,560))
    Label(fig[0,1:2], title, fontsize=18, padding=(0,0,8,0))

    # left: AM vs BAM by geometry
    axL = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Relative ΔBSH (mean consumers)",
               title="AM (dashed) vs BAM (solid)")
    for g in (:random,:clustered,:front)
        r = relcurves[g]; p = pal[g]
        lines!(axL, r.x, r.relAM;  color=p, linestyle=:dash,  label="AM $g")
        lines!(axL, r.x, r.relBAM; color=p, linewidth=2,      label="BAM $g")
        if with_placebo && placebo !== nothing
            lp = placebo[g]
            lines!(axL, lp.x, lp.rel; color=p, linestyle=:dot, linewidth=2, label="placebo $g")
        end
    end
    axislegend(axL; position=:lb, framevisible=false)

    # right: P_fail curves (BAM)
    axR = Axis(fig[1,2], xlabel="Area lost (fraction)", ylabel="P_fail (A-suitable kept cells)",
               title="Biotic gate failure")
    for g in (:random,:clustered,:front)
        r = pfail[g]
        lines!(axR, r.x, r.y; color=pal[g], linewidth=2, label=string(g))
    end
    axislegend(axR; position=:rt, framevisible=false)
    fig
end

# Fig.2 — CDFs of per-species relative loss at f*
function fig2_cdfs(cdfs::Dict{Symbol,Any}; fstar::Float64=0.6, title::String="")
    pal = palette_geom()
    fig = Figure(; size=(900,560))
    Label(fig[0,1], title * @sprintf("  (f* = %.2f)", fstar),
          fontsize=18, padding=(0,0,8,0))
    ax = Axis(fig[1,1],
              xlabel="relative loss (−ΔBSH/BSH₀)",
              ylabel="Fraction of consumers",
              title="Per-species loss CDFs")

    for g in (:random,:clustered,:front)
        d = cdfs[g]
        stairs!(ax, d.xAM,  d.yAM;  color=pal[g], linestyle=:dash, linewidth=2,
                label=string(g)*" AM")
        stairs!(ax, d.xBAM, d.yBAM; color=pal[g], linewidth=2,
                label=string(g)*" BAM")
    end

    axislegend(ax; position=:rb, framevisible=false)
    return fig
end

# Fig 2.1
"Show KS (row 1) and tail-gap (row 2) per geometry as big numbers."
function fig2_metric_numbers(metrics_by_geom::Dict{Symbol,NamedTuple};
                             fstar::Real, tail::Real=0.8, title::AbstractString="")
    geoms = (:random, :clustered, :front)
    fig = Figure(; size=(1800, 650))
    Label(fig[0,1], title == "" ?
        "Dist. metrics at f* = $(round(fstar,digits=2)) (tail=$(round(tail,digits=2)))" : title,
        fontsize=24, tellwidth=false)

    for (j,g) in enumerate(geoms)
        ks   = get(metrics_by_geom, g, (KS=NaN, tail_gap=NaN)).KS
        tgap = get(metrics_by_geom, g, (KS=NaN, tail_gap=NaN)).tail_gap

        ax1 = Axis(fig[1,j]); hidespines!(ax1); hidedecorations!(ax1)
        text!(ax1, 0.5, 0.5, text=isnan(ks) ? "NA" : @sprintf("%.3f", ks),
              align=(:center,:center), fontsize=64)

        ax2 = Axis(fig[2,j]); hidespines!(ax2); hidedecorations!(ax2)
        text!(ax2, 0.5, 0.5, text=isnan(tgap) ? "NA" : @sprintf("%.3f", tgap),
              align=(:center,:center), fontsize=64)
    end
    fig
end

# Fig.3 — worst-geometry markers (small multiples, flexible grid)
function fig3_rank_smallmultiples(ranks::Vector{NamedTuple}; 
                                  title="Worst geometry at f*",
                                  ncols::Int=4)
    pal = palette_geom()
    n = length(ranks)
    nrows = cld(n, ncols)  # number of rows needed

    fig = Figure(; size=(300*ncols, 200*nrows))  # scale size to number of panels
    Label(fig[0,1:ncols], title, fontsize=18, padding=(0,0,8,0))

    for (k,r) in enumerate(ranks)
        i = cld(k, ncols)       # row index
        j = ((k-1) % ncols) + 1 # column index
        ax = Axis(fig[i,j], title=r.label, xgridvisible=false, ygridvisible=false)
        hidedecorations!(ax)

        text!(ax, 0.1, 0.6, text="AM :", color=:gray, fontsize=14, align=(:left,:center))
        scatter!(ax, [0.6], [0.6]; color=pal[r.am], markersize=16)

        text!(ax, 0.1, 0.3, text="BAM:", color=:gray, fontsize=14, align=(:left,:center))
        scatter!(ax, [0.6], [0.3]; color=pal[r.bam], markersize=16)
    end
    return fig
end

# Fig.4 — regime maps (heatmaps) over D and R
function fig4_regimemap(Dgrid::Vector{Float64}, Rweb::Vector{Float64},
                        Zdid::Matrix{Float64}, Zflip::Matrix{Float64};
                        title::String="")
    fig = Figure(; size=(1200,520))
    Label(fig[0,1:2], title, fontsize=18, padding=(0,0,8,0))

    # --- panel 1
    ax1 = Axis(fig[1,1], xlabel="Climate tail index D", ylabel="Redundancy R",
               title="Mean DiD at f*")
    hm1 = heatmap!(ax1, Dgrid, Rweb, Zdid'; colormap=:viridis, interpolate=false)
    Colorbar(fig[1,2], hm1, width=12)

    # --- panel 2
    ax2 = Axis(fig[1,3], xlabel="Climate tail index D", ylabel="Redundancy R",
               title="Rank-flip probability")
    hm2 = heatmap!(ax2, Dgrid, Rweb, Zflip'; colormap=:plasma, interpolate=false)
    Colorbar(fig[1,4], hm2, width=12)

    return fig
end

"Plot absolute BSH curves (AM dashed, BAM solid) for several geometries."
function fig_abs_bsh_vs_loss(curves::Dict{Symbol,NamedTuple}; title::AbstractString="")
    geoms = collect(keys(curves))
    fig = Figure(; size=(1050, 340))
    for (j,g) in enumerate(geoms)
        dat = curves[g]
        ax = Axis(fig[1,j], title=String(g),
                  xlabel="area lost (fraction)",
                  ylabel="BSH (mean over consumers / original area)")
        lines!(ax, dat.loss, dat.AM;  color=:gray35, linestyle=:dash,  linewidth=3, label="AM")
        lines!(ax, dat.loss, dat.BAM; color=:black,  linestyle=:solid, linewidth=3, label="BAM")
        hlines!(ax, [0.0], color=(:gray,0.4), linestyle=:dot)
        axislegend(ax, position=:lb)
    end
    isempty(title) || Label(fig[0, 1:length(geoms)], title, fontsize=18, tellwidth=false)
    return fig
end

function fig_realities_Aonly(curves::Dict{Symbol,NamedTuple}; title::AbstractString="")
    geoms = collect(keys(curves))
    fig = Figure(; size=(1050,340))
    for (j,g) in enumerate(geoms)
        dat = curves[g]
        ax = Axis(fig[1,j], title=String(g),
                  xlabel="area lost (fraction)",
                  ylabel="A-only (mean over consumers / original area)")
        lines!(ax, dat.loss, dat.ABM; color=:dodgerblue, linewidth=3, label="ABM")
        lines!(ax, dat.loss, dat.MAB; color=:darkorange, linewidth=3, label="MAB")
        lines!(ax, dat.loss, dat.BAM; color=:seagreen,  linewidth=3, label="BAM")
        axislegend(ax, position=:lb)
    end
    isempty(title) || Label(fig[0,1:length(geoms)], title, fontsize=18, tellwidth=false)
    return fig
end

end # module