module Figs

using CairoMakie
using Statistics, Printf

export fig1_core, fig2_cdfs, fig3_rank_smallmultiples, fig4_regimemap

palette_geom() = Dict(:random => RGBf(0.26,0.47,0.96),
                      :clustered => RGBf(0.96,0.59,0.12),
                      :front => RGBf(0.13,0.69,0.47))

# Fig.1 — AM vs BAM curves by geometry + placebo + P_fail
function fig1_core(relcurves::Dict{Symbol,Any}, placebo::Dict{Symbol,Any}, pfail::Dict{Symbol,Any};
                   title::String="")
    pal = palette_geom()
    fig = Figure(; size=(1600,560))
    Label(fig[0,1:2], title, fontsize=18, padding=(0,0,8,0))

    # left: AM vs BAM by geometry
    axL = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Relative ΔBSH (mean consumers)",
               title="AM (dashed) vs BAM (solid)")
    for g in (:random,:clustered,:front)
        r = relcurves[g]; p = pal[g]
        lines!(axL, r.x, r.relAM; color=p, linestyle=:dash, label="AM $g")
        lines!(axL, r.x, r.relBAM; color=p, linewidth=2, label="BAM $g")
        lp = placebo[g]
        lines!(axL, lp.x, lp.rel; color=p, linestyle=:dot, linewidth=2, label="placebo $g")
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
    fig = Figure(resolution=(900,560))
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

# Fig.3 — worst-geometry markers (small multiples)
function fig3_rank_smallmultiples(ranks::Vector{NamedTuple}; title="Worst geometry at f*")
    pal = palette_geom()
    fig = Figure(; size=(1200,600))
    Label(fig[0,1:3], title, fontsize=18, padding=(0,0,8,0))
    n = length(ranks)
    for (k,r) in enumerate(ranks)
        ax = Axis(fig[1, k], title=r.label, xgridvisible=false, ygridvisible=false)
        hidedecorations!(ax)
        text!(ax, 0.1, 0.6, text="AM :", color=:gray, fontsize=16, align=(:left,:center))
        scatter!(ax, [0.6], [0.6]; color=pal[r.am], markersize=18)
        text!(ax, 0.1, 0.3, text="BAM:", color=:gray, fontsize=16, align=(:left,:center))
        scatter!(ax, [0.6], [0.3]; color=pal[r.bam], markersize=18)
    end
    fig
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

end # module
