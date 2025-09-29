# visualise_metaweb_range.jl
# Simple, *posterior* visualizations of the original metawebs we use.
function _adjacency_heatmap!(ax, A::BitMatrix)
    image!(ax, permutedims(Float64.(A)); colormap=:viridis)
    ax.xlabel = "prey"; ax.ylabel = "predators"
end

# Histogram helper that drops zeros
# Histogram helper that drops zeros
function _hist_nozeros!(ax, vals; nbins=25, lab="")
    v = filter(!iszero, vals)
    isempty(v) && return
    h = fit(Histogram, v; nbins=nbins)
    xs = midpoints(h.edges[1])
    ys = h.weights ./ sum(h.weights)
    barplot!(ax, xs, ys; label=lab, transparency=true)
end


# Compute in/out degree (full matrix, not just consumers)
indeg(A)  = vec(sum(A; dims=1))   # how many predators eat q
outdeg(A) = vec(sum(A; dims=2))   # how many prey predator p eats

"""
    visualise_metaweb_range(; S=175, basal_frac=0.30, R95=5,
                             motif_types=[:chains,:mixed,:omnivory],
                             Cs=[0.04, 0.10, 0.18], seed=1,
                             figdir="figs_metaweb")

Makes three figures (3×3 each, rows=motifs, cols=connectance):
1) adjacency matrices,
2) in- & out-degree histograms overlapped (zeros excluded),
3) total degree histograms (zeros excluded).
"""
function visualise_metaweb_range(; S=175, basal_frac=0.30, R95=5,
                                 motif_types=[:chains,:mixed,:omnivory],
                                 Cs=[0.001, 0.01,0.04, 0.1, 0.18], seed=1,
                                 figdir=joinpath(@__DIR__, "figs_metaweb"))
    mkpath(figdir)
    rng = MersenneTwister(seed)

    # --------- (1) Adjacencies ------------
    figA = Figure(size=(1200,1000))
    for (i, motif) in enumerate(motif_types), (j, Ctar) in enumerate(Cs)
        mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                                   connectance=Ctar, R95=R95, motif_mix=motif)
        A = mw.A
        Creal = count(==(true), A) / (S^2)
        ax = Axis(
            figA[i,j],
            title="$(String(motif)) — target C=$(Ctar) | realized C=$(round(Creal, digits=3))",
            titlesize=10
        )
        _adjacency_heatmap!(ax, A)
    end
    save(joinpath(figdir, "MW_Adjacency_3x3.png"), figA)
    display(figA)

    # --------- (2) In/Out degree -----------
    figB = Figure(size=(1200,1000))
    for (i, motif) in enumerate(motif_types), (j, Ctar) in enumerate(Cs)
        mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                                   connectance=Ctar, R95=R95, motif_mix=motif)
        A = mw.A
        ax = Axis(figB[i,j], title="$(String(motif)) — C=$(Ctar)",
                  xlabel="degree", ylabel="frequency",
                  titlesize=10)
        _hist_nozeros!(ax, indeg(A); nbins=25, lab="in-degree")
        _hist_nozeros!(ax, outdeg(A); nbins=25, lab="out-degree")
        axislegend(ax; position=:rt, framevisible=false, labelsize=10)
    end
    save(joinpath(figdir, "MW_InOutDegree_3x3.png"), figB)
    display(figB)

    # --------- (3) Total degree ------------
    figC = Figure(size=(1200,1000))
    for (i, motif) in enumerate(motif_types), (j, Ctar) in enumerate(Cs)
        mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                                   connectance=Ctar, R95=R95, motif_mix=motif)
        A = mw.A
        d = indeg(A) .+ outdeg(A)
        ax = Axis(figC[i,j], title="$(String(motif)) — C=$(Ctar)",
                  xlabel="total degree", ylabel="frequency",
                  titlesize=10)
        _hist_nozeros!(ax, d; nbins=25, lab="degree")
    end
    save(joinpath(figdir, "MW_TotalDegree_3x3.png"), figC)
    display(figC)

    return (; fig_adj=figA, fig_inout=figB, fig_deg=figC)
end

# If you want to run directly:
visualise_metaweb_range()
