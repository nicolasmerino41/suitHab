# visualise_metaweb_range.jl
# Posterior visualizations of metawebs (R95 free, sweep connectance C)

using Random, Statistics, StatsBase, CairoMakie
using .MetaWeb

function _adjacency_heatmap!(ax, A::BitMatrix)
    image!(ax, permutedims(Float64.(A)); colormap=:viridis)
    ax.xlabel = "prey"; ax.ylabel = "predators"
end

# Histogram helper that drops zeros
function _hist_nozeros!(ax, vals; nbins=25, lab="")
    v = filter(!iszero, vals)
    isempty(v) && return
    h = fit(Histogram, v; nbins=nbins)
    xs = midpoints(h.edges[1])
    ys = h.weights ./ sum(h.weights)
    barplot!(ax, xs, ys; label=lab, transparency=true)
end

# Degree utilities
indeg(A)  = vec(sum(A; dims=1))   # how many predators eat q
outdeg(A) = vec(sum(A; dims=2))   # how many prey predator p eats

# Realized R95 (consumers only)
r95_from_A(A::AbstractMatrix, nb::Int) = begin
    ks = outdeg(A)
    quantile(ks[(nb+1):end], 0.95)
end

"""
    visualise_metaweb_range(; S=175, basal_frac=0.30,
                             motif_types=[:chains,:mixed,:omnivory],
                             Cs=[0.001, 0.01, 0.04, 0.10, 0.18],
                             seed=1, figdir="figs_metaweb")

Builds metawebs with **control = :C** (exact C), leaving **R95 free**.
Makes three figures (rows=motifs, cols=connectance):
1) adjacency matrices,
2) in- & out-degree histograms,
3) total degree histograms.
"""
function visualise_metaweb_range(; S=175, basal_frac=0.30,
                                 motif_types=[:chains,:mixed,:omnivory],
                                 Cs=[0.001, 0.01, 0.04, 0.10, 0.18],
                                 seed=1,
                                 figdir=joinpath(@__DIR__, "figs_metaweb"))
    mkpath(figdir)
    rng = MersenneTwister(seed)
    nb  = max(1, round(Int, basal_frac*S))

    nrows = length(motif_types); ncols = length(Cs)
    figsize = (260*ncols, 260*nrows)

    # --------- (1) Adjacencies ------------
    figA = Figure(size=figsize)
    for (i, motif) in enumerate(motif_types), (j, Ctar) in enumerate(Cs)
        mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                                   control=:C, connectance=Ctar,
                                   R95=nothing, motif_mix=motif)
        A = mw.A
        Creal  = MetaWeb.global_connectance(A)
        R95real = r95_from_A(A, nb)
        ax = Axis(figA[i,j],
            title="$(String(motif)) — target C=$(round(Ctar,digits=3)) | " *
                  "Cᵣ=$(round(Creal,digits=3)), R95ᵣ=$(round(R95real,digits=1))",
            titlesize=10)
        _adjacency_heatmap!(ax, A)
    end
    save(joinpath(figdir, "MW_Adjacency_Cs.png"), figA)
    display(figA)

    # --------- (2) In/Out degree -----------
    figB = Figure(size=figsize)
    for (i, motif) in enumerate(motif_types), (j, Ctar) in enumerate(Cs)
        mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                                   control=:C, connectance=Ctar,
                                   R95=nothing, motif_mix=motif)
        A = mw.A
        ax = Axis(figB[i,j],
                  title="$(String(motif)) — C=$(round(Ctar,digits=3))",
                  xlabel="degree", ylabel="frequency", titlesize=10)
        _hist_nozeros!(ax, indeg(A); nbins=25, lab="in-degree")
        _hist_nozeros!(ax, outdeg(A); nbins=25, lab="out-degree")
        axislegend(ax; position=:rt, framevisible=false, labelsize=10)
    end
    save(joinpath(figdir, "MW_InOutDegree_Cs.png"), figB)
    display(figB)

    # --------- (3) Total degree ------------
    figC = Figure(size=figsize)
    for (i, motif) in enumerate(motif_types), (j, Ctar) in enumerate(Cs)
        mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                                   control=:C, connectance=Ctar,
                                   R95=nothing, motif_mix=motif)
        A = mw.A
        d = indeg(A) .+ outdeg(A)
        ax = Axis(figC[i,j],
                  title="$(String(motif)) — C=$(round(Ctar,digits=3))",
                  xlabel="total degree", ylabel="frequency", titlesize=10)
        _hist_nozeros!(ax, d; nbins=25, lab="degree")
    end
    save(joinpath(figdir, "MW_TotalDegree_Cs.png"), figC)
    display(figC)

    return (; fig_adj=figA, fig_inout=figB, fig_deg=figC)
end

# If you want to run directly:
visualise_metaweb_range()
