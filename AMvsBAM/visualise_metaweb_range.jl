######## visualise_metaweb_range.jl ########################################
using CairoMakie, Statistics
include("../SetUp.jl")
include("src/metaweb.jl");   using .MetaWeb

# --- helpers --------------------------------------------------------------
function _adjacency_heatmap!(ax, mw)
    # sort species: basal first, then by in-degree to reveal structure
    roles = mw.trophic_role
    indeg = vec(sum(mw.A; dims=1))
    order = sortperm(1:length(roles); by=i -> (roles[i] == :basal ? 0 : 1, indeg[i]))
    Aord = mw.A[order, order]
    heatmap!(ax, Aord; interpolate=false)
    hidespines!(ax); ax.xlabel = "prey"; ax.ylabel = "predators"
end

function _degree_histogram!(ax, mw; normalize=true)
    A = mw.A
    indeg  = vec(sum(A; dims=1))      # consumers' diet breadth
    outdeg = vec(sum(A; dims=2))      # prey's vulnerability
    k = indeg .+ outdeg               # total degree
    # histogram of total degree (you can add subpanels if you prefer split)
    bins = 0:maximum(k)
    counts = zeros(Float64, length(bins))
    for v in k
        counts[Int(v)+1] += 1
    end
    if normalize
        counts ./= sum(counts)
        ax.ylabel = "frequency"
    else
        ax.ylabel = "count"
    end
    barplot!(ax, collect(bins), counts)
    ax.xlabel = "degree"
    hidespines!(ax)
end

# --- main entry -----------------------------------------------------------
"""
visualise_metaweb_range(; S=175, basal_frac=0.30, R95=5,
                         motif_types=[:chains,:mixed,:omnivory],
                         Cs=[0.04, 0.10, 0.18], seed=1,
                         figdir=joinpath(@__DIR__, "figs"))

Builds 3×3 examples of metawebs across motif types (rows) and connectances (columns),
and saves two figures:
  - MW_Adjacency_3x3.png  (adjacency matrices)
  - MW_Degree_3x3.png     (degree distributions)
"""
function visualise_metaweb_range(; S=175, basal_frac=0.30, R95=5,
                                 motif_types=[:chains,:mixed,:omnivory],
                                 Cs=[0.04, 0.10, 0.18], seed=1,
                                 figdir=joinpath(@__DIR__, "figs"))
    mkpath(figdir)
    rng = MersenneTwister(seed)

    # ---- 3×3 adjacency ----
    figA = Figure(size=(1200, 1000))
    for (i, motif) in enumerate(motif_types), (j, C) in enumerate(Cs)
        mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                                   connectance=C, R95=R95, motif_mix=motif)
        ax = Axis(figA[i,j], title="$(String(motif)) — C=$(round(C,digits=2))")
        _adjacency_heatmap!(ax, mw)
    end
    save(joinpath(figdir, "MW_Adjacency_3x3.png"), figA)
    display(figA)

    # ---- 3×3 degree distributions ----
    figB = Figure(size=(1200, 1000))
    for (i, motif) in enumerate(motif_types), (j, C) in enumerate(Cs)
        mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                                   connectance=C, R95=R95, motif_mix=motif)
        ax = Axis(figB[i,j], title="$(String(motif)) — C=$(round(C,digits=2))")
        _degree_histogram!(ax, mw; normalize=true)
    end
    save(joinpath(figdir, "MW_Degree_3x3.png"), figB)
    display(figB)

    return (; fig_adj=figA, fig_deg=figB)
end

# run if called directly
visualise_metaweb_range()
############################################################################
