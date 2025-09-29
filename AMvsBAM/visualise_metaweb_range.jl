######## visualise_metaweb_range.jl ########################################
using CairoMakie, Statistics
include("../SetUp.jl")
include("src/metaweb.jl");   using .MetaWeb

# --- helpers --------------------------------------------------------------
# ---------- consistent connectance metrics ----------
function connectance_metrics(mw)
    A      = mw.A
    S      = size(A,1)
    iscons = mw.trophic_role .!= :basal
    consix = findall(iscons)
    basix  = findall(.!iscons)

    links_total = sum(A)
    C_realized  = links_total / (S * max(length(consix), 1))

    # split into consumer→basal and consumer→consumer parts
    # rows = predators, cols = prey (A[row, col] = 1 means row eats col)
    # If your generator sometimes flips orientation, this is still fine
    # because we use biological roles to index blocks.
    C_cb = length(consix) == 0 || length(basix) == 0 ? 0.0 :
           sum(A[consix, basix]) / (length(consix) * length(basix))
    C_cc = length(consix) < 2 ? 0.0 :
           sum(A[consix, consix]) / (length(consix) * length(consix))

    return (; C_realized, C_cb, C_cc)
end

# ---------- adjacency heatmap with numbers in title ----------
function adjacency_panel!(ax, mw; clim=(0,1))
    A      = mw.A
    iscons = mw.trophic_role .!= :basal
    consix = findall(iscons)
    basix  = findall(.!iscons)

    # reorder so predators (consumers) are rows; prey (basal+consumers) are columns
    # rows: all species, but we will show consumers top-most for clarity
    row_order = vcat(consix, basix)                  # consumers first (predators)
    col_order = vcat(basix, consix)                  # basal first, then consumer prey

    Ar = A[row_order, col_order]                     # reordered adjacency

    image!(ax, float.(Ar); interpolate=false, colormap=:viridis, colorrange=clim)
    ax.xlabel = "prey"
    ax.ylabel = "predators"
    ax.xticksvisible = false
    ax.yticksvisible = false
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax
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

    figA = Figure(size=(1200, 1000))
    for (i, motif) in enumerate(motif_types), (j, Ctar) in enumerate(Cs)
        mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                                   connectance=Ctar, R95=R95, motif_mix=motif)
        # (optional) call your calibrator here if you want exact targets:
        # calibrate_connectance!(mw; C_target=Ctar, tol=0.001)

        mets = connectance_metrics(mw)
        ax = Axis(figA[i, j], titlesize=10)
        adjacency_panel!(ax, mw)
        ax.title = string(motif, " — C_target=", round(Ctar, digits=2),
                          " | realized=", round(mets.C_realized, digits=3),
                          " | cb=", round(mets.C_cb, digits=3),
                          " | cc=", round(mets.C_cc, digits=3))
    end
    save(joinpath(figdir, "MW_Adjacency_3x3_with_numbers.png"), figA)
    display(figA)
    return figA
end

# run if called directly
visualise_metaweb_range()
############################################################################

