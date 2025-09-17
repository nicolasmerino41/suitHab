module FigSweeps

using CairoMakie, Printf
using ..BSH
using ..Metawebs
using ..Grids

export sweep_fig1_3x3

# -- helper: draw a compact Fig.1-style panel into an Axis
function _mini_panel!(ax::Axis, rel::Dict{Symbol,Any}, plac::Dict{Symbol,Any}; title::String="")
    geoms  = (:random, :clustered, :front)
    colors = (:dodgerblue, :darkorange, :seagreen)

    for (g, col) in zip(geoms, colors)
        r = rel[g]
        lines!(ax, r.x, r.relAM;  color=col, linestyle=:dash,  linewidth=2)  # AM
        lines!(ax, r.x, r.relBAM; color=col, linestyle=:solid, linewidth=3)  # BAM
        p = plac[g]
        lines!(ax, p.x, p.rel;    color=col, linestyle=:dot,   linewidth=2)  # placebo
    end
    ax.title   = title
    ax.xlabel  = "area lost (fraction)"
    ax.ylabel  = "relative ΔBSH"
    hidespines!(ax, :t, :r)
    return ax
end

"""
sweep_fig1_3x3(; rng, grid, archetype=:mid,
               tauB_list=[0.40,0.50,0.60], T_list=[6,8,12],
               S=175, basal_frac=0.35, loss_fracs=0.2:0.1:0.8, seed_A=1,
               τA=0.5, τocc=0.2, γ=3.0, movement::Symbol=:component,
               outpath="figs/fig1_sweep_mid.png")

Builds a 3×3 grid of compact Fig.1-style panels:
rows = τB (biotic threshold), columns = movement gate T (component size).

- Movement is ON when `movement=:component` (recommended for contrast).
- One call produces ONE figure for the chosen `archetype`.
"""
function sweep_fig1_3x3(; rng, archetype::Symbol=:mid,
    grids_vec::Vector{Tuple{String, Grid}}, grid_type::String="patchy",
    tauB_list=[0.40, 0.50, 0.60], T_list=[6, 8, 12],
    S::Int=175, basal_frac::Float64=0.35, loss_fracs=0.2:0.1:0.8, seed_A::Int=1,
    τA::Float64=0.5, τocc::Float64=0.2, γ::Float64=3.0,
    movement::Symbol=:component,
    outpath::AbstractString="figs/fig1_sweep_$(String(archetype)).png")

    # one metaweb per archetype
    pool = Metawebs.build_metaweb_archetype(rng; S=S, basal_frac=basal_frac, archetype=archetype)
    
    # pick the grid from the vector by matching the name
    idx = findfirst(x -> x[1] == grid_type, grids_vec)
    name, grid = grids_vec[idx]
    
    fig = Figure(; size=(1800, 1500))
    for (i, τB) in enumerate(tauB_list), (j, T) in enumerate(T_list)
        pars = BSH.BAMParams(; τA=τA, τB=τB, τocc=τocc, γ=γ, movement=movement, T=T)

        rel   = BSH.relative_loss_curves(rng, pool, grid, pars; loss_fracs, seed_A)
        plac  = BSH.placebo_curves(; rng, pool, grid, pars, loss_fracs, seed_A)

        ax = Axis(fig[i, j])
        _mini_panel!(ax, rel, plac; title=@sprintf("τB=%.2f,  T=%d", τB, T))
    end

    Label(fig[0, 1:3],
        "Archetype = $(String(archetype)) — movement $(movement) — $name grid",
        fontsize=20, tellwidth=false)

    save(outpath, fig)
    
    return fig
end

end # module
