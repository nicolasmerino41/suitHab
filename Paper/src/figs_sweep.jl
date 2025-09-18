module FigSweeps

using CairoMakie, Printf
using ..BSH
using ..Metawebs
using ..Grids

export sweep_fig1_3x3, sweep_fig1_pair_3x3, sweep_fig1_3x3_for_pool


# -- helper: draw a compact Fig.1-style panel into an Axis
function _mini_panel!(ax::Axis, rel::Dict{Symbol,Any}, plac::Union{Nothing,Dict{Symbol,Any}};
                      title::String="", with_placebo::Bool=true)
    geoms  = (:random, :clustered, :front)
    colors = (:dodgerblue, :darkorange, :seagreen)

    for (g, col) in zip(geoms, colors)
        r = rel[g]
        lines!(ax, r.x, r.relAM;  color=col, linestyle=:dash,  linewidth=2)  # AM
        lines!(ax, r.x, r.relBAM; color=col, linestyle=:solid, linewidth=3)  # BAM
        if with_placebo && plac !== nothing
            p = plac[g]
            lines!(ax, p.x, p.rel; color=col, linestyle=:dot, linewidth=2)   # placebo
        end
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
               outpath="figs/fig1_sweep_mid.png", with_placebo=true)

Builds a 3×3 grid of compact Fig.1-style panels:
rows = τB (biotic threshold), columns = movement gate T (component size).

- Movement is ON when `movement=:component` (recommended for contrast).
- One call produces ONE figure for the chosen `archetype`.
- If `with_placebo=false`, placebo curves are skipped.
"""
function sweep_fig1_3x3(; rng, archetype::Symbol=:mid,
    grids_vec::Vector{Tuple{String, Grid}}, grid_type::String="patchy",
    tauB_list=[0.40, 0.50, 0.60], T_list=[6, 8, 12],
    S::Int=175, basal_frac::Float64=0.35, loss_fracs=0.2:0.1:0.8, seed_A::Int=1,
    τA::Float64=0.5, τocc::Float64=0.2, γ::Float64=3.0,
    movement::Symbol=:component, with_placebo::Bool=true,
    outpath::AbstractString="figs/fig1_sweep_$(String(archetype)).png")

    # one metaweb per archetype
    pool = Metawebs.build_metaweb_archetype(rng; S=S, basal_frac=basal_frac, archetype=archetype)

    # pick the grid from the vector by matching the name
    idx = findfirst(x -> x[1] == grid_type, grids_vec)
    name, grid = grids_vec[idx]

    fig = Figure(; size=(1800, 1500))
    for (i, τB) in enumerate(tauB_list), (j, T) in enumerate(T_list)
        pars = BSH.BAMParams(; τA=τA, τB=τB, τocc=τocc, γ=γ, movement=movement, T=T)

        rel  = BSH.relative_loss_curves(rng, pool, grid, pars; loss_fracs, seed_A)
        plac = with_placebo ? BSH.placebo_curves(; rng, pool, grid, pars, loss_fracs, seed_A) : nothing

        ax = Axis(fig[i, j])
        _mini_panel!(ax, rel, plac; title=@sprintf("τB=%.2f,  T=%d", τB, T), with_placebo=with_placebo)
    end

    Label(fig[0, 1:3],
        "Archetype = $(String(archetype)) — movement $(movement) — $name grid",
        fontsize=20, tellwidth=false)

    save(outpath, fig)
    return fig
end

# === NEW: Fig.1 + P_fail, 3×3 sweep =========================================
export sweep_fig1_pair_3x3

# compact pair inside one cell: left (AM dashed vs BAM solid [+ placebo opt]),
# right (P_fail among A-suitable kept cells, BAM only)
function _mini_pair!(
    gl::GridLayout, rel::Dict{Symbol,Any},
    pf::Dict{Symbol,Any}, plac::Union{Nothing,Dict{Symbol,Any}};
    with_placebo::Bool, title::String=""
)
    geoms  = (:random, :clustered, :front)
    colors = (to_color(:dodgerblue), to_color(:darkorange), to_color(:seagreen))

    axL = Axis(gl[1,1]); axR = Axis(gl[1,2])

    # left: AM (dash) vs BAM (solid) (+ placebo dotted if requested)
    for (g, col) in zip(geoms, colors)
        r = rel[g]
        lines!(axL, r.x, r.relAM;  color=col, linestyle=:dash,  linewidth=2)
        lines!(axL, r.x, r.relBAM; color=col, linestyle=:solid, linewidth=3)
        if with_placebo && plac !== nothing
            p = plac[g]
            lines!(axL, p.x, p.rel; color=col, linestyle=:dot, linewidth=2)
        end
    end
    axL.ylabel = "relative ΔBSH"
    hidespines!(axL, :t, :r)

    # right: P_fail
    for (g, col) in zip(geoms, colors)
        p = pf[g]
        lines!(axR, p.x, p.y; color=col, linewidth=2.5)
    end
    axR.ylabel = "P_fail (A-suitable kept cells)"
    hidespines!(axR, :t, :r)

    Label(gl[0,1:2], title; fontsize=14, tellwidth=false)
    return nothing
end

"""
sweep_fig1_pair_3x3(; rng, archetype=:mid, grids_vec, grid_type="patchy",
                    tauB_list=[0.40,0.50,0.60], T_list=[6,8,12],
                    S=175, basal_frac=0.35, loss_fracs=0.2:0.1:0.8, seed_A=1,
                    τA=0.5, τocc=0.2, γ=3.0, movement=:component,
                    with_placebo=false,
                    outpath="figs/fig1_pair_sweep_mid.png")
Build a 3×3 figure; each cell shows **(left)** AM vs BAM curves and
**(right)** the corresponding P_fail lines for the same parameter combo.
"""
function sweep_fig1_pair_3x3(;
    rng, archetype::Symbol=:mid,
    grids_vec::Vector{Tuple{String,Grid}}, grid_type::String="patchy",
    tauB_list=[0.40,0.50,0.60], T_list=[6,8,12],
    S::Int=175, basal_frac::Float64=0.35, loss_fracs=0.2:0.1:0.8, seed_A::Int=1,
    τA::Float64=0.5, τocc::Float64=0.2, γ::Float64=3.0,
    movement::Symbol=:component, with_placebo::Bool=false,
    outpath::AbstractString="figs/fig1_pair_sweep_$(String(archetype)).png"
)
    # metaweb for this archetype
    pool = Metawebs.build_metaweb_archetype(rng; S=S, basal_frac=basal_frac, archetype=archetype)

    # pick grid by name from provided vector
    idx = findfirst(x -> x[1] == grid_type, grids_vec)
    idx === nothing && error("grid_type='$grid_type' not found in grids_vec")
    gname, grid = grids_vec[idx]

    # layout: each cell reserves two axes (left/right)
    fig = Figure(; size=(2200, 1700))
    for (i, τB) in enumerate(tauB_list), (j, T) in enumerate(T_list)
        pars = BSH.BAMParams(; τA=τA, τB=τB, τocc=τocc, γ=γ, movement=movement, T=T)

        rel  = BSH.relative_loss_curves(rng, pool, grid, pars; loss_fracs, seed_A)
        pf   = BSH.pfail_curve(; rng, pool, grid, pars, loss_fracs, seed_A)
        plac = with_placebo ? BSH.placebo_curves(; rng, pool, grid, pars, loss_fracs, seed_A) : nothing

        gl = GridLayout(fig[i, j])
        _mini_pair!(gl, rel, pf, plac;
                    with_placebo, title=@sprintf("τB=%.2f,  T=%d", τB, T))
    end

    Label(fig[0, 1:3],
        "Archetype=$(String(archetype)) — movement=$(String(movement)) — $gname grid",
        fontsize=20, tellwidth=false)

    save(outpath, fig)
    return fig
end

"""
sweep_fig1_3x3_for_pool(pool; rng, grid, tauB_list=[0.40,0.50,0.60], T_list=[6,8,12],
                        τA=0.5, τocc=0.2, γ=3.0, movement=:component,
                        loss_fracs=0.2:0.1:0.8, seed_A=1,
                        with_placebo=false, title="",
                        outpath="figs/fig1_sweep_custom.png")

Make the same 3×3 Fig-1 style grid **for a provided SpeciesPool**.
"""
function sweep_fig1_3x3_for_pool(pool::Metawebs.SpeciesPool; rng,
    grid::Grids.Grid,
    tauB_list=[0.40,0.50,0.60], T_list=[6,8,12],
    τA=0.5, τocc=0.2, γ=3.0, movement::Symbol=:component,
    loss_fracs=0.2:0.1:0.8, seed_A::Int=1,
    with_placebo::Bool=false, title::AbstractString="",
    outpath::AbstractString="Paper/figs/fig1_sweep_custom.png")

    fig = Figure(; size=(1800,1500))
    for (i, τB) in enumerate(tauB_list), (j, T) in enumerate(T_list)
        pars = BSH.BAMParams(; τA=τA, τB=τB, τocc=τocc, γ=γ, movement=movement, T=T)
        rel  = BSH.relative_loss_curves(rng, pool, grid, pars; loss_fracs, seed_A)
        plac = with_placebo ? BSH.placebo_curves(; rng, pool, grid, pars, loss_fracs, seed_A) : nothing
        _mini_panel!(Axis(fig[i,j]), rel, plac; title=@sprintf("τB=%.2f,  T=%d", τB, T),
                     with_placebo=with_placebo)
    end
    Label(fig[0, 1:3], isempty(title) ? "Custom metaweb" : title, fontsize=20, tellwidth=false)
    save(outpath, fig)
    fig
end

end # module