# src/metaweb_variety_sweep.jl
module MetawebVarietySweep

using Random, Printf
using CairoMakie
using ..Grids
using ..Metawebs
using ..BSH
using ..FigSweeps

export run_metaweb_variety_sweeps

"""
run_metaweb_variety_sweeps(; rng, grids_vec, grid_type="patchy",
    S=175, tauB_list=[0.40,0.50,0.60], T_list=[6,8,12],
    movement=:component, τA=0.5, τocc=0.2, γ=3.0,
    loss_fracs=0.2:0.1:0.8, seed_A=1, outdir="figs/metaweb_variety")

Build a *wider* set of metawebs (niche / modular / powerlaw + your low/mid/high),
and save a 3×3 mini-panel for each, on the chosen grid_type.
"""
function run_metaweb_variety_sweeps(; rng::AbstractRNG,
    grids_vec::Vector{Tuple{String,Grids.Grid}}, grid_type::String="patchy",
    S::Int=175,
    tauB_list=[0.40,0.50,0.60], T_list=[6,8,12],
    movement::Symbol=:component, τA::Float64=0.5, τocc::Float64=0.2, γ::Float64=3.0,
    loss_fracs=0.2:0.1:0.8, seed_A::Int=1,
    outdir::AbstractString="Paper/figs/metaweb_variety",
    display_figs::Bool=true)

    mkpath(outdir)
    # pick grid
    idx  = findfirst(x->x[1]==grid_type, grids_vec)
    grid = grids_vec[idx][2]

    pools = Vector{Tuple{String,Metawebs.SpeciesPool}}()

    # existing archetypes
    push!(pools, ("arche_low",  Metawebs.build_metaweb_archetype(rng; S=S, basal_frac=0.35, archetype=:low)))
    push!(pools, ("arche_mid",  Metawebs.build_metaweb_archetype(rng; S=S, basal_frac=0.35, archetype=:mid)))
    push!(pools, ("arche_high", Metawebs.build_metaweb_archetype(rng; S=S, basal_frac=0.35, archetype=:high)))

    # niche-like
    for β in (1.2, 1.8, 2.6)
        push!(pools, (@sprintf("niche_beta%.1f", β),
                      Metawebs.build_metaweb_niche(rng; S=S, beta=β)))
    end

    # modular
    push!(pools, ("mod_K3_pin0.30_pout0.04",
                  Metawebs.build_metaweb_modular(rng; S=S, K=3, p_in=0.30, p_out=0.04)))
    push!(pools, ("mod_K4_pin0.25_pout0.02",
                  Metawebs.build_metaweb_modular(rng; S=S, K=4, p_in=0.25, p_out=0.02)))

    # powerlaw
    push!(pools, ("pow_α2.2", Metawebs.build_metaweb_powerlaw(rng; S=S, alpha=2.2, kmax=8)))
    push!(pools, ("pow_α2.8", Metawebs.build_metaweb_powerlaw(rng; S=S, alpha=2.8, kmax=8)))

    # sweep each pool
    for (name, pool) in pools
        fig = FigSweeps.sweep_fig1_3x3_for_pool(pool; rng, grid,
            tauB_list=tauB_list, T_list=T_list,
            τA=τA, τocc=τocc, γ=γ, movement=movement,
            loss_fracs=loss_fracs, seed_A=seed_A, with_placebo=false,
            title="metaweb = $(name) — movement $(movement) — $(grid_type) grid")

        outpath_file = joinpath(outdir, @sprintf("fig1sweep_%s_%s.png", grid_type, name))
        save(outpath_file, fig)
        display_figs && display(fig)
    end

    return nothing
end

end # module
