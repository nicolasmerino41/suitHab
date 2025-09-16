
# --- small helper: consumer-level BSH for one pool/mask
function _per_species_AM_ABM(pool::SpeciesPool, grid::Grid, A::Matrix{Float64},
                             keep::BitVector; M_level::Symbol, τA::Float64, τocc::Float64)
    parsB = bam_from_axes(; B_level=:soft, M_level=M_level, τA=τA, τocc=τocc)  # B param here is ignored below
    # Full (use the caller's B-level via bam passed in separately if you prefer)
    # We'll pass bam explicitly from caller to be precise.

    error("Use the wrapper below; this internal stub is bypassed")
end

# Build distributions across seeds for one geometry at loss f
function collect_distributions_at_loss(; grid::Grid, S::Int, basal_frac::Float64,
        A_level::Symbol, B_level::Symbol, M_level::Symbol,
        hl_kind::Symbol, loss_pick::Float64,
        seeds_pool, seeds_mask, sim_seed::Int,
        nseeds_cluster::Int, front_axis::Symbol, front_noise::Float64,
        τA::Float64, τocc::Float64)

    keepfrac = 1 - loss_pick
    vals_AM  = Float64[]   # concatenated across consumers & seeds
    vals_ABM = Float64[]

    make_keepmask = function(rng::AbstractRNG, kind::Symbol, keep_frac::Float64)
        kind === :random   && return random_mask(rng, grid.C, keep_frac)
        kind === :clustered&& return clustered_mask(rng, grid, keep_frac; nseeds=nseeds_cluster)
        kind === :front    && return frontlike_mask(rng, grid, keep_frac; axis=front_axis, noise=front_noise)
        error("Unknown hl_kind=$kind")
    end

    for ps in seeds_pool, ms in seeds_mask
        rng_pool = MersenneTwister(hash((sim_seed,:pool,ps)))
        pool     = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
        A        = abiotic_matrix(pool, grid)
        pars     = bam_from_axes(; B_level=B_level, M_level=M_level, τA=τA, τocc=τocc)
        bam, mp  = pars.bam, pars.mp

        rng_mask = MersenneTwister(hash((sim_seed,:mask,ms,hl_kind,loss_pick)))
        keep     = make_keepmask(rng_mask, hl_kind, keepfrac)

        # Full
        Pfull, _, _, _ = assemble_BAM(pool, grid, A, keep; bam=bam, mp=mp)
        bsh_full       = bsh1_per_species(Pfull, (A .>= τA), pool)[consumer_mask(pool)]

        # Climate-only (Biotic OFF)
        pars_AM = bam_from_axes(; B_level=:none, M_level=M_level, τA=τA, τocc=τocc)
        P_am, _, _, _ = assemble_BAM(pool, grid, A, keep; bam=pars_AM.bam, mp=pars_AM.mp)
        bsh_AM        = bsh1_per_species(P_am, (A .>= τA), pool)[consumer_mask(pool)]

        append!(vals_AM,  bsh_AM)
        append!(vals_ABM, bsh_full)
    end
    return vals_AM, vals_ABM
end

# Simple boxplot (pure Makie; no external deps)
function simple_boxplot!(ax::Axis, xpos::Float64, data::Vector{Float64};
                         width::Float64=0.28, color=:gray, strokecolor=:black, label=nothing)
    isempty(data) && return nothing
    q1, med, q3 = quantile(data, (0.25, 0.50, 0.75))
    iqr = q3 - q1
    lo  = maximum([minimum(data), q1 - 1.5*iqr])
    hi  = minimum([maximum(data), q3 + 1.5*iqr])
    # box
    poly!(ax, Point2f[(xpos-width/2,q1),(xpos+width/2,q1),(xpos+width/2,q3),(xpos-width/2,q3)],
          color=color, strokecolor=strokecolor)
    # median
    lines!(ax, [xpos-width/2, xpos+width/2], [med, med]; color=strokecolor, linewidth=2)
    # whiskers
    lines!(ax, [xpos, xpos], [q3, hi]; color=strokecolor)
    lines!(ax, [xpos, xpos], [q1, lo]; color=strokecolor)
    # caps
    lines!(ax, [xpos-0.12, xpos+0.12], [hi, hi]; color=strokecolor)
    lines!(ax, [xpos-0.12, xpos+0.12], [lo, lo]; color=strokecolor)
    return nothing
end

"""
plot_bsh_distributions_at_loss(; grid, S, basal_frac, A_level, B_level, M_level,
                               loss_pick, seeds_pool, seeds_mask, sim_seed,
                               nseeds_cluster, front_axis, front_noise, τA, τocc)

Produces a 1×1 figure:
  X axis = geometry (Random/Clustered/Front). For each geometry, two side-by-side boxes:
  left = A×M (climate-only), right = A×B×M (full BSH). Values are per-consumer.
"""
function plot_bsh_distributions_at_loss(; grid::Grid, S::Int=150, basal_frac::Float64=0.45,
        A_level::Symbol=:intermediate, B_level::Symbol=:soft, M_level::Symbol=:on,
        loss_pick::Float64=0.6, seeds_pool=1:6, seeds_mask=1:6, sim_seed::Int=1234,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        τA::Float64=0.5, τocc::Float64=0.35)

    R_AM, R_ABM = collect_distributions_at_loss(; grid,S,basal_frac,A_level,B_level,M_level,
                                                hl_kind=:random, loss_pick,
                                                seeds_pool,seeds_mask,sim_seed,
                                                nseeds_cluster,front_axis,front_noise,τA,τocc)
    C_AM, C_ABM = collect_distributions_at_loss(; grid,S,basal_frac,A_level,B_level,M_level,
                                                hl_kind=:clustered, loss_pick,
                                                seeds_pool,seeds_mask,sim_seed,
                                                nseeds_cluster,front_axis,front_noise,τA,τocc)
    F_AM, F_ABM = collect_distributions_at_loss(; grid,S,basal_frac,A_level,B_level,M_level,
                                                hl_kind=:front, loss_pick,
                                                seeds_pool,seeds_mask,sim_seed,
                                                nseeds_cluster,front_axis,front_noise,τA,τocc)

    fig = Figure(; size=(1200, 500))
    Label(fig[0,:],
          @sprintf("Per-species BSH at loss=%.2f — A=%s, B=%s, M=%s", loss_pick,
                   String(A_level), String(B_level), String(M_level));
          fontsize=16, padding=(0,0,8,0))
    ax = Axis(fig[1,1], ylabel="BSH (fraction)", xticks=(1:3, ["Random","Clustered","Front-like"]))

    # positions (pair per geometry)
    offs = 0.18
    # Random
    simple_boxplot!(ax, 1-offs, R_AM;  color=RGBAf(0.30,0.45,0.85,0.30), strokecolor=:steelblue)
    simple_boxplot!(ax, 1+offs, R_ABM; color=RGBAf(0.20,0.65,0.35,0.30), strokecolor=:seagreen)
    # Clustered
    simple_boxplot!(ax, 2-offs, C_AM;  color=RGBAf(0.30,0.45,0.85,0.30), strokecolor=:steelblue)
    simple_boxplot!(ax, 2+offs, C_ABM; color=RGBAf(0.20,0.65,0.35,0.30), strokecolor=:seagreen)
    # Front-like
    simple_boxplot!(ax, 3-offs, F_AM;  color=RGBAf(0.30,0.45,0.85,0.30), strokecolor=:steelblue)
    simple_boxplot!(ax, 3+offs, F_ABM; color=RGBAf(0.20,0.65,0.35,0.30), strokecolor=:seagreen)

    # legend proxies
    l1 = lines!(ax, [NaN,NaN], [NaN,NaN]; color=:steelblue)
    l2 = lines!(ax, [NaN,NaN], [NaN,NaN]; color=:seagreen)
    axislegend(ax, [l1,l2], ["Climate-only (A×M)", "Full BSH (A×B×M)"];
               position=:rt, framevisible=false)

    display(fig)
    return fig
end

# (2) Distributions at a chosen loss
_ = plot_bsh_distributions_at_loss(; grid,
       S=150, basal_frac=0.45,
       A_level=:divergent, B_level=:strong, M_level=:on,
       loss_pick=0.0,
       seeds_pool=1:6, seeds_mask=1:6, sim_seed=1234,
       nseeds_cluster=6, front_axis=:x, front_noise=0.04,
       τA=0.5, τocc=0.35)