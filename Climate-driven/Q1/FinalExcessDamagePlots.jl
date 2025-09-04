# === excess curves vs random ===============================================
# Assumes you already defined:
#   - decomp_at_mask_fraction_safe(pool, grid; τ, keepmask)
#   - random_mask, clustered_mask, consumer_hotspot_mask, prey_hotspot_mask

# One decomposition curve across losses for a given scenario
function decomp_curve(pool::SpeciesPool, grid::Grid; τ::Float64,
                      keeps::Vector{Float64}, kind::Symbol,
                      nseeds_cluster::Int=1, hotspot_power::Float64=2.5,
                      seed::Int=101)
    clim = Float64[]; inter = Float64[]; syn = Float64[]
    for kf in keeps
        km = if kind === :random
            random_mask(grid.C, kf; seed=seed)
        elseif kind === :clustered
            clustered_mask(grid, kf; nseeds=nseeds_cluster, seed=seed)
        elseif kind === :hotspot_cons
            consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=kf, power=hotspot_power, seed=seed)
        elseif kind === :hotspot_prey
            prey_hotspot_mask(grid, pool; τ=τ, keep_frac=kf, power=hotspot_power, seed=seed)
        elseif kind === :hotspot_prey_resid
            prey_hotspot_residual_mask(grid, pool; τ=τ, keep_frac=kf, power=hotspot_power, seed=seed)
        elseif kind === :hotspot_cons_resid
            consumer_hotspot_residual_mask(grid, pool; τ=τ, keep_frac=kf, power=hotspot_power, seed=seed)
        else
            error("Unknown kind = $kind")
        end
        d = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=km)
        push!(clim, d.clim); push!(inter, d.inter); push!(syn, d.syn)
    end
    return (; loss = 1 .- keeps, clim, inter, syn)
end

# Excess = scenario – random (component-wise)
function excess_curves(pool::SpeciesPool, grid::Grid; τ::Float64,
                       keeps::Vector{Float64}, hotspot_power::Float64=2.5,
                       nseeds_cluster::Int=1, seed::Int=101)
    base = decomp_curve(pool, grid; τ=τ, keeps=keeps, kind=:random, seed=seed)
    make_diff(kind) = let C = decomp_curve(pool, grid; τ=τ, keeps=keeps,
                                           kind=kind, hotspot_power=hotspot_power,
                                           nseeds_cluster=nseeds_cluster, seed=seed)
        (; loss = base.loss,
           clim  = C.clim  .- base.clim,
           inter = C.inter .- base.inter,
           syn   = C.syn   .- base.syn)
    end
    clustered     = make_diff(:clustered)
    hotspot_cons  = make_diff(:hotspot_cons)
    hotspot_prey  = make_diff(:hotspot_prey)
    hotspot_prey_resid  = make_diff(:hotspot_prey_resid)
    hotspot_cons_resid  = make_diff(:hotspot_cons_resid)
    return (; clustered, hotspot_cons, hotspot_prey, hotspot_prey_resid, hotspot_cons_resid)
end

# ---- run & plot ------------------------------------------------------------
begin
    grid = make_grid(60,60; seed=11)
    pool = build_pool(200;
        basal_frac=0.35, seed=1,
        sigma=0.25, density=0.18, pmax=0.75,
        niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.04,
        b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04
    )
    τ      = 0.55
    keeps  = collect(0.2:0.05:0.85)     # loss = 1 - keep
    pow    = 2.5

    E = excess_curves(pool, grid; τ=τ, keeps=keeps, hotspot_power=pow)

    function panel!(ax, C; title="")
        lines!(ax, C.loss, C.clim,  color=:dodgerblue, label="Climate-only")
        lines!(ax, C.loss, C.inter, color=:orange,     label="Interaction-only")
        lines!(ax, C.loss, C.syn,   color=:forestgreen, label="Synergy")
        hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
        ax.title = title
        ax.xlabel = "Area lost (fraction)"
        ax.ylabel = "Excess ΔBSH (scenario − random)"
    end

    fig = Figure(; size=(1080,360))
    ax1 = Axis(fig[1,1])
    panel!(ax1, E.clustered; title="Clustered")
    ax2 = Axis(fig[1,2])
    panel!(ax2, E.hotspot_cons; title="Hotspot — Consumers")
    ax3 = Axis(fig[1,3])
    panel!(ax3, E.hotspot_prey; title="Hotspot — Prey")
    ax4 = Axis(fig[1,4])
    panel!(ax4, E.hotspot_cons_resid; title="Hotspot — Consumers (residual)")
    ax5 = Axis(fig[1,5])
    panel!(ax5, E.hotspot_prey_resid; title="Hotspot — Prey (residual)")

    axislegend(ax1; position=:lb)
    # panel!(Axis(fig[1,1]), E.clustered;    title="Clustered")
    # panel!(Axis(fig[1,2]), E.hotspot_cons; title="Hotspot — Consumers")
    # panel!(Axis(fig[1,3]), E.hotspot_prey; title="Hotspot — Prey")
    # panel!(Axis(fig[1,4]), E.hotspot_prey_resid; title="Hotspot — Prey (residual)")
    # axislegend(fig[1,1], position=:lt)
    display(fig)
end
