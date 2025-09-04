# --- helpers you already defined elsewhere (assumed available) -------------
# make_grid, build_pool, climate_pass, assemble
# random_mask, clustered_mask, consumer_hotspot_mask
# decomp_at_mask_fraction_safe (fraction-basis, NaN-safe)

# --- choose mask for a scenario at a given loss --------------------------------
"""
choose_mask(scenario; grid, pool, τ, keep, nseeds_cluster, hotspot_power, mode)

scenario ∈ (:clustered, :hotspot)
mode     ∈ (:remove, :protect)  # for hotspot only (default :remove)
"""
function choose_mask(scenario::Symbol;
    grid::Grid, pool::SpeciesPool, τ::Float64, keep::Float64,
    nseeds_cluster::Int=1, hotspot_power::Float64=2.5, mode::Symbol=:remove, seed::Int=101)

    if scenario === :clustered
        return clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=seed)
    elseif scenario === :hotspot
        # build a KEEP mask that keeps top-share cells
        keep_hot = consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep, power=hotspot_power, what=:cons, seed=seed)
        if mode === :protect
            return keep_hot                    # protect-top: keep hotspot cells
        elseif mode === :remove
            return .!keep_hot                  # remove-top: keep everything BUT hotspot cells
        else
            error("mode must be :protect or :remove")
        end
    else
        error("unknown scenario $scenario")
    end
end

# --- Control-1 curves for one (pool, τ) and one scenario -----------------------
"""
control1_curves(pool, grid; τ, losses, scenario, mode=:remove)

Returns:
(losses, ex_with[], ex_noInt[])
where ex_* is (scenario − random) at each loss fraction.
"""
function control1_curves(pool::SpeciesPool, grid::Grid;
    τ::Float64, losses=0.05:0.05:0.8, scenario::Symbol=:clustered,
    nseeds_cluster::Int=1, hotspot_power::Float64=2.5, mode::Symbol=:remove)

    ex_with  = Float64[]
    ex_noInt = Float64[]
    for f in losses
        keep = 1.0 - f
        kmR = random_mask(grid.C, keep; seed=111)
        kmS = choose_mask(scenario; grid, pool, τ, keep, nseeds_cluster, hotspot_power, mode, seed=112)

        dR = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmR)
        dS = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmS)

        # total ΔB with interactions vs random
        Δ_with  = (dS.clim + dS.inter + dS.syn) - (dR.clim + dR.inter + dR.syn)
        # climate-only control (no interactions)
        Δ_noInt = dS.clim - dR.clim

        push!(ex_with,  Δ_with)
        push!(ex_noInt, Δ_noInt)
    end
    return (; losses=collect(losses), ex_with, ex_noInt)
end

# --- parameter grid to explore -------------------------------------------------
function param_sweep_control1(; 
    grid = make_grid(60,60; seed=11),
    S=200, basal_frac=0.35,
    τ_list = [0.50, 0.55, 0.60],
    sigma_list   = [0.22, 0.35, 0.45],
    density_list = [0.12, 0.25, 0.35],
    pmax_list    = [0.70, 0.80, 0.90],
    niche_modes  = [:bimodal, :uniform],
    mu_basal_sd  = 0.04,
    b0_basal=0.08, bspread_basal=0.02,
    b0_cons=0.12,  bspread_cons=0.04,
    hotspot_power=2.5,
    losses = 0.05:0.05:0.8,
    scenarios = [:clustered, :hotspot],
    hotspot_mode = :remove,   # :remove (targeted loss) or :protect
)
    rows = DataFrame()
    fig = Figure(; size=(980,420))
    axs = (Axis(fig[1,1], title="Clustered (with interactions)"),
           Axis(fig[1,2], title="Hotspot (with interactions)"))
    for ax in axs
        ax.xlabel = "Area lost (fraction)"
        ax.ylabel = "Excess (scenario − random)"
        hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dot)
    end

    # spaghetti storage for median (optional)
    store = Dict(:clustered => Float64[], :hotspot => Float64[])
    L = length(losses)

    for τ in τ_list, σ in sigma_list, dens in density_list, pmax in pmax_list, nm in niche_modes
        pool = build_pool(S;
            basal_frac=basal_frac, seed=1,
            sigma=σ, density=dens, pmax=pmax,
            niche_mode=nm, mu_basal_centers=(0.25,0.75), mu_basal_sd=mu_basal_sd,
            b0_basal=b0_basal, bspread_basal=bspread_basal,
            b0_cons=b0_cons,   bspread_cons=bspread_cons)

        for scen in scenarios
            cur = control1_curves(pool, grid; τ=τ, losses=losses,
                                  scenario=scen, mode=(scen==:hotspot ? hotspot_mode : :remove),
                                  hotspot_power=hotspot_power)
            # thin spaghetti line (with interactions)
            ax = scen==:clustered ? axs[1] : axs[2]
            lines!(ax, cur.losses, cur.ex_with; linewidth=1, color=(scen==:clustered ? :steelblue : :darkorange), transparency=true)

            # summary at 50% loss (or nearest)
            idx = findmin(abs.(cur.losses .- 0.5))[2]
            exW = cur.ex_with[idx]
            exN = cur.ex_noInt[idx]

            # qualitative labels
            signlab(x; eps=1e-3) = x > eps ? "pos" : x < -eps ? "neg" : "≈0"
            wlab = signlab(exW); nlab = signlab(exN)
            driver = wlab != nlab ? "interaction-driven" :
                     (wlab=="≈0" ? "neutral" : "climate-aligned")

            push!(rows, (
                τ=τ, sigma=σ, density=dens, pmax=pmax, niche=String(nm),
                scenario=String(scen),
                loss_at_eval=cur.losses[idx],
                excess_with=exW, excess_noInt=exN,
                with_label=wlab, noInt_label=nlab, driver=driver
            ))
        end
    end

    # decorate & show
    axislegend(axs[1], [LineElement(color=:steelblue, width=2)], ["with interactions"]; position=:lt)
    axislegend(axs[2], [LineElement(color=:darkorange, width=2)], ["with interactions"]; position=:lt)
    display(fig)

    df = DataFrame(rows)
    CSV.write("control1_sweep.csv", df)
    @info "Saved sweep to control1_sweep.csv"
    return df
end

# ---- run it ------------------------------------------------------------------
df = param_sweep_control1();

# quick tally table
combine(groupby(df, [:scenario, :with_label, :driver]), nrow => :count)

using AlgebraOfGraphics
plt = data(df) * mapping(:scenario, :with_label, color=:driver) * visual(Scatter; markersize=20)
draw(plt; figure=(;; size=(600,320)))
