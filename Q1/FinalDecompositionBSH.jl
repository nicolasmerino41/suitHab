#
# === CELL SCORES WE’LL USE FOR HOTSPOTS ======================================
# - consumer_score[c] = fraction of consumers present in cell c
# - prey_score[c]      = fraction of consumers that have ≥1 prey present in c
#
function _consumer_and_prey_scores(pool::SpeciesPool, grid::Grid; τ::Float64)
    Z = climate_pass(pool, grid; τ=τ)
    P = assemble(Z, pool)                  # S×C presence matrix (Bool)
    S, C = size(P)
    cons = .!pool.basal

    # consumer richness (fraction of all consumers present)
    ncons = count(cons)
    cons_count = vec(sum(@view P[cons, :]; dims=1))
    cons_score = ncons > 0 ? cons_count ./ ncons : zeros(C)

    # for each consumer s, do they have ≥1 prey present in cell c?
    prey_support = zeros(Float64, C)
    for c in 1:C
        good = 0
        for s in 1:S
            cons[s] || continue
            prey = pool.E[s]
            if !isempty(prey)
                # any prey present in this cell?
                hit = false
                @inbounds for q in prey
                    if P[q, c]; hit = true; break; end
                end
                good += hit ? 1 : 0
            end
        end
        prey_support[c] = ncons > 0 ? good / ncons : 0.0
    end
    return (cons_score, prey_support)  # both length C in [0,1]
end

# Utility: select the top keep_frac of cells by a score (with optional “power” skew)
function _top_mask_by_score(score::AbstractVector{<:Real}, keep_frac::Float64; power::Float64=2.5, seed::Int=0)
    C = length(score)
    nkeep = clamp(round(Int, keep_frac * C), 1, C)
    # skew and tiny jitter to break ties
    w = (score .^ power) .+ 1e-12 .* rand(MersenneTwister(seed), C)
    idx = partialsortperm(w, rev=true, 1:nkeep)
    keepmask = falses(C); keepmask[idx] .= true
    return keepmask
end

"Hotspot of consumer presence (cells rich in consumers). Returns keepmask."
function consumer_hotspot_mask(grid::Grid, pool::SpeciesPool; τ::Float64, keep_frac::Float64, power::Float64=2.5, seed::Int=0)
    cons_score, _ = _consumer_and_prey_scores(pool, grid; τ=τ)
    return _top_mask_by_score(cons_score, keep_frac; power=power, seed=seed)
end

"Hotspot of prey support (cells where many consumers find ≥1 prey). Returns keepmask."
function prey_hotspot_mask(grid::Grid, pool::SpeciesPool; τ::Float64, keep_frac::Float64, power::Float64=2.5, seed::Int=0)
    _, prey_score = _consumer_and_prey_scores(pool, grid; τ=τ)
    return _top_mask_by_score(prey_score, keep_frac; power=power, seed=seed)
end

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

# --- example run + plot
begin
    grid = make_grid(60,60; seed=11)
    S = 200; basal_frac = 0.35
    # pick any regime you like here:
    pool = build_pool(S;
        basal_frac=basal_frac, seed=1,
        sigma=0.25, density=0.18, pmax=0.75,
        niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.04,
        b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04
    )
    τ = 0.55
    keeps = collect(0.2:0.05:0.85)      # area kept; loss = 1 - keep
    pow  = 2.5

    C_rand  = decomp_curve(pool, grid; τ=τ, keeps=keeps, kind=:random)
    C_clust = decomp_curve(pool, grid; τ=τ, keeps=keeps, kind=:clustered, nseeds_cluster=1)
    C_hc    = decomp_curve(pool, grid; τ=τ, keeps=keeps, kind=:hotspot_cons, hotspot_power=pow)
    C_hp    = decomp_curve(pool, grid; τ=τ, keeps=keeps, kind=:hotspot_prey, hotspot_power=pow)
    C_hc_r = decomp_curve(pool, grid; τ=τ, keeps=keeps, kind=:hotspot_cons_resid, hotspot_power=pow, seed=101)
    C_hp_r = decomp_curve(pool, grid; τ=τ, keeps=keeps, kind=:hotspot_prey_resid, hotspot_power=pow, seed=101)

    function plot_panel!(ax, C; title="")
        lines!(ax, C.loss, C.clim,  color=:dodgerblue, label="Climate-only")
        lines!(ax, C.loss, C.inter, color=:orange,     label="Interaction-only")
        lines!(ax, C.loss, C.syn,   color=:forestgreen,label="Synergy")
        hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
        ax.title = title
        ax.xlabel = "Area lost (fraction)"
        ax.ylabel = "ΔBSH decomposition"
    end

    fig = Figure(; size=(1120,600))
    plot_panel!(Axis(fig[1,1]), C_rand;  title="Random")
    plot_panel!(Axis(fig[1,2]), C_clust; title="Clustered")
    plot_panel!(Axis(fig[1,3]), C_hc;    title="Hotspot — Consumers")
    plot_panel!(Axis(fig[2,1]), C_hp;    title="Hotspot — Prey")
    plot_panel!(Axis(fig[2,2]), C_hc_r;  title="Hotspot — Consumers Resid")
    plot_panel!(Axis(fig[2,3]), C_hp_r;  title="Hotspot — Prey Resid")
    # axislegend(fig[1,1], position=:lt)
    display(fig)
end

