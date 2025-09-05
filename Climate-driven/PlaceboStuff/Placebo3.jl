# Map any of the possible field names to a single canonical name
_get(nt::NamedTuple, alts::Tuple) = begin
    for k in alts
        if hasproperty(nt, k); return getfield(nt, k); end
    end
    error("Expected one of $(alts) in $(propertynames(nt))")
end

"""
collect_decomp_series(vec_of_namedtuples)

Accepts a Vector of NamedTuples that can use either:
  (:clim, :inter, :syn, :tot)  OR  (:dA_only, :dInt_only, :synergy, :dB)

Returns (clim, inter, syn, tot) as Float64 vectors.
"""
function collect_decomp_series(vec::Vector{<:NamedTuple})
    clim   = Float64[_get(nt, (:clim,   :dA_only))   for nt in vec]
    inter  = Float64[_get(nt, (:inter,  :dInt_only)) for nt in vec]
    syn    = Float64[_get(nt, (:syn,    :synergy))   for nt in vec]
    tot    = Float64[_get(nt, (:tot,    :dB))        for nt in vec]  # total ΔBSH
    return (; clim, inter, syn, tot)
end

"""
excess_from_decomps(R, S)

Given two series (e.g., R = random, S = scenario) returned by
collect_decomp_series, compute excess components and summaries.
"""
function excess_from_decomps(R, S)
    ED_clim = S.clim  .- R.clim
    ED_int  = S.inter .- R.inter
    ED_syn  = S.syn   .- R.syn
    # total excess either from parts or from totals (both should match up to fp noise)
    ED_tot  = ED_clim .+ ED_int .+ ED_syn
    RED_tot = ED_tot ./ (abs.(R.tot) .+ 1e-12)                # relative excess (unit-free)
    share   = [abs(ED_tot[i]) < 1e-12 ? NaN : ED_int[i]/ED_tot[i] for i in eachindex(ED_tot)]
    return (; ED_clim, ED_int, ED_syn, ED_tot, RED_tot, share)
end

"Left: excess components; Right: interaction share."
function plot_excess_and_share(losses, ex; title="")
    fig = Figure(; size=(950,380))

    ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)",
               ylabel="Excess ΔBSH (scenario − random)", title=title)
    lines!(ax1, losses, ex.ED_clim, label="Climate-only")
    lines!(ax1, losses, ex.ED_int,  label="Interaction-only")
    lines!(ax1, losses, ex.ED_syn,  label="Synergy")
    hlines!(ax1, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax1, position=:lt)

    ax2 = Axis(fig[1,2], xlabel="Area lost (fraction)",
               ylabel="Interaction share of excess")
    lines!(ax2, losses, ex.share)
    hlines!(ax2, [0.0, 1.0], color=(:gray,0.4), linestyle=:dash)

    display(fig); fig
end

"Left: total excess; Right: relative excess (unit-free)."
function plot_total_and_relative(losses, ex; title="")
    fig = Figure(resolution=(900,360))
    ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Excess ΔBSH (total)", title=title)
    lines!(ax1, losses, ex.ED_tot); hlines!(ax1, [0.0], color=(:gray,0.5), linestyle=:dash)

    ax2 = Axis(fig[1,2], xlabel="Area lost (fraction)", ylabel="Relative excess (vs random)")
    lines!(ax2, losses, ex.RED_tot); hlines!(ax2, [0.0], color=(:gray,0.5), linestyle=:dash)

    display(fig); fig
end

"Plot excess for scenario and its rewired placebo (uses *_rw series in results.data)."
function plot_placebo(results; scen::Symbol=:hotspot, title="")
    losses = results.losses
    R  = collect_decomp_series(results.data[:random])
    S  = collect_decomp_series(results.data[scen])
    exS = excess_from_decomps(R, S).ED_tot

    Rw = collect_decomp_series(results.data[:random_rw])
    Sw = collect_decomp_series(results.data[Symbol(scen, "_rw")])
    exP = excess_from_decomps(Rw, Sw).ED_tot

    fig = Figure(resolution=(760,320))
    ax  = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Excess ΔBSH vs random", title=title)
    lines!(ax, losses, exS, label="with interactions")
    lines!(ax, losses, exP, linestyle=:dash, label="rewired placebo")
    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax, position=:lt)
    display(fig); fig
end

begin
    grid    = make_grid(60,60; seed=11)
    S       = 220
    τ       = 0.58
    losses  = collect(0.15:0.05:0.80)
    regimes = RS_regimes_for_max()

    best_res = nothing
    best_name = ""
    best_share = -Inf

    for (name, kw) in regimes
        pool = build_pool(S;
            basal_frac=0.35, seed=1,
            sigma=kw[:sigma], density=kw[:density], pmax=kw[:pmax],
            niche_mode=kw[:niche_mode], mu_basal_centers=(0.25,0.75),
            mu_basal_sd=kw[:mu_basal_sd],
            b0_basal=0.08, bspread_basal=0.02,
            b0_cons=0.12,  bspread_cons=0.04)

        res = run_excess_with_controls(pool, grid; τ=τ, losses=losses,
                                       hotspot=:A_matched, seed=101)

        # interaction share of excess for hotspot at ~50% loss
        dR  = collect_decomp_series(res.data[:random])
        dH  = collect_decomp_series(res.data[:hotspot])
        exH = excess_from_decomps(dR, dH)
        i50 = argmin(abs.(losses .- 0.50))
        share50 = exH.share[i50]

        if share50 > best_share
            best_share = share50
            best_name  = String(name)
            best_res   = res
        end
    end

    println("Best regime @50% loss (interaction share of excess) → ",
            best_name, " ; share = ", round(best_share, digits=3))

    # ---- HOTSPOT (best regime) ----
    dR  = collect_decomp_series(best_res.data[:random])
    dH  = collect_decomp_series(best_res.data[:hotspot])
    exH = excess_from_decomps(dR, dH)

    plot_excess_and_share(best_res.losses, exH; title="A-matched hotspot — $(best_name)")
    plot_total_and_relative(best_res.losses, exH;  title="A-matched hotspot — $(best_name)")
    plot_placebo(best_res; scen=:hotspot, title="Placebo (rewired) — $(best_name)")

    # ---- CLUSTERED (same regime) ----
    dC  = collect_decomp_series(best_res.data[:cluster])
    exC = excess_from_decomps(dR, dC)

    plot_excess_and_share(best_res.losses, exC; title="Clustered — $(best_name)")
    plot_total_and_relative(best_res.losses, exC;  title="Clustered — $(best_name)")
    plot_placebo(best_res; scen=:cluster, title="Placebo (clustered) — $(best_name)")
end