"""
Run all computations for the BAM dashboard and return results as a NamedTuple.
"""
function run_dashboard_analysis(; nx::Int=60, ny::Int=60, S::Int=150, basal_frac::Float64=0.45,
        A_level::Symbol=:intermediate, B_level::Symbol=:soft, M_level::Symbol=:on,
        loss_fracs=0.2:0.1:0.8, show_loss_for_fingerprints::Float64=0.6,
        seeds_pool=1:6, seeds_mask=1:6, sim_seed::Int=1234,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        τA::Float64=0.5, τocc::Float64=0.35, T_frac_on::Float64=0.02)

    grid = make_grid(nx, ny; seed=42)

    # ΔBSH curves
    dR = sweep_dBSH_axes(; grid,S,basal_frac,A_level,B_level,M_level,hl_kind=:random,
                         loss_fracs,seeds_pool,seeds_mask,sim_seed,nseeds_cluster,
                         front_axis,front_noise,τA,τocc,T_frac_on)
    dC = sweep_dBSH_axes(; grid,S,basal_frac,A_level,B_level,M_level,hl_kind=:clustered,
                         loss_fracs,seeds_pool,seeds_mask,sim_seed,nseeds_cluster,
                         front_axis,front_noise,τA,τocc,T_frac_on)
    dF = sweep_dBSH_axes(; grid,S,basal_frac,A_level,B_level,M_level,hl_kind=:front,
                         loss_fracs,seeds_pool,seeds_mask,sim_seed,nseeds_cluster,
                         front_axis,front_noise,τA,τocc,T_frac_on)

    # fingerprints (separate figure, not embedded)
    raw  = summarize_all_combos(; nx, ny, S, basal_frac,
                                          loss_pick=show_loss_for_fingerprints,
                                          seeds_pool, seeds_mask, sim_seed)
    elas = elasticity_summarize_all(; nx, ny, S, basal_frac,
                                              loss_pick=show_loss_for_fingerprints,
                                              seeds_pool, seeds_mask, sim_seed)

    # fragmentation diagnostics
    fR = frag_diagnostic_axes(; grid,S,basal_frac,A_level,B_level,M_level,hl_kind=:random,
                               loss_fracs,seeds_pool,seeds_mask,sim_seed,nseeds_cluster,
                               front_axis,front_noise,τA,τocc,T_frac_on)
    fC = frag_diagnostic_axes(; grid,S,basal_frac,A_level,B_level,M_level,hl_kind=:clustered,
                               loss_fracs,seeds_pool,seeds_mask,sim_seed,nseeds_cluster,
                               front_axis,front_noise,τA,τocc,T_frac_on)
    fF = frag_diagnostic_axes(; grid,S,basal_frac,A_level,B_level,M_level,hl_kind=:front,
                               loss_fracs,seeds_pool,seeds_mask,sim_seed,nseeds_cluster,
                               front_axis,front_noise,τA,τocc,T_frac_on)

    # baseline + prediction
    base = baseline_indicators(grid, S, basal_frac;
                               A_level,B_level,M_level,τA,τocc,T_frac_on)
    pred = predict_geometry_ranking(base.D, base.R, base.Θ;
                                    movement_on=(M_level===:on))

    return (; grid, A_level, B_level, M_level,
            dR, dC, dF, fR, fC, fF,
            raw, elas, base, pred,
            loss_fracs, show_loss_for_fingerprints)
end

"""
Plot BAM dashboard panels given results from `run_dashboard_analysis`.
"""
function plot_dashboard(res)
    # color palette (stable across panels)
    colors = Dict(
        :random    => RGBf(0.25, 0.45, 0.85),
        :clustered => RGBf(0.95, 0.60, 0.15),
        :front     => RGBf(0.20, 0.65, 0.35),
    )

    fig = Figure(; size = (1400, 900))
    Label(fig[0, :],
          "BAM Dashboard — A=$(res.A_level), B=$(res.B_level), M=$(res.M_level)";
          fontsize = 18, padding = (0, 0, 8, 0))

    # ── Top-left: ΔBSH vs loss ──────────────────────────────────────────────
    ax_d = Axis(fig[1, 1];
        title = "ΔBSH vs loss",
        xlabel = "Area lost (fraction)",
        ylabel = "ΔBSH (mean consumers)")

    hR = lines!(ax_d, res.dR.x, res.dR.y; color = colors[:random],    linewidth = 2)
    hC = lines!(ax_d, res.dC.x, res.dC.y; color = colors[:clustered], linewidth = 2)
    hF = lines!(ax_d, res.dF.x, res.dF.y; color = colors[:front],     linewidth = 2)

    axislegend(ax_d, [hR, hC, hF], ["Random", "Clustered", "Front-like"];
               position = :lt, framevisible = false, labelsize = 10, padding = (2,2,2,2))

    # ── Top-right: Fragmentation diagnostic ─────────────────────────────────
    ax_f = Axis(fig[1, 2];
        title = "Fragmentation diagnostic",
        xlabel = "Area lost (fraction)",
        ylabel = "Share of suitable cells with M = 1")

    gR = lines!(ax_f, res.fR.x, res.fR.y; color = colors[:random],    linewidth = 2)
    gC = lines!(ax_f, res.fC.x, res.fC.y; color = colors[:clustered], linewidth = 2)
    gF = lines!(ax_f, res.fF.x, res.fF.y; color = colors[:front],     linewidth = 2)
    hlines!(ax_f, [0.5]; color = RGBAf(0,0,0,0.35), linestyle = :dash)

    axislegend(ax_f, [gR, gC, gF], ["Random", "Clustered", "Front-like"];
               position = :rb, framevisible = false, labelsize = 10, padding = (2,2,2,2))

    # ── Bottom: Baseline indices & prediction (tidy label stack) ────────────
    info = fig[2, 1:2] = GridLayout()
    Label(info[1,1], "Baseline indices & prediction"; halign = :left, fontsize = 14,
          padding = (0,0,6,0))

    gl = info[2,1] = GridLayout()
    Label(gl[1,1], @sprintf("D (climate divergence): %.2f", res.base.D); halign = :left)
    Label(gl[2,1], @sprintf("R (diet redundancy): mean=%.2f, p95=%.2f",
                             res.base.R.mean, res.base.R.p95); halign = :left)
    Label(gl[3,1], @sprintf("Θ (movement stiffness): %.2f   (T=%d)",
                             res.base.Θ, res.base.T); halign = :left)
    Label(gl[4,1], @sprintf("LPRI (tail prey): mean=%.2f, min=%.2f",
                             res.base.LPRI.mean, res.base.LPRI.min); halign = :left)
    Label(gl[5,1], "Predicted ranking (worst → best):"; halign = :left,
          padding = (6,0,0,0))

    pr = res.pred
    rank_str = string(pr.worst, "  >  ", get(pr, :middle, :clustered), "  >  ", pr.best)
    Label(gl[6,1], rank_str; halign = :left)
    if haskey(pr, :rationale)
        Label(gl[7,1], "Rationale: "*pr.rationale; halign = :left,
              fontsize = 10, color = RGBAf(0,0,0,0.7))
    end

    display(fig)
    # return fig
end

res = run_dashboard_analysis(
    ; A_level=:divergent, B_level=:strong, M_level=:on,
    loss_fracs=0.2:0.1:0.8, show_loss_for_fingerprints=0.6,
    T_frac_on=0.03
)
plot_dashboard(res)