# compute Δ’s across geometries for one (A,B)
function _delta_curves_at_f_for_all_geoms(; grid::Grid, S::Int, basal_frac::Float64,
        A_level::Symbol, B_level::Symbol, M_level::Symbol,
        loss_pick::Float64, seeds_pool, seeds_mask, sim_seed::Int,
        nseeds_cluster::Int, front_axis::Symbol, front_noise::Float64,
        τA::Float64, τocc::Float64)

    geoms = (:random,:clustered,:front)
    Δnon = Dict{Symbol,Float64}(); Δbio = Dict{Symbol,Float64}()

    for hk in geoms
        # get AM/ABM at f and at 0
        c  = _sweep_AM_ABM_onegeom(; grid,S,basal_frac,A_level,B_level,M_level,
                                   hl_kind=hk, loss_fracs=[loss_pick],
                                   seeds_pool,seeds_mask,sim_seed,
                                   nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on=0.0)
        c0 = _sweep_AM_ABM_onegeom(; grid,S,basal_frac,A_level,B_level,M_level,
                                   hl_kind=hk, loss_fracs=[0.0],
                                   seeds_pool,seeds_mask,sim_seed,
                                   nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on=0.0)
        Δnon[hk] = c.AM[1]  - c0.AM[1]
        Δbio[hk] = c.ABM[1] - c0.ABM[1]
    end
    return Δnon, Δbio
end

"""
spread_metrics_panel(; grid, S, basal_frac, M_level, loss_pick, A_levels, B_levels, ...)

Scatter: x = Δspread = spread_bio - spread_non ; y = Dmax = max_g Δ_bio - max_g Δ_non
Points labeled by (A,B). Table printed to REPL and a small text block added.
"""
function spread_metrics_panel(; grid::Grid, S::Int=150, basal_frac::Float64=0.45,
        M_level::Symbol=:on, loss_pick::Float64=0.6,
        A_levels=(:neutral,:intermediate,:divergent), B_levels=(:none,:soft,:strong),
        seeds_pool=1:6, seeds_mask=1:6, sim_seed::Int=1234,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        τA::Float64=0.5, τocc::Float64=0.35)

    xs = Float64[]  # Δspread
    ys = Float64[]  # Dmax
    labs = String[]

    for Alev in A_levels, Blev in B_levels
        Δnon, Δbio = _delta_curves_at_f_for_all_geoms(; grid,S,basal_frac,
                          A_level=Alev, B_level=Blev, M_level,
                          loss_pick, seeds_pool,seeds_mask,sim_seed,
                          nseeds_cluster,front_axis,front_noise,τA,τocc)

        non_vals = collect(values(Δnon)); bio_vals = collect(values(Δbio))
        spread_non = maximum(non_vals) - minimum(non_vals)
        spread_bio = maximum(bio_vals) - minimum(bio_vals)
        Δspread = spread_bio - spread_non

        Dmax = maximum(bio_vals) - maximum(non_vals)

        push!(xs, Δspread); push!(ys, Dmax)
        push!(labs, "A=$(String(Alev)), B=$(String(Blev))")
    end

    fig = Figure(; size=(900, 600))
    Label(fig[0,:], @sprintf("Sanity metrics at loss=%.2f — M=%s", loss_pick, String(M_level));
          fontsize=18, padding=(0,0,8,0))
    ax = Axis(fig[1,1], xlabel="Δspread = spread_bio − spread_non",
                        ylabel="Dmax = max_g Δ_bio − max_g Δ_non",
                        title="Positive (down) ⇒ biotic increases worst-geometry loss / geometry dependence")
    scatter!(ax, xs, ys; color=:black, markersize=8)
    for i in eachindex(xs)
        text!(ax, xs[i], ys[i]; text=labs[i], align=(:left,:bottom), fontsize=9)
    end
    vlines!(ax, [0.0]; color=:gray, linestyle=:dash)
    hlines!(ax, [0.0]; color=:gray, linestyle=:dash)

    display(fig)
    return (; fig, table=(Δspread=xs, Dmax=ys, label=labs))
end

# 4) Sanity metrics panel (Δspread and Dmax)
_ = spread_metrics_panel(; grid,
      S=150, basal_frac=0.45, M_level=:on, loss_pick=0.60,
      A_levels=(:neutral,:intermediate,:divergent),
      B_levels=(:none,:soft,:strong),
      seeds_pool=1:6, seeds_mask=1:6, sim_seed=1234,
      nseeds_cluster=6, front_axis=:x, front_noise=0.04,
      τA=0.5, τocc=0.35)