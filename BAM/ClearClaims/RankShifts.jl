# helper: worst geometry at a given f (uses same baseline at f=0)
function _worst_geom_at_f(; grid::Grid, S::Int, basal_frac::Float64,
        A_level::Symbol, B_level::Symbol, M_level::Symbol,
        loss_pick::Float64, seeds_pool, seeds_mask, sim_seed::Int,
        nseeds_cluster::Int, front_axis::Symbol, front_noise::Float64,
        τA::Float64, τocc::Float64)

    geoms = (:random, :clustered, :front)

    # Δ(f) := BSH(f) - BSH(0) for AM and ABM, by geometry
    Δnon = Dict{Symbol,Float64}();  Δbio = Dict{Symbol,Float64}()

    for hk in geoms
        c_f  = _sweep_AM_ABM_onegeom(; grid,S,basal_frac,A_level,
                                     B_level, M_level,
                                     hl_kind=hk, loss_fracs=[loss_pick],
                                     seeds_pool,seeds_mask,sim_seed,
                                     nseeds_cluster,front_axis,front_noise,τA,τocc,
                                     T_frac_on=0.0)

        c_0  = _sweep_AM_ABM_onegeom(; grid,S,basal_frac,A_level,
                                     B_level, M_level,
                                     hl_kind=hk, loss_fracs=[0.0],
                                     seeds_pool,seeds_mask,sim_seed,
                                     nseeds_cluster,front_axis,front_noise,τA,τocc,
                                     T_frac_on=0.0)

        Δnon[hk] = c_f.AM[1]  - c_0.AM[1]     # climate-only change
        Δbio[hk] = c_f.ABM[1] - c_0.ABM[1]    # biotic change
    end

    # choose worst (most negative) with a *vector aligned* to geoms
    v_non = [Δnon[g] for g in geoms]
    v_bio = [Δbio[g] for g in geoms]
    worst_off_g = geoms[argmin(v_non)]
    worst_on_g  = geoms[argmin(v_bio)]

    return (; worst_off_g, worst_on_g, Δnon, Δbio)
end

function rank_shift_grid(; grid::Grid, S::Int=150, basal_frac::Float64=0.45,
        M_level::Symbol=:on, loss_pick::Float64=0.6,
        A_levels=(:neutral,:intermediate,:divergent),
        B_levels=(:none,:soft,:strong),
        seeds_pool=1:6, seeds_mask=1:6, sim_seed::Int=1234,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        τA::Float64=0.5, τocc::Float64=0.35)

    geoms = (:random,:clustered,:front)

    fig = Figure(; size=(1100, 720))
    Label(fig[0,:], @sprintf("Rank shifts at loss=%.2f — M=%s", loss_pick, String(M_level));
          fontsize=18, padding=(0,0,8,0))

    # grid of small axes (rows=A, cols=B)
    gl = fig[1,1] = GridLayout()
    for (i,Alev) in enumerate(A_levels)
        for (j,Blev) in enumerate(B_levels)
            ax = Axis(gl[i,j]; width=260, height=180,
                      limits=(0,1,0,1))  # fixed canvas
            hidedecorations!(ax); hidespines!(ax)
            text!(ax, 0.04, 0.94,
                  text = @sprintf("A=%s, B=%s", String(Alev), String(Blev)),
                  space=:relative, align=(:left,:top), fontsize=11)

            w = _worst_geom_at_f(; grid,S,basal_frac,
                    A_level=Alev, B_level=Blev, M_level,
                    loss_pick, seeds_pool,seeds_mask,sim_seed,
                    nseeds_cluster,front_axis,front_noise,τA,τocc)

            # left: B OFF; right: B ON
            scatter!(ax, [0.35],[0.55]; marker=G_MARK[w.worst_off_g],
                    color=A_COL[Alev], strokecolor=:black,
                    strokewidth=B_STROKE[Blev], markersize=16)
            scatter!(ax, [0.65],[0.55]; marker=G_MARK[w.worst_on_g],
                    color=A_COL[Alev], strokecolor=:black,
                    strokewidth=B_STROKE[Blev], markersize=16)
            text!(ax, 0.35, 0.30, text="B OFF", align=(:center,:top), fontsize=10)
            text!(ax, 0.65, 0.30, text="B ON",  align=(:center,:top), fontsize=10)

            if w.worst_off_g != w.worst_on_g
                arrows!(ax, Point2f(0.42,0.55) => Point2f(0.58,0.55);
                        linewidth=2, arrowsize=9)
            end
        end
    end

    # compact legend
    leg = fig[1,2] = GridLayout()
    Label(leg[1,1], "Geometry markers:", halign=:left)
    for (r,(g,mk)) in enumerate(G_MARK)
        axg = Axis(leg[1+r,1]; width=26, height=26)
        hidedecorations!(axg); hidespines!(axg)
        scatter!(axg, [0],[0]; marker=mk, markersize=12, color=:gray)
        Label(leg[1+r,2], _glabel[g]; halign=:left)
    end
    Label(leg[5,1], "A color / B stroke encode A,B", halign=:left)

    display(fig)
    return fig
end

_ = rank_shift_grid(
    ; grid,
    S=150, basal_frac=0.45, M_level=:on, loss_pick=0.40,
    A_levels=(:neutral,:intermediate,:divergent),
    B_levels=(:none,:soft,:strong),
    seeds_pool=1:6, seeds_mask=1:6, sim_seed=1234,
    nseeds_cluster=6, front_axis=:x, front_noise=0.04,
    τA=0.5, τocc=0.35
)