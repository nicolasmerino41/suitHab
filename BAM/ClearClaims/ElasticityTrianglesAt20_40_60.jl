const A_COL = Dict(:neutral=>RGBf(0.55,0.55,0.55),
                   :intermediate=>RGBf(0.25,0.45,0.85),
                   :divergent=>RGBf(0.85,0.35,0.25))
const G_MARK = Dict(:random=>:circle, :clustered=>:rect, :front=>:utriangle)
const B_STROKE = Dict(:none=>0.6, :soft=>1.4, :strong=>2.4)

"""
elasticity_triangles_all_scenarios(; grid, S, basal_frac, M_level,
                                   A_levels, B_levels, loss_list,
                                   seeds_pool, seeds_mask, sim_seed,
                                   nseeds_cluster, front_axis, front_noise, τA, τocc)
"""
function elasticity_triangles_all_scenarios(; grid::Grid, S::Int=150, basal_frac::Float64=0.45,
        M_level::Symbol=:on, A_levels=(:neutral,:intermediate,:divergent),
        B_levels=(:none,:soft,:strong), loss_list=[0.2,0.4,0.6],
        seeds_pool=1:4, seeds_mask=1:4, sim_seed::Int=1234,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        τA::Float64=0.5, τocc::Float64=0.35)

    geoms = (:random,:clustered,:front)

    function ternary_frame!(ax, title)
        xs=[0.0,1.0,0.5,0.0]; ys=[0.0,0.0,√3/2,0.0]
        lines!(ax, xs, ys; color=:black)
        text!(ax, 0.02, -0.04, text="Abiotic", align=(:left,:bottom))
        text!(ax, 0.98, -0.04, text="Biotic",  align=(:right,:bottom))
        text!(ax, 0.50,  √3/2+0.03, text="Movement", align=(:center,:bottom))
        hidespines!(ax, :r, :t); ax.title=title
    end

    fig = Figure(; size=(1500, 600))
    Label(fig[0,:], "Elasticity shares across scenarios — M=$(String(M_level))"; fontsize=18, padding=(0,0,8,0))

    for (col, f) in enumerate(loss_list)
        ax = Axis(fig[1,col])
        ternary_frame!(ax, @sprintf("loss = %.2f", f))
        keepfrac = 1 - f

        for Alev in A_levels, Blev in B_levels, hk in geoms
            # shares at this f for this (A,B) and geometry
            shares = geometry_elasticity_shares_all_geoms(; grid,S,basal_frac,
                        A_level=Alev, B_level=Blev, M_level,
                        loss_pick=f, seeds_pool,seeds_mask,sim_seed,
                        nseeds_cluster,front_axis,front_noise,τA,τocc)
            s = shares[hk]
            x,y = (s.S_B + 0.5*s.S_M, (√3/2)*s.S_M)
            scatter!(ax, [x],[y];
                     color=A_COL[Alev], marker=G_MARK[hk],
                     strokecolor=:black, strokewidth=B_STROKE[Blev],
                     markersize=10)
        end
    end

    # custom legend
    lg = fig[1, length(loss_list)+1] = GridLayout()
    # A-level colors
    Label(lg[1,1], "A-levels:", halign=:left)
    for (i,(a,c)) in enumerate(A_COL)
        axc = Axis(lg[1+i,1]; width=22, height=22); hidedecorations!(axc); hidespines!(axc)
        scatter!(axc, [0],[0]; color=c, marker=:circle, markersize=10)
        Label(lg[1+i,2], String(a); halign=:left)
    end
    # Geometry markers
    rowoff = 1 + length(A_COL) + 1
    Label(lg[rowoff,1], "Geometry:", halign=:left)
    for (j,(g,mk)) in enumerate(G_MARK)
        axg = Axis(lg[rowoff+j,1]; width=22, height=22); hidedecorations!(axg); hidespines!(axg)
        scatter!(axg, [0],[0]; color=:gray, marker=mk, markersize=10)
        Label(lg[rowoff+j,2], _glabel[g]; halign=:left)
    end
    # B-level stroke
    rowoff2 = rowoff + length(G_MARK) + 1
    Label(lg[rowoff2,1], "B-levels:", halign=:left)
    for (k,(b,sw)) in enumerate(B_STROKE)
        axb = Axis(lg[rowoff2+k,1]; width=22, height=22); hidedecorations!(axb); hidespines!(axb)
        scatter!(axb, [0],[0]; color=:white, strokecolor=:black, strokewidth=sw, markersize=12)
        Label(lg[rowoff2+k,2], String(b); halign=:left)
    end

    display(fig)
    return fig
end

_ = elasticity_triangles_all_scenarios(; grid,
      S=150, basal_frac=0.45, M_level=:on,
      A_levels=(:neutral,:intermediate,:divergent),
      B_levels=(:none,:soft,:strong),
      loss_list=[0.2,0.4,0.6],
      seeds_pool=1:4, seeds_mask=1:4, sim_seed=1234,
      nseeds_cluster=6, front_axis=:x, front_noise=0.04,
      τA=0.5, τocc=0.35)