curves = Metrics.abs_bsh_vs_loss(; rng, pool=pool_low, grid=grid_mosa, pars=pars,
                                 A_fn=A_fn_punch, loss_fracs=0.2:0.1:0.8, seed_A=1)
figABS = Figs.fig_abs_bsh_vs_loss(curves; title="Absolute BSH vs loss — gradient")
save("Paper/figs/ABS_BSH_vs_loss.png", figABS)

begin 
    curves = Metrics.abs_bsh_vs_loss(; rng, pool=pool_low, grid=grid_grad, pars=pars, loss_fracs=0.2:0.05:0.8, seed_A=1)
    fig = Figure(size=(900,300))
    for (k,geom) in enumerate((:front,:clustered,:random))
        ax = Axis(fig[1,k], title=String(geom), xlabel="area lost", ylabel="BAM − AM")
        Δ = curves[geom].BAM .- curves[geom].AM
        lines!(ax, curves[geom].loss, Δ)
        hlines!(ax, [0], color=(:gray,0.5), linestyle=:dash)
    end
    display(fig)
end

###############################

realA = run_realities_Aonly(; rng, grid=grid_grad,
         loss_fracs=0.2:0.1:0.8, S=S, basal_frac=basal_frac, seed_A=1,
         Npools=12, geoms=(:random,:clustered,:front), τA=pars.τA)

figR = Figs.fig_realities_Aonly(realA; title="A-only response under ABM / MAB / BAM realities")
save("Paper/figs/Realities_Aonly.png", figR)

###################################
# choose “dominance” strengths; tweak to taste
pars_ABM = BSH.BAMParams(; τA=0.6, movement=:off)
pars_MAB = BSH.BAMParams(; τA=0.50, movement=:component, T=18)
pars_BAM = BSH.BAMParams(; τA=0.50, τB=0.65, movement=:off)

realSuit = MetawebSweep.run_realities_suitable_area(
    ; rng, grid=grid_grad, loss_fracs=0.2:0.05:0.8, S=S, basal_frac=basal_frac,
      Npools=12, seed_A=1, geoms=(:front,:clustered,:random),
      pars_ABM=pars_ABM, pars_MAB=pars_MAB, pars_BAM=pars_BAM,
      agg=:kofn, kreq=2
)

# quick plot (reuse your Figs style)
begin
    fig = Figure(; size=(1050,330))
    for (k,geom) in enumerate((:front,:clustered,:random))
        ax = Axis(fig[1,k], title=String(geom),
                xlabel="area lost (fraction)",
                ylabel="suitable area (mean over consumers / original area)",
                ylabelsize=10)
        x = realSuit[geom].loss
        lines!(ax, x, realSuit[geom].ABM; label="ABM")
        lines!(ax, x, realSuit[geom].MAB; label="MAB")
        lines!(ax, x, realSuit[geom].BAM; label="BAM")
        axislegend(ax, position=:rt)
    end
    display(fig)
end
save("Paper/figs/Realities_suitable.png", fig); display(fig)

###########################################################
weights = MetawebSweep.simplex_grid(0.25)  # coarse: 0, .25, .5, .75, 1
blend = MetawebSweep.run_blend_realities_suitable_area(
    ; rng, grid=grid_grad, S=S, basal_frac=basal_frac,
    loss_fracs=0.2:0.1:0.8, geoms=(:random,:clustered,:front),
    pars_ABM=BAMParams(; τA=0.6, movement=:off),
    pars_MAB=BAMParams(; τA=0.50, movement=:component, T=18),
    pars_BAM=BAMParams(; τA=0.50, τB=0.65, movement=:off),
    weights=weights, Nrep=10
)

# Example: plot a few slices
begin
    fig = Figure(; size=(900,320))
    picks = weights
    hl_scenario = :front

    # helper to find nearest available key
    function nearest_key(dict::Dict, target::Tuple{Float64,Float64,Float64})
        keys_list = collect(keys(dict))
        dists = [sum(abs.(k .- target)) for k in keys_list]
        return keys_list[argmin(dists)]
    end

    ax = Axis(fig[1,1], title="$(hl_scenario)", xlabel="area lost", ylabel="suitable / area0")

    for w in picks
        k = nearest_key(blend[hl_scenario], w)
        lines!(ax, blend[:random][k].loss, blend[hl_scenario][k].y;
               label="A=$(w[1]) M=$(w[2]) B=$(w[3])")
    end

    axislegend(ax, position=:rt, labelsize=10)
    display(fig)
end

save("Paper/figs/Blend_random.png", fig); display(fig)

############################################################
function fig_blend_triangles(blend; fstar=0.6, geoms=(:random,:clustered,:front))
    # barycentric → 2D triangle coords
    function bary_to_cart(w)
        wA,wM,wB = w
        x = 0.5 * (2wM + wB) / (wA+wM+wB)
        y = (√3/2) * wB / (wA+wM+wB)
        return (x,y)
    end

    # collect global value range across all geoms
    all_vals = Float64[]
    for g in geoms, res in values(blend[g])
        if !isempty(res.loss) && !isempty(res.y)
            idx = argmin(abs.(res.loss .- fstar))
            push!(all_vals, res.y[idx])
        end
    end
    vmin, vmax = extrema(all_vals)

    fig = Figure(; size=(1600,600))

    for (j,g) in enumerate(geoms)
        ax = Axis(fig[1,j]; title=string(g), aspect=1,
                  xgridvisible=false, ygridvisible=false)
        hidedecorations!(ax)

        xs, ys, vals = Float64[], Float64[], Float64[]
        for (w,res) in blend[g]
            if isempty(res.loss) || isempty(res.y)
                continue
            end
            idx = argmin(abs.(res.loss .- fstar))
            x,y = bary_to_cart(w)
            push!(xs,x); push!(ys,y); push!(vals,res.y[idx])
        end

        # Triangle outline
        poly!(ax, Point2f[(0,0),(1,0),(0.5,√3/2)];
              color=:transparent, strokecolor=:black, strokewidth=1.5)

        scatter!(ax, xs, ys; color=vals, colormap=:viridis,
                 markersize=28, colorrange=(vmin,vmax))

        # Corner labels
        text!(ax, 0,0; text="A", align=(:right,:top), fontsize=18)
        text!(ax, 1,0; text="M", align=(:left,:top), fontsize=18)
        text!(ax, 0.5,√3/2; text="B", align=(:center,:bottom), fontsize=18)
    end

    Colorbar(fig[1,4]; colormap=:viridis, limits=(vmin,vmax),
             label="Retained suitable area")

    Label(fig[0,:], "Retained suitable area at f* = $fstar",
          fontsize=24, tellwidth=false)

    return fig
end

fig = fig_blend_triangles(blend)

############################################
"""
fig_blend_delta(Δrf; title="Random vs Front", clim=nothing)

Plot a ternary map where color = delta = (random – front) at f*.
Positive ⇒ random retains more suitable area (less damaging).
Negative ⇒ random loses more (worse than front).
"""
function fig_blend_delta(Δrf; title="Random – Front", clim=nothing)
    # barycentric → 2D triangle coords
    bary_to_cart(wA,wM,wB) = (
        0.5 * (2wM + wB) / (wA+wM+wB),
        (√3/2) * wB / (wA+wM+wB)
    )

    xs, ys, vals = Float64[], Float64[], Float64[]
    for (wA,wM,wB,δ) in Δrf
        x,y = bary_to_cart(wA,wM,wB)
        push!(xs,x); push!(ys,y); push!(vals,δ)
    end

    vmin, vmax = isnothing(clim) ? extrema(vals) : clim

    fig = Figure(; size=(600,550))
    ax  = Axis(fig[1,1]; aspect=1, title=title,
               xgridvisible=false, ygridvisible=false)
    hidedecorations!(ax)

    # Triangle outline
    poly!(ax, Point2f[(0,0),(1,0),(0.5,√3/2)];
          color=:transparent, strokecolor=:black, strokewidth=1.5)

    # Scatter points
    scatter!(ax, xs, ys; color=vals, colormap=:balance,
             markersize=28, colorrange=(vmin,vmax))

    # Corner labels
    text!(ax, 0,0; text="A", align=(:right,:top), fontsize=18)
    text!(ax, 1,0; text="M", align=(:left,:top), fontsize=18)
    text!(ax, 0.5,√3/2; text="B", align=(:center,:bottom), fontsize=18)

    Colorbar(fig[1,2]; colormap=:balance, limits=(vmin,vmax),
             label="Δ (Random – Front)")

    return fig
end

Δrf = MetawebSweep.blend_geometry_delta_at_fstar(blend, 0.6; geomA=:random, geomB=:clustered)
figΔ = fig_blend_delta(Δrf; title="Random vs Clustered at f* = 0.6")

#############################################
"""
fig_blend_elasticity(blend_elast; geom=:random, fstar=0.6, clim=nothing)

Plot a ternary map where color = slope (dY/df) at f* for a given geometry.
"""
function fig_blend_elasticity(blend_elast::Dict{Symbol,<:Any};
                              geom::Symbol=:random, fstar::Float64=0.6,
                              clim=nothing)

    # pick geometry’s slope tuples
    data = blend_elast[geom]

    # barycentric → 2D
    bary_to_cart(wA,wM,wB) = (
        0.5 * (2wM + wB) / (wA+wM+wB),
        (√3/2) * wB / (wA+wM+wB)
    )

    xs, ys, vals = Float64[], Float64[], Float64[]
    for (wA,wM,wB,slope) in data
        x,y = bary_to_cart(wA,wM,wB)
        push!(xs,x); push!(ys,y); push!(vals,slope)
    end

    vmin, vmax = isnothing(clim) ? extrema(vals) : clim

    fig = Figure(; size=(600,550))
    ax  = Axis(fig[1,1]; aspect=1,
               title="Elasticity at f*=$(fstar), geom=$(String(geom))",
               xgridvisible=false, ygridvisible=false)
    hidedecorations!(ax)

    # Triangle outline
    poly!(ax, Point2f[(0,0),(1,0),(0.5,√3/2)];
          color=:transparent, strokecolor=:black, strokewidth=1.5)

    # Scatter points
    scatter!(ax, xs, ys; color=vals, colormap=:plasma,
             markersize=28, colorrange=(vmin,vmax))

    # Corner labels
    text!(ax, 0,0; text="A", align=(:right,:top), fontsize=18)
    text!(ax, 1,0; text="M", align=(:left,:top), fontsize=18)
    text!(ax, 0.5,√3/2; text="B", align=(:center,:bottom), fontsize=18)

    Colorbar(fig[1,2]; colormap=:plasma, limits=(vmin,vmax),
             label="∂y/∂f at f*")

    return fig
end

elast = blend_elasticity_at_fstar(blend, 0.6)
figE = fig_blend_elasticity(elast; geom=:front, fstar=0.6)
