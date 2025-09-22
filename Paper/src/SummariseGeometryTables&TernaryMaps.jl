# ========= TABLE: summarize geometry effects =========
"""
make_geometry_effects_table!(iofile; scenario_name, fstar, boot, weights)

Writes/returns a one-row DataFrame with:
- mean_abs_RF, lo_abs_RF, hi_abs_RF
- mean_abs_RC, lo_abs_RC, hi_abs_RC
- share_R_gt_F, share_R_gt_C
"""
function make_geometry_effects_row(; scenario_name::String, fstar::Float64, boot, weights)
    # Mean absolute deltas across mixtures
    μRF, loRF, hiRF = boot.dRF_mean, boot.dRF_lo, boot.dRF_hi
    μRC, loRC, hiRC = boot.dRC_mean, boot.dRC_lo, boot.dRC_hi

    mean_abs_RF = mean(abs.(μRF))
    lo_abs_RF   = mean(abs.(loRF))
    hi_abs_RF   = mean(abs.(hiRF))

    mean_abs_RC = mean(abs.(μRC))
    lo_abs_RC   = mean(abs.(loRC))
    hi_abs_RC   = mean(abs.(hiRC))

    # Share of mixtures where Δ>0 AND CI excludes 0 (i.e., lo>0)
    share_R_gt_F = mean(boot.dRF_lo .> 0.0)
    share_R_gt_C = mean(boot.dRC_lo .> 0.0)

    DataFrame(; scenario=scenario_name, fstar=round(fstar,digits=2),
              mean_abs_RF=mean_abs_RF, lo_abs_RF=lo_abs_RF, hi_abs_RF=hi_abs_RF,
              mean_abs_RC=mean_abs_RC, lo_abs_RC=lo_abs_RC, hi_abs_RC=hi_abs_RC,
              share_R_gt_F=share_R_gt_F, share_R_gt_C=share_R_gt_C)
end

"""
build_geometry_effects_table(scenarios, fstars; A_seeds, pool, grid, weights)

Runs the bootstrap for each (scenario, f*) pair and collects rows into a table.
Saves CSV and also returns the DataFrame.
"""
function build_geometry_effects_table(scenarios, fstars; A_seeds=1:24, pool, grid, weights)
    results = Vector{DataFrame}(undef, Threads.nthreads())

    Threads.@threads for tid in 1:Threads.nthreads()
        local_rows = DataFrame()
        for (sname, pars, A_fn, agg, kreq) in scenarios
            for fstar in fstars
                boot = Metrics.bootstrap_mixture_maps(;
                    rng=MersenneTwister(11 + tid),   # thread-local RNG seed
                    pool, grid, pars, weights, fstar,
                    A_fn=A_fn, agg=agg, kreq=kreq, A_seeds=A_seeds
                )
                row = make_geometry_effects_row(; scenario_name=sname, fstar=fstar, boot=boot, weights=weights)
                append!(local_rows, row)
            end
        end
        results[tid] = local_rows
    end

    rows = vcat(results...)  # merge all thread-local results

    mkpath("Paper/figs/ternary")
    CSV.write("Paper/figs/ternary/geometry_effects_summary.csv", rows)

    return rows
end

# ========= FILLED TERNARY MAPS =========
# Simple inverse-distance interpolation from scattered mixtures to a lattice
function _idw_interp(x::Vector{Float64}, y::Vector{Float64}, z::Vector{Float64}, xp, yp; p=2.0)
    wsum = 0.0; zsum = 0.0
    for i in eachindex(x)
        dx = xp - x[i]; dy = yp - y[i]
        d2 = dx*dx + dy*dy
        if d2 < 1e-9
            return z[i]     # exact hit
        end
        wi = 1.0 / (d2^(p/2))
        wsum += wi; zsum += wi * z[i]
    end
    return zsum / (wsum + 1e-12)
end

# Build an equilateral ternary lattice and polygons
function filled_ternary_map(weights::Vector{NTuple{3,Float64}},
                            values::Vector{Float64};
                            title="filled map", colormap=:balance, clim=nothing,
                            res::Int=28, fname::String="filled.png")

    # known samples → (x,y)
    xs = Float64[]; ys = Float64[]
    for (a,b,m) in weights
        c = TernaryPlots.ternary_coords(a,b,m)
        push!(xs, c.x); push!(ys, c.y)
    end
    zs = values

    # lattice in barycentric space, then to (x,y)
    step = 1.0/res
    Agrid = Float64[]; Bgrid = Float64[]; X = Float64[]; Y = Float64[]; Z = Float64[]
    for i in 0:res
        a = i*step
        for j in 0:(res-i)
            b = j*step
            m = 1 - a - b
            c = TernaryPlots.ternary_coords(a,b,m)
            push!(Agrid,a); push!(Bgrid,b); push!(X,c.x); push!(Y,c.y)
            push!(Z, _idw_interp(xs, ys, zs, c.x, c.y))
        end
    end

    fig = Figure(; size=(800,720))
    ax  = Axis(fig[1,1]; title=title, xticksvisible=false, yticksvisible=false)
    ax.aspect = DataAspect()
    # draw triangle frame
    TernaryPlots.draw_frame!(ax)

    # Fill by little triangles: for each (i,j) cell we split into two tris
    # Mapping (i,j) in barycentric grid to flat index:
    index(i,j) = (i*(2res+3-i))÷2 + j + 1   # closed-form row indexing
    for i in 0:res-1
        for j in 0:(res-1-i)
            # corners (a,b) indices
            k00 = index(i,j)
            k10 = index(i+1,j)
            k01 = index(i,j+1)
            # poly 1: (i,j) → (i+1,j) → (i,j+1)
            poly!(ax, Point2f[(X[k00],Y[k00]), (X[k10],Y[k10]), (X[k01],Y[k01])];
                 color=mean((Z[k00],Z[k10],Z[k01])), colormap=colormap, strokewidth=0)
            if j+1 <= res-1-i
                # the second triangle exists only if the square is complete
                k11 = index(i+1,j+1)
                poly!(ax, Point2f[(X[k10],Y[k10]), (X[k11],Y[k11]), (X[k01],Y[k01])];
                     color=mean((Z[k10],Z[k11],Z[k01])), colormap=colormap, strokewidth=0)
            end
        end
    end

    # Colorbar: sample a scatter to register a color scale
    sc = scatter!(ax, [0.0], [0.0]; markersize=0.1, color=Z, colormap=colormap)
    if isnothing(clim); clim = (minimum(Z), maximum(Z)); end
    sc[1].colorrange[] = clim
    Colorbar(fig[1,2], sc, label=title)

    xlims!(ax, -0.05, 1.05); ylims!(ax, -0.05, √3/2 + 0.05)
    save(fname, fig); display(fig)
end

# 1) Build the table (writes CSV and returns a DataFrame)
tbl = build_geometry_effects_table(scenarios, fstars;
        A_seeds=1:24, pool=pool_mid, grid=grid_grad, weights=weights)
first(tbl, 6)  # peek in REPL
# CSV -> Paper/figs/ternary/geometry_effects_summary.csv

# 2) Make filled ternaries for one scenario + f*
sname, pars, A_fn, agg, kreq = scenarios[1]      # e.g., :Baseline
fstar = 0.6
boot = Metrics.bootstrap_mixture_maps(; rng=MersenneTwister(11),
            pool=pool_mid, grid=grid_grad, pars, weights, fstar,
            A_fn=A_fn, agg=agg, kreq=kreq, A_seeds=1:24)

# Elasticity filled maps (one per geometry)
for g in (:random, :clustered, :front)
    filled_ternary_map(weights, boot.elastic_mean[g];
        title="Elasticity ∂y/∂f at f*=$(fstar) — $(g)",
        colormap=:plasma,
        fname="Paper/figs/ternary/filled_elasticity_$(sname)_$(g)_f$(Int(100fstar)).png")
end

# Δ filled maps
maxabs = maximum(abs, vcat(boot.dRF_mean, boot.dRC_mean))
filled_ternary_map(weights, boot.dRF_mean;
    title="Δ (Random − Front) at f*=$(fstar) — $(sname)",
    colormap=:balance, clim=(-maxabs, maxabs),
    fname="Paper/figs/ternary/filled_delta_RF_$(sname)_f$(Int(100fstar)).png")

filled_ternary_map(weights, boot.dRC_mean;
    title="Δ (Random − Clustered) at f*=$(fstar) — $(sname)",
    colormap=:balance, clim=(-maxabs, maxabs),
    fname="Paper/figs/ternary/filled_delta_RC_$(sname)_f$(Int(100fstar)).png")
