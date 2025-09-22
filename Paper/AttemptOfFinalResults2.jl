#### RealityMix.jl (paste in your driver or a tiny src file you include) ####
using Random, Statistics
using CairoMakie
using .Grids, .HL, .BSH
mkpath("Paper/AttemptOfFinalResults2")
# --- barycentric <-> plot coords (equilateral triangle, base = A–M, top = B)
_tri_x(a,b,m) = a*0.0 + b*0.5 + m*1.0        # unnormalized; we’ll scale
_tri_y(a,b,m) = a*0.0 + b*√3/2 + m*0.0

"Even grid of integer barycentric coordinates with step K (e.g., K=6 → 28 pts)."
function triangle_grid(K::Int)
    pts = NTuple{3,Float64}[]
    for i in 0:K, j in 0:(K-i)
        k = K - i - j
        push!(pts, (i/K, j/K, k/K))  # (A,B,M)
    end
    pts
end

"Evaluate A-only, AM, BAM responses at loss f for a keep mask geometry."
function _responses_at_f(; rng, pool, grid, pars, A_seed::Int, geom::Symbol, f::Float64)
    keepfrac = 1 - f
    keep = geom === :random    ? HL.random_mask(rng, grid.C, keepfrac) :
           geom === :clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                 HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)
    # one abiotic field (like everywhere else)
    A = BSH.abiotic_matrix(pool, grid; seed=A_seed)

    yA = BSH.mean_Aonly_over_consumers_area0(view(A, :, keep), pool, grid.C)
    yAM = BSH.mean_BSH_over_consumers_area0(rng, pool, grid, pars; A=A, keepmask=keep, mode=:AM)
    yBM = BSH.mean_BSH_over_consumers_area0(rng, pool, grid, pars; A=A, keepmask=keep, mode=:BAM)
    (; yA, yAM, yBM)
end

"Mixture response y = wA*yA + wB*yBM + wM*yAM."
_mixture_y(yA, yBM, yAM, wA, wB, wM) = wA*yA + wB*yBM + wM*yAM

"""
Bootstrap the mixture field at a single f (e.g., f*):
- weights :: vector of (A,B,M) tuples
- A_seeds :: range (e.g., 1:16) used as bootstrap draws
Returns NamedTuple with mean & CI (per weight) for each geometry.
"""
function mixture_field_boot(; rng, pool, grid, pars, weights::Vector{NTuple{3,Float64}},
                            geoms::Tuple{Vararg{Symbol}}, f::Float64,
                            A_seeds = 1:12, alpha::Float64=0.10)
    out = Dict{Symbol,NamedTuple}()
    for geom in geoms
        vals = [Float64[] for _ in 1:length(weights)]
        for seedA in A_seeds
            resp = _responses_at_f(; rng, pool, grid, pars, A_seed=seedA, geom, f)
            for (i,(wA,wB,wM)) in enumerate(weights)
                y = _mixture_y(resp.yA, resp.yBM, resp.yAM, wA,wB,wM)
                push!(vals[i], y)
            end
        end
        μ  = map(mean, vals)
        lo = map(v -> quantile(v, alpha/2),  vals)
        hi = map(v -> quantile(v, 1-alpha/2), vals)
        out[geom] = (; mean=μ, lo=lo, hi=hi)
    end
    out
end

"Finite-diff elasticity at f: slope dy/df using ±δf around f."
function mixture_elasticity_boot(; rng, pool, grid, pars, weights, geoms, f::Float64,
                                 δf::Float64=0.02, A_seeds=1:12)
    fld_lo = mixture_field_boot(; rng, pool, grid, pars, weights, geoms, f=f-δf, A_seeds)
    fld_hi = mixture_field_boot(; rng, pool, grid, pars, weights, geoms, f=f+δf, A_seeds)
    out = Dict{Symbol,Vector{Float64}}()
    for g in geoms
        e = (fld_hi[g].mean .- fld_lo[g].mean) ./ (2δf)
        out[g] = e
    end
    out
end

"Convenience: two-geometry contrast (Random–Front or Clustered–Front)."
contrast_field(mean1::Vector{Float64}, mean2::Vector{Float64}) = mean1 .- mean2

# ---- plotting helpers (filled triangles via fast IDW) -----------------------

"Map barycentric (A,B,M) to 2D plot coords."
function bary_to_xy(w::NTuple{3,Float64})
    a,b,m = w
    return _tri_x(a,b,m), _tri_y(a,b,m)
end

"Interpolate scattered triangle values (x,y,v) onto a dense mesh for a 'filled' look."
function idw_fill(xs, ys, vs; res::Int=200, p::Float64=2.0)
    xmin, xmax = minimum(xs), maximum(xs)
    ymin, ymax = minimum(ys), maximum(ys)
    gx = range(xmin, xmax; length=res)
    gy = range(ymin, ymax; length=res)
    Z  = zeros(length(gx), length(gy))
    for (i,xx) in enumerate(gx), (j,yy) in enumerate(gy)
        wsum = 0.0; vsum = 0.0
        for (x,y,v) in zip(xs,ys,vs)
            d = hypot(xx-x, yy-y)
            w = d < 1e-6 ? 1e9 : 1/(d^p)
            wsum += w; vsum += w*v
        end
        Z[i,j] = vsum / wsum
    end
    (; gx, gy, Z)
end

"Paint a filled triangle from weights + values. Returns the Makie figure."
function triangle_filled_plot(weights, values; title::String="", colormap=:viridis,
                              clim=nothing, res::Int=220, size=(540,520))
    xs = Float64[]; ys = Float64[]; for w in weights
        x,y = bary_to_xy(w); push!(xs,x); push!(ys,y)
    end
    fld = idw_fill(xs, ys, values; res=res)

    fig = Figure(; size=size)
    ax  = Axis(fig[1,1], title=title, xticksvisible=false, yticksvisible=false,
               xlabel="A … M", ylabel="", xgridvisible=false, ygridvisible=false)
    heatmap!(ax, fld.gx, fld.gy, fld.Z; colormap, colorrange=clim)
    # draw triangle edges + vertices
    tri = [(0.0,0.0), (0.5,√3/2), (1.0,0.0), (0.0,0.0)]
    lines!(ax, first.(tri), last.(tri), color=:black, linewidth=1.2)
    text!(ax, 0.0, 0.0, text="A", align=(:left,:top), fontsize=14)
    text!(ax, 1.0, 0.0, text="M", align=(:right,:top), fontsize=14)
    text!(ax, 0.5, √3/2, text="B", align=(:center,:bottom), fontsize=14)
    Colorbar(fig[1,2])
    fig
end

#### DRIVERMIX.jl — main A–B–M figure maker ####
# === shared context (reuse your setup) ===
rng = MersenneTwister(123)
nx, ny = 120, 120
S = 200; basal_frac = 0.25

G = GridsInit.build_default_grids(nx, ny; seeds=(11,12,13,14))
grid = G.grad                   # headline grid; change if desired
pool = Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:mid)

pars = BSH.BAMParams(; τA=0.5, τB=0.5, movement=:component, T=8)
loss_fracs = 0.2:0.05:0.8
best = Metrics.ks_best_fstar(; rng, pool, grid, pars, loss_fracs, geometry=:front, seed_A=1)
fstar = best.f
println("Using f* = ", fstar)

# === mixture design: a regular grid inside the triangle (balanced coverage) ===
weights = triangle_grid(6)   # 6 → 28 points (fast); 8 → 45; 10 → 66

# === core bootstrap once; reuse everywhere ===
A_seeds = 1:12
fld = mixture_field_boot(; rng, pool, grid, pars, weights,
                         geoms=(:random,:clustered,:front), f=fstar, A_seeds)

# ---------- F-A: Reality triangle per geometry (filled, mean) ----------
for g in (:random,:clustered,:front)
    μ = fld[g].mean
    clim = (minimum(μ), maximum(μ))
    fig = triangle_filled_plot(weights, μ;
        title="Retained suitable area at f*=$(round(fstar,digits=2)) — $(String(g))",
        colormap=:viridis, clim=clim, res=260, size=(560,540))
    save("Paper/AttemptOfFinalResults2/FA_triangle_$(String(g)).png", fig)
end

# ---------- F-B: Geometry contrasts (Δ Random–Front, Δ Clustered–Front) ----------
dRF = contrast_field(fld[:random].mean,   fld[:front].mean)
dCF = contrast_field(fld[:clustered].mean,fld[:front].mean)
maxabs = maximum(abs, vcat(dRF, dCF))
for (name,dat) in [("Random_minus_Front", dRF), ("Clustered_minus_Front", dCF)]
    fig = triangle_filled_plot(weights, dat;
        title="Δ $(replace(name,'_' => ' ')) at f*=$(round(fstar,digits=2))",
        colormap=:balance, clim=(-maxabs, maxabs), res=260, size=(560,540))
    save("Paper/AttemptOfFinalResults2/FB_contrast_$(name).png", fig)
end

# ---------- F-C: Elasticity triangles (∂y/∂f at f*) ----------
elas = mixture_elasticity_boot(; rng, pool, grid, pars, weights,
                               geoms=(:random,:clustered,:front), f=fstar, δf=0.02,
                               A_seeds=A_seeds)
for g in (:random,:clustered,:front)
    e = elas[g]
    fig = triangle_filled_plot(weights, e;
        title="Elasticity ∂y/∂f at f*=$(round(fstar,digits=2)) — $(String(g))",
        colormap=:plasma, clim=(minimum(e), maximum(e)), res=260, size=(560,540))
    save("Paper/AttemptOfFinalResults2/FC_elasticity_$(String(g)).png", fig)
end

# ---------- F-D: Triangle-strip curves (4 paths × 3 geometries, ribbons) ----------
function pick_path(sym::Symbol)
    if sym === :A_edge      # A→M edge at B=0 (11 points)
        return [(a, 0.0, 1-a) for a in range(0,1; length=11)]
    elseif sym === :B_edge  # B→A (M=0)
        return [(1-b, b, 0.0) for b in range(0,1; length=11)]
    elseif sym === :M_edge  # M→B (A=0)
        return [(0.0, b, 1-b) for b in range(0,1; length=11)]
    else # center ray A→center→M
        return [(0.33,0.34,0.33)]
    end
end

paths = Dict(
    :A_edge => pick_path(:A_edge),
    :B_edge => pick_path(:B_edge),
    :M_edge => pick_path(:M_edge),
    :center => [(0.33,0.34,0.33)]
)

losses = collect(loss_fracs)
begin
    fig = Figure(; size=(1180,680))
    gi = 1
    for (pname, pts) in pairs(paths)
        row = gi; gi += 1
        for (col, geom) in enumerate((:random,:clustered,:front))
            ax = Axis(fig[row,col], title="$(String(pname)) — $(String(geom))",
                    xlabel="area lost f", ylabel="suitable area / original")
            for w in pts
                # fast: bootstrap over A_seeds for each f
                μ = Float64[]; lo = Float64[]; hi = Float64[]
                for f in losses
                    tmp = mixture_field_boot(; rng, pool, grid, pars, weights=[w],
                                            geoms=(geom,), f=f, A_seeds=1:8) # smaller boot here
                    push!(μ,  tmp[geom].mean[1])
                    push!(lo, tmp[geom].lo[1]); push!(hi, tmp[geom].hi[1])
                end
                band!(ax, losses, lo, hi; color=RGBAf(0.2, 0.2, 0.2, 0.15))
                lines!(ax, losses, μ; linewidth=2)
            end
        end
    end
    display(fig)
end
save("Paper/AttemptOfFinalResults2/FD_triangle_strips.png", fig)
