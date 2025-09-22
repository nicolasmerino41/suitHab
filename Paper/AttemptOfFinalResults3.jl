# === run_story.jl ============================================================
# Uses your modules only
include("../SetUp.jl")
include("src/grids.jl");      using .Grids
include("src/metawebs.jl");   using .Metawebs
include("src/hl.jl");         using .HL
include("src/bsh.jl");        using .BSH
include("src/metrics.jl");    using .Metrics
include("src/figs.jl");       using .Figs

# ---------- SETTINGS ----------
rng        = MersenneTwister(123)
nx, ny     = 120, 120
S          = 200
basal_frac = 0.25
fstar      = 0.60              # show at 60% loss
geom_list  = (:front, :clustered, :random)
Npool      = 10                 # metaweb reps for CIs
loss_fracs = 0.2:0.05:0.8       # used in Fig 5
outdir     = "Paper/story"; isdir(outdir) || mkpath(outdir)

# Build grids (single source of truth)
include("src/grids_init.jl"); using .GridsInit
G        = GridsInit.build_default_grids(nx, ny; seeds=(11,12,13,14))
grid     = G.grad               # pick one grid for the story; change if desired

# BAM params (movement ON so M is meaningful)
pars = BSH.BAMParams(; τA=0.5, τB=0.5, τocc=0.2, γ=3.0, movement=:component, T=8)

# ---------- MIXTURE SAMPLER (barycentric) ----------
# mixtures are (wA, wB, wM), wA+wB+wM=1
function mix_points(nside::Int=17)
    pts = Tuple{Float64,Float64,Float64}[]
    for i in 0:nside, j in 0:(nside-i)
        k = nside - i - j
        push!(pts, (i/nside, j/nside, k/nside))
    end
    pts
end

# Reality builder: A–B–M mixture
# - wA: abiotic suitability (A≥τA & kept)
# - wB: biotic coupling (τB / k-of-n effect)
# - wM: movement component constraint (component size ≥ T)
# We create a *realized* occupancy mask P_real per species by blending the three gates.
# Then measure “suitable area” as mean over consumers / original area.
function realized_area_mean(; rng, grid, pars::BSH.BAMParams,
        wA, wB, wM, f::Float64, geom::Symbol, pool::Metawebs.SpeciesPool, seed_A::Int)

    @assert isapprox(wA+wB+wM, 1.0; atol=1e-8)

    # Build abiotic matrix once per pool/seed to keep paths consistent
    A = BSH.abiotic_matrix(pool, grid; seed=seed_A)

    keepfrac = 1 - f
    keep = geom === :random    ? HL.random_mask(rng, grid.C, keepfrac) :
           geom === :clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                 HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)

    # Gating pieces (per species × cell)
    # M gate: component constraint among kept & A≥τA
    M = pars.movement === :off ? trues(pool.S, grid.C) :
        BSH.movement_gate(grid, A, keep; τA=pars.τA, T=pars.T)

    # AM occupancy (A gate + movement)
    P_AM, _ = BSH.assemble_AM(pool, grid, A, keep; pars=pars)

    # BAM occupancy (A + movement + prey sufficiency)
    bam      = BSH.assemble_BAM(pool, grid, A, keep; pars=pars)
    P_BAM    = bam.P

    # “A only” boolean (just climate pass & kept & movement gate)
    P_Aonly  = falses(pool.S, grid.C)
    @inbounds for s in 1:pool.S, i in 1:grid.C
        P_Aonly[s,i] = keep[i] & (A[s,i] ≥ pars.τA) & M[s,i]
    end

    # Blend the three realized maps:
    #   wA -> P_Aonly,  wM -> M∩(A≥τA∩kept),  wB -> P_BAM
    # For M, we need an occupancy-like map; use (A≥τA & kept & M) = P_M
    P_M = P_Aonly
    P_real = falses(pool.S, grid.C)
    @inbounds for s in 1:pool.S, i in 1:grid.C
        val = wA*(P_Aonly[s,i] ? 1.0 : 0.0) +
              wM*(P_M[s,i]     ? 1.0 : 0.0) +
              wB*(P_BAM[s,i]   ? 1.0 : 0.0)
        # realized presence if weighted sum crosses 0.5 (soft mixture → hard occupancy)
        P_real[s,i] = (val ≥ 0.5)
    end

    # Mean over consumers / original area
    cons = findall(!, pool.basal)
    isempty(cons) && return 0.0
    m = mean([sum(@view P_real[s, :]) for s in cons]) / grid.C
    return m
end

# Cache metawebs for CIs
pools = [Metawebs.build_metaweb_archetype(MersenneTwister(10+i); S=S, basal_frac=basal_frac, archetype=:mid)
         for i in 1:Npool]

# Compute surfaces (one pass; reused for Figs 1–3)
function build_surfaces(; fstar, geom_list, nside=25)
    pts   = mix_points(nside)
    surf  = Dict{Symbol, Tuple{Vector{Float64},Vector{Float64},Vector{Float64},Matrix{Float64}}}()
    # coords for triangular grid (for plotting)
    xs = Float64[]; ys = Float64[]
    for (wA,wB,wM) in pts
        # map barycentric → cartesian (equilateral triangle)
        x = wA*0.0 + wB*0.5 + wM*1.0
        y = wA*0.0 + wB*√3/2 + wM*0.0
        push!(xs,x); push!(ys,y)
    end
    for geom in geom_list
        Z = zeros(length(pts))
        for (t,(wA,wB,wM)) in enumerate(pts)
            vals = Float64[]
            for (k,pool) in enumerate(pools)
                push!(vals, realized_area_mean(; rng, grid, pars, wA, wB, wM, f=fstar,
                                               geom, pool, seed_A=100+k))
            end
            Z[t] = mean(vals)
        end
        surf[geom] = (xs, ys, [getindex(p,1) for p in pts], reshape(Z, :, 1))  # keep pts too
    end
    return (pts=pts, surf=surf, xs=xs, ys=ys)
end

@info "Computing surfaces…"
SFC = build_surfaces(; fstar=fstar, geom_list=geom_list, nside=25)

# ---------- Fig 1: Reality surface per geometry ----------
function plot_surface(geom::Symbol; title="")
    xs, ys, _, ZZ = SFC.surf[geom]
    # Make a dense scatter heat (Makie handles it well)
    fig = Figure(resolution=(820,820))
    ax  = Axis(fig[1,1], title=title, xlabel="A … M", ylabel="")
    sc  = scatter!(ax, xs, ys; markersize=6, color=vec(ZZ), colormap=:viridis)
    # draw triangle frame
    lines!(ax, [0, 0.5, 1, 0], [0, √3/2, 0, 0], color=:black)
    Colorbar(fig[1,2], sc, label="retained area")
    save(joinpath(outdir, "Fig1_surface_$(Symbol(geom)).png"), fig); fig
end

fig1a = plot_surface(:front;     title="Retained suitable area at f*=$(fstar) — front")
fig1b = plot_surface(:clustered; title="Retained suitable area at f*=$(fstar) — clustered")
fig1c = plot_surface(:random;    title="Retained suitable area at f*=$(fstar) — random")

# ---------- Fig 2: Geometry contrast (Random−Front, Clustered−Front) ----------
function plot_contrast(g1::Symbol, g0::Symbol; title="")
    xs, ys, _, Z1 = SFC.surf[g1]
    _,  _, _, Z0  = SFC.surf[g0]
    Δ = vec(Z1) .- vec(Z0)
    fig = Figure(resolution=(820,820))
    ax  = Axis(fig[1,1], title=title, xlabel="A … M", ylabel="")
    sc  = scatter!(ax, xs, ys; markersize=6, color=Δ, colormap=:balance)
    lines!(ax, [0, 0.5, 1, 0], [0, √3/2, 0, 0], color=:black)
    Colorbar(fig[1,2], sc, label="Δ retained area")
    save(joinpath(outdir, "Fig2_contrast_$(g1)_minus_$(g0).png"), fig); fig
end

fig2a = plot_contrast(:random,    :front;     title="Δ Random − Front at f*=$(fstar)")
fig2b = plot_contrast(:clustered, :front;     title="Δ Clustered − Front at f*=$(fstar)")

# ---------- Fig 3: Elasticity map ∂y/∂f at f* ----------
function elasticity_surface(geom::Symbol; df=0.02)
    xs = SFC.xs; ys = SFC.ys; pts = SFC.pts
    Z  = zeros(length(pts))
    for (t,(wA,wB,wM)) in enumerate(pts)
        f1, f2 = clamp(fstar-df, 0.0, 0.95), clamp(fstar+df, 0.0, 0.95)
        vals1 = Float64[]; vals2 = Float64[]
        for (k,pool) in enumerate(pools)
            push!(vals1, realized_area_mean(; rng, grid, pars, wA, wB, wM, f=f1, geom, pool, seed_A=300+k))
            push!(vals2, realized_area_mean(; rng, grid, pars, wA, wB, wM, f=f2, geom, pool, seed_A=300+k))
        end
        Z[t] = (mean(vals2) - mean(vals1)) / (f2 - f1)   # ∂y/∂f (negative)
    end
    fig = Figure(resolution=(820,820))
    ax  = Axis(fig[1,1], title="Elasticity ∂y/∂f at f*=$(fstar) — $(String(geom))", xlabel="A … M", ylabel="")
    sc  = scatter!(ax, xs, ys; markersize=6, color=Z, colormap=:plasma)
    lines!(ax, [0, 0.5, 1, 0], [0, √3/2, 0, 0], color=:black)
    Colorbar(fig[1,2], sc, label="∂y/∂f")
    save(joinpath(outdir, "Fig3_elasticity_$(Symbol(geom)).png"), fig); fig
end

fig3a = elasticity_surface(:front)
fig3b = elasticity_surface(:clustered)
fig3c = elasticity_surface(:random)

# ---------- Fig 4: Edge slices with CIs ----------
function edge_slices(; N=11)
    # three edges: A–M (B=0), A–B (M=0), B–M (A=0)
    edges = [(:AM, [(x, 0.0, 1-x) for x in range(0,1;length=N)]),
             (:AB, [(x, 1-x, 0.0) for x in range(0,1;length=N)]),
             (:BM, [(0.0, x, 1-x) for x in range(0,1;length=N)])]
    fig = Figure(resolution=(1200,350))
    for (col,(name, pts)) in enumerate(edges)
        ax = Axis(fig[1,col], title=String(name), xlabel="mixture coordinate", ylabel="retained area")
        for (gcol, geom) in zip([:black, :orange, :steelblue], geom_list)
            μ = Float64[]; lo = Float64[]; hi = Float64[]
            for (wA,wB,wM) in pts
                vals = Float64[]
                for (k,pool) in enumerate(pools)
                    push!(vals, realized_area_mean(; rng, grid, pars, wA, wB, wM, f=fstar, geom, pool, seed_A=500+k))
                end
                push!(μ, mean(vals)); qs = quantile(vals, [0.1, 0.9]); push!(lo, qs[1]); push!(hi, qs[2])
            end
            x = range(0,1; length=N)
            lines!(ax, x, μ; color=gcol, label=String(geom))
            band!(ax, x, lo, hi; color=(gcol, 0.2))
        end
        axislegend(ax, position=:lt)
    end
    save(joinpath(outdir, "Fig4_edges_CI.png"), fig); fig
end
fig4 = edge_slices()

# ---------- Fig 5: Tail-risk vs f for representative mixtures ----------
function tail_risk_panel()
    mixes = Dict(
        "A-heavy" => (0.7, 0.1, 0.2),
        "Mixed"   => (0.33, 0.34, 0.33),
        "B-heavy" => (0.1, 0.8, 0.1)
    )
    fig = Figure(resolution=(1200,350))
    for (i,(label,(wA,wB,wM))) in enumerate(pairs(mixes))
        ax = Axis(fig[1,i], title=label, xlabel="area lost f", ylabel="P10 per-species retained area")
        for (gcol, geom) in zip([:black, :orange, :steelblue], geom_list)
            P10 = Float64[]
            for f in loss_fracs
                vals = Float64[]
                for (k,pool) in enumerate(pools)
                    # per-species retained areas from P_real (reuse blend rule)
                    A = BSH.abiotic_matrix(pool, grid; seed=700+k)
                    keepfrac = 1 - f
                    keep = geom === :random    ? HL.random_mask(rng, grid.C, keepfrac) :
                           geom === :clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                                 HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)
                    M = pars.movement === :off ? trues(pool.S, grid.C) :
                        BSH.movement_gate(grid, A, keep; τA=pars.τA, T=pars.T)
                    P_AM, _ = BSH.assemble_AM(pool, grid, A, keep; pars=pars)
                    P_BAM   = BSH.assemble_BAM(pool, grid, A, keep; pars=pars).P
                    P_Aonly = falses(pool.S, grid.C)
                    @inbounds for s in 1:pool.S, i in 1:grid.C
                        P_Aonly[s,i] = keep[i] & (A[s,i] ≥ pars.τA) & M[s,i]
                    end
                    P_M = P_Aonly
                    P_real = falses(pool.S, grid.C)
                    @inbounds for s in 1:pool.S, i in 1:grid.C
                        val = wA*(P_Aonly[s,i] ? 1.0 : 0.0) + wM*(P_M[s,i] ? 1.0 : 0.0) + wB*(P_BAM[s,i] ? 1.0 : 0.0)
                        P_real[s,i] = (val ≥ 0.5)
                    end
                    cons = findall(!, pool.basal)
                    spvals = [sum(@view P_real[s, :]) / grid.C for s in cons]
                    push!(vals, quantile(spvals, 0.10))
                end
                push!(P10, mean(vals))
            end
            lines!(ax, collect(loss_fracs), P10; color=gcol, label=String(geom))
        end
        axislegend(ax, position=:rt)
    end
    save(joinpath(outdir, "Fig5_tailrisk.png"), fig); fig
end
fig5 = tail_risk_panel()

@info "All figs saved in $(outdir)"
# ============================================================================
