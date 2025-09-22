include("../SetUp.jl")
include("src/grids.jl");      using .Grids
include("src/grids_init.jl"); using .GridsInit
include("src/metawebs.jl");   using .Metawebs
include("src/hl.jl");         using .HL
include("src/bsh.jl");        using .BSH
include("src/metrics.jl");    using .Metrics
include("src/ternary_plots.jl"); using .TernaryPlots

# ---------------- Settings ----------------
rng = MersenneTwister(123)
nx, ny = 120, 120
S = 200
basal_frac = 0.25

# Grids
G = GridsInit.build_default_grids(nx, ny; seeds=(11,12,13,14))
grid_grad = G.grad   # use gradient for clarity of front effects

# Pool (fixed across scenarios)
pool = Metawebs.build_metaweb_archetype(rng; S=S, basal_frac=basal_frac, archetype=:mid)

# Loss levels
fstars = [0.4, 0.6, 0.8]

# Mixture lattice on simplex (step = 0.25)
weights = [(a,b,1-a-b) for a in 0:0.25:1 for b in 0:0.25:(1-a)]

# ---------- Scenario configs ----------
# (a) Baseline
pars_base = BSH.BAMParams(; τA=0.5, τB=0.50, movement=:component, T=8, τocc=0.2, γ=3.0)
A_base   = (pool, grid; seed) -> BSH.abiotic_matrix(pool, grid; niche_width=0.12, seed=seed)
agg_base, kreq_base = :mean, 1

# (b) B-fragile (structured & stringent biotics)
pars_Bfrag = BSH.BAMParams(; τA=0.5, τB=0.50, movement=:component, T=8)
A_Bfrag = (pool, grid; seed) -> BSH.abiotic_matrix_aligned(
    pool, grid; seed=seed, bias_basal=0.80, align=0.85, niche_basal=0.08, niche_cons=0.12)
agg_Bfrag, kreq_Bfrag = :kofn, 2   # require ≥2 prey

# (c) M-strong
pars_Mstrong = BSH.BAMParams(; τA=0.5, τB=0.50, movement=:component, T=14)
A_Mstrong = A_base
agg_M, kreq_M = :mean, 1

# (d) A-anisotropic (front should be worse under A)
# (we reuse grid_grad which already has a strong axis; same A as baseline)
pars_Aaniso = pars_base
A_Aaniso = (pool, grid; seed) -> BSH.abiotic_matrix(pool, grid; niche_width=0.10, seed=seed)
agg_A, kreq_A = :mean, 1

scenarios = [
    ("baseline",  pars_base,   A_base,   agg_base,   kreq_base),
    ("Bfragile",  pars_Bfrag,  A_Bfrag,  agg_Bfrag,  kreq_Bfrag),
    ("Mstrong",   pars_Mstrong,A_Mstrong,agg_M,      kreq_M),
    ("Aaniso",    pars_Aaniso, A_Aaniso, agg_A,      kreq_A),
]

# ------------ Runners -------------
function compute_maps(; name, pars, A_fn, agg, kreq, fstar, grid=grid_grad)
    geomlist = (:random, :clustered, :front)
    # Elasticities per geometry
    elast = Dict{Symbol,Vector{Float64}}()
    for g in geomlist
        elast[g] = [Metrics.mixture_elasticity_at_f(; rng, pool, grid, pars,
                        wA=a, wB=b, wM=m, fstar=fstar, geometry=g,
                        A_fn=A_fn, agg=agg, kreq=kreq, seed_A=1)
                    for (a,b,m) in weights]
    end
    # Geometry contrasts
    dRF = [Metrics.mixture_geom_delta(; rng, pool, grid, pars,
                wA=a, wB=b, wM=m, fstar=fstar, other=:front,
                A_fn=A_fn, agg=agg, kreq=kreq, seed_A=1)
           for (a,b,m) in weights]
    dRC = [Metrics.mixture_geom_delta(; rng, pool, grid, pars,
                wA=a, wB=b, wM=m, fstar=fstar, other=:clustered,
                A_fn=A_fn, agg=agg, kreq=kreq, seed_A=1)
           for (a,b,m) in weights]
    (; elast, dRF, dRC)
end

function plot_ternary_set(; title, z_by_geom::Dict{Symbol,Vector{Float64}}, fname)
    fig = Figure(; size=(1500,420))
    geoms = (:random, :clustered, :front)
    for (i,g) in enumerate(geoms)
        ax = Axis(fig[1,i], title=string(g), xticksvisible=false, yticksvisible=false)
        s = TernaryPlots.ternary_scatter!(ax, weights, z_by_geom[g]; colormap=:plasma, markersize=16)
        Colorbar(fig[1,i+1], s, label="∂y/∂f at f*"; width=12)
    end
    fig[0,1] = Label(fig, title, fontsize=20, tellwidth=false)
    save(fname, fig)
    display(fig)
end

function plot_ternary_delta(; title, vals::Vector{Float64}, fname)
    # single ternary with diverging map
    fig = Figure(; size=(520,520))
    ax = Axis(fig[1,1], title=title, xticksvisible=false, yticksvisible=false)
    s = TernaryPlots.ternary_scatter!(ax, weights, vals; colormap=:balance, markersize=16)
    Colorbar(fig[1,2], s, label="Δ (Random − Other)")
    save(fname, fig); display(fig)
end

mkpath("Paper/figs/ternary")

# ------------ MAIN -------------
for (sname, pars, A_fn, agg, kreq) in scenarios
    for fstar in fstars
        res = compute_maps(; name=sname, pars, A_fn, agg, kreq, fstar)
        # Elasticities per geometry
        plot_ternary_set(
            ; title="Elasticity at f*=$(fstar), scenario=$(sname)",
              z_by_geom=res.elast,
              fname="Paper/figs/ternary/elasticity_$(sname)_f$(Int(round(100*fstar))).png")

        # Geometry contrasts (Random vs Front / Clustered)
        plot_ternary_delta(
            ; title="Random − Front at f*=$(fstar)  ($(sname))",
              vals=res.dRF,
              fname="Paper/figs/ternary/delta_RF_$(sname)_f$(Int(round(100*fstar))).png")

        plot_ternary_delta(
            ; title="Random − Clustered at f*=$(fstar)  ($(sname))",
              vals=res.dRC,
              fname="Paper/figs/ternary/delta_RC_$(sname)_f$(Int(round(100*fstar))).png")

        # Effect-size summaries
        open("Paper/figs/ternary/effects_summary.txt", "a") do io
            mRF = mean(res.dRF); rRF = quantile(res.dRF, [0.1,0.9])
            mRC = mean(res.dRC); rRC = quantile(res.dRC, [0.1,0.9])
            println(io,
                "scenario=$(sname) f*=$(fstar):  " *
                "meanΔ(R−F)=$(round(mRF; digits=4))  [p10=$(round(rRF[1]; digits=4)), p90=$(round(rRF[2]; digits=4))];  " *
                "meanΔ(R−C)=$(round(mRC; digits=4))  [$(round(rRC[1]; digits=4)), $(round(rRC[2]; digits=4))]"
            )
        end
    end
end

# -----------------------------------
# ------------ BOOTSTRAP ------------
# -----------------------------------
# draw ternary with mean color and CI width mapped to marker edge alpha
function plot_ternary_set_CI(; title, mean_by_geom::Dict{Symbol,Vector{Float64}},
        lo_by_geom::Dict{Symbol,Vector{Float64}},
        hi_by_geom::Dict{Symbol,Vector{Float64}}, fname)

    fig = Figure(; size=(1600,480))
    geoms = (:random, :clustered, :front)

    for (i,g) in enumerate(geoms)
        ax = Axis(fig[1,i]; title=string(g), xticksvisible=false, yticksvisible=false)
        ax.aspect = DataAspect()

        # CI width -> alpha
        w = hi_by_geom[g] .- lo_by_geom[g]
        w_norm = (w .- minimum(w)) ./ (maximum(w) - minimum(w) + 1e-12)
        alphas = .35 .+ .6 .* (1 .- w_norm)
        stroke_colors = RGBAf.(0, 0, 0, alphas)

        xs = Float64[]; ys = Float64[]
        for (a,b,m) in weights
            c = TernaryPlots.ternary_coords(a,b,m)
            push!(xs,c.x); push!(ys,c.y)
        end

        sc = scatter!(ax, xs, ys;
            color=mean_by_geom[g], colormap=:plasma,
            strokewidth=1.2, strokecolor=stroke_colors,
            markersize=18)

        Colorbar(fig[1,i+1], sc, label="∂y/∂f at f*"; width=12)
        xlims!(ax, -0.05, 1.05); ylims!(ax, -0.05, √3/2+0.05)
    end

    fig[0,1] = Label(fig, title, fontsize=20, tellwidth=false)
    save(fname, fig); display(fig)
end


function plot_ternary_delta_CI(; title, μ::Vector{Float64}, lo::Vector{Float64}, hi::Vector{Float64}, fname)
    fig = Figure(; size=(620,560))
    ax = Axis(fig[1,1]; title=title, xticksvisible=false, yticksvisible=false)
    ax.aspect = DataAspect()

    xs = Float64[]; ys = Float64[]
    for (a,b,m) in weights
        c = TernaryPlots.ternary_coords(a,b,m)
        push!(xs,c.x); push!(ys,c.y)
    end

    # Uncertainty -> alpha
    w = hi .- lo
    w_norm = (w .- minimum(w)) ./ (maximum(w)-minimum(w)+1e-12)
    alphas = .35 .+ .6 .* (1 .- w_norm)
    stroke_colors = RGBAf.(0, 0, 0, alphas)

    sc = scatter!(ax, xs, ys;
        color=μ, colormap=:balance,
        markersize=18, strokewidth=1.2,
        strokecolor=stroke_colors)

    Colorbar(fig[1,2], sc, label="Δ (Random − Other)")
    xlims!(ax, -0.05, 1.05); ylims!(ax, -0.05, √3/2+0.05)
    save(fname, fig); display(fig)
end

# --- bootstrap compute + plots for each scenario & f* ---
A_seeds = 1:24                # << number of bootstrap draws (increase if you want)
for (sname, pars, A_fn, agg, kreq) in scenarios
    for fstar in fstars
        boot = Metrics.bootstrap_mixture_maps(; rng, pool, grid=grid_grad, pars,
                    weights, fstar, A_fn, agg, kreq, A_seeds=A_seeds)

        # Elasticities (mean + CI) per geometry
        plot_ternary_set_CI(
            ; title="Elasticity at f*=$(fstar), scenario=$(sname)",
              mean_by_geom=boot.elastic_mean,
              lo_by_geom=boot.elastic_lo,
              hi_by_geom=boot.elastic_hi,
              fname="Paper/figs/ternary/elasticity_$(sname)_f$(Int(round(100fstar)))_boot.png")

        # Δ maps with CI
        plot_ternary_delta_CI(
            ; title="Random − Front at f*=$(fstar)  ($(sname))",
              μ=boot.dRF_mean, lo=boot.dRF_lo, hi=boot.dRF_hi,
              fname="Paper/figs/ternary/delta_RF_$(sname)_f$(Int(round(100fstar)))_boot.png")

        plot_ternary_delta_CI(
            ; title="Random − Clustered at f*=$(fstar)  ($(sname))",
              μ=boot.dRC_mean, lo=boot.dRC_lo, hi=boot.dRC_hi,
              fname="Paper/figs/ternary/delta_RC_$(sname)_f$(Int(round(100fstar)))_boot.png")

        # Text summaries
        open("Paper/figs/ternary/effects_summary_boot.txt","a") do io
            mRF = mean(boot.dRF_mean); wRF = mean(boot.dRF_hi .- boot.dRF_lo)
            mRC = mean(boot.dRC_mean); wRC = mean(boot.dRC_hi .- boot.dRC_lo)
            println(io, "scenario=$(sname) f*=$(fstar):  meanΔ(R−F)=$(round(mRF,digits=4))  mean-CIwidth=$(round(wRF,digits=4)); ",
                        "meanΔ(R−C)=$(round(mRC,digits=4))  mean-CIwidth=$(round(wRC,digits=4))")
        end
    end
end

