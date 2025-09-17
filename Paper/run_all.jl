# loading modules
include("src/grids.jl");      using .Grids
include("src/metawebs.jl");   using .Metawebs
include("src/hl.jl");         using .HL
include("src/bsh.jl");        using .BSH
include("src/metrics.jl");    using .Metrics
include("src/figs.jl");       using .Figs
include("src/grids_init.jl"); using .GridsInit
include("src/metaweb_sweep.jl");       using .MetawebSweep
include("src/metaweb_descriptive.jl"); using .MetawebDescriptive
include("src/figs_sweep.jl"); using .FigSweeps
mkpath("Paper/figs")

# ----------------------------
# 1) Settings (keep simple)
# ----------------------------
rng = MersenneTwister(123)
nx, ny = 60, 60
S = 175
basal_frac = 0.35

# --- BUILD GRIDS HERE (single source of truth) -------------------------------
G = GridsInit.build_default_grids(nx, ny; seeds=(11,12,13,14))
GridsInit.print_grid_diagnostics(G)

# preserve your existing variable names
grid_grad  = G.grad
grid_patch = G.patch
grid_mosa  = G.mosaic
grid_ridge = G.ridge
grids_vec  = GridsInit.as_pairs(G)   # for sweeps
# ---------------------------------------------------------------------------

# metaweb archetypes (accepted as-is; diagnostics printed)
pool_high = Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:high)
pool_mid  = Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:mid)
pool_low  = Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:low)

for (name,p) in [("high",pool_high),("mid",pool_mid),("low",pool_low)]
    d = Metawebs.metaweb_diagnostics(p)
    @info "metaweb $name  C=$(round(d.C,digits=3))  med_diet=$(d.meddiet)  basal=$(round(d.basal_share,digits=2))"
end

# Parameters
# movement OFF
pars_off = BSH.BAMParams(; τA=0.5, τB=0.55, τocc=0.2, γ=3.0, movement=:off, T=6)
# movement ON (connected component gate with size T)
pars_on  = BSH.BAMParams(; τA=0.5, τB=0.55, τocc=0.2, γ=3.0, movement=:component, T=6)
# STRONG COMBINATION
pars_strong = BSH.BAMParams(; τA=0.55, τB=0.50, movement=:component, T=8)

pars = pars_strong # TODO CHOOSE DEPENDING ON YOUR NEEDS

loss_fracs = 0.2:0.1:0.8
fstar      = 0.6
at_index   = findfirst(==(fstar), collect(loss_fracs))

# ------------------------------------
# 2) Fig.1 — core (grid=gradient, metaweb=mid)
# ------------------------------------
pool = pool_low
grid = grid_ridge

relcurves = BSH.relative_loss_curves(rng, pool, grid, pars; loss_fracs, seed_A=1)
rel_strong = BSH.relative_loss_curves(rng, pool_low, grid_ridge, pars_strong; loss_fracs, seed_A=2)
placebo   = BSH.placebo_curves(; rng, pool, grid, pars, loss_fracs, seed_A=1)
pfail     = BSH.pfail_curve(; rng, pool, grid, pars, loss_fracs, seed_A=1)

fig1 = Figs.fig1_core(relcurves, placebo, pfail;
    title = "AM vs BAM (baseline equalized) — " *
            (pool === pool_high ? "high" : pool === pool_mid ? "mid" : "low") * " redundancy, " *
            (grid === grid_grad ? "gradient" : grid === grid_patch ? "patchy" :
             grid === grid_mosa ? "mosaic" : grid === grid_ridge ? "ridge" : "unknown") * " grid, and movement " *
            (pars.movement === :off ? "OFF" : "ON")

)

save("Paper/figs/Fig1_core.png", fig1); display(fig1)

# ------------------------------------
# 2.1) Fig.1.1 — FigsSweeps.sweep_fig1_3x3
# ------------------------------------
# choose which grid shows contrasts best (ridge often pops); swap to grid_grad/patch/mosa as desired
grid_for_sweep = "ridge"

# 3×3 per archetype (movement ON). Change tauB_list/T_list inside if you want.
FigSweeps.sweep_fig1_3x3(; rng, archetype=:low, grids_vec=grids_vec, grid_type=grid_for_sweep,
    outpath="Paper/figs/Fig1_sweep_low_$(grid_for_sweep).png")

FigSweeps.sweep_fig1_3x3(; rng, archetype=:mid, grids_vec=grids_vec, grid_type=grid_for_sweep,
    outpath="Paper/figs/Fig1_sweep_mid_$(grid_for_sweep).png")

FigSweeps.sweep_fig1_3x3(; rng, archetype=:high, grids_vec=grids_vec, grid_type=grid_for_sweep,
    outpath="Paper/figs/Fig1_sweep_high_$(grid_for_sweep).png")

# ------------------------------------
# 3) Fig.2 — per-species CDFs at f*
# ------------------------------------
function cdf_xy(v)
    s = sort(v); n = length(s)
    (s, (1:n)./n)
end
cdfs = Dict{Symbol,Any}()
for g in (:random,:clustered,:front)
    rAM, rBAM = BSH.per_species_relative_loss(rng, pool, grid, pars; fstar, geometry=g)
    x1,y1 = cdf_xy(-rAM)   # use -Δ/BSH0 as positive losses
    x2,y2 = cdf_xy(-rBAM)
    cdfs[g] = (; xAM=x1, yAM=y1, xBAM=x2, yBAM=y2)
end
fig2 = Figs.fig2_cdfs(cdfs; fstar, title="Distributional impact (per-species)")
save("Paper/figs/Fig2_cdfs.png", fig2); display(fig2)

# ------------------------------------
# 4) Fig.3 — rank small multiples across archetypes × grids at f*
# ------------------------------------
combos = [("High-R, gradient", pool_high, grid_grad),
          ("Low-R, gradient",  pool_low,  grid_grad),
          ("Mid-R, patchy",    pool_mid,  grid_patch),
          ("Mid-R, ridge",     pool_mid,  grid_ridge)]

ranks = NamedTuple[]
for (label, p, g) in combos
    rel = BSH.relative_loss_curves(rng, p, g, pars; loss_fracs)
    am, bm, flip = Metrics.rank_flip(rel; at_index)
    push!(ranks, (; label, am, bam=bm))
end
fig3 = Figs.fig3_rank_smallmultiples(ranks; title="Worst geometry at f* = $(fstar)")
save("Paper/figs/Fig3_rank.png", fig3); display(fig3)

# ------------------------------------
# 5) Fig.4 — regime maps over (D, R)
# ------------------------------------
Dvals = [Grids.climate_tail_index(g) for g in (grid_grad, grid_patch, grid_mosa, grid_ridge)]
Rvals = Float64[]
for p in (pool_low, pool_mid, pool_high)
    diets = [length(p.prey[s]) for s in 1:p.S if !p.basal[s]]
    push!(Rvals, isempty(diets) ? 0.0 : quantile(diets, 0.95))
end

Zdid  = zeros(length(Rvals), length(Dvals))
Zflip = zeros(length(Rvals), length(Dvals))
for (j,(g, D)) in enumerate(zip((grid_grad,grid_patch,grid_mosa,grid_ridge), Dvals))
    for (i,(p, R)) in enumerate(zip((pool_low,pool_mid,pool_high), Rvals))
        rel = BSH.relative_loss_curves(rng, p, g, pars; loss_fracs)
        sumry = Metrics.regime_summary(rel; at_index)
        Zdid[i,j]  = sumry.meanDiD
        _,_,fl = Metrics.rank_flip(rel; at_index)
        Zflip[i,j] = fl ? 1.0 : 0.0
    end
end
fig4 = Figs.fig4_regimemap(Dvals, Rvals, Zdid, Zflip; title="Regime map at f* = $(fstar)")
save("Paper/figs/Fig4_regimes.png", fig4); display(fig4)

@info "Done. Figures saved in ./figs"

# === grid maps (already fixed by you) =========================================
include("src/plot_grids.jl")
plot_all_grids(; nx=nx, ny=ny, out="out/grids.png")

# === metaweb sweep & descriptive =============================================
sweep_res60 = MetawebSweep.run_metaweb_sweep(
    ; grids=grids_vec,                                 # <- uses the pairs from GridsInit
    S=S, basal_frac=basal_frac, archetype=:mid,
    loss_fracs=loss_fracs, fstar=fstar,
    keep_probs=[1.0, 0.7, 0.5, 0.3], caps=[9999, 6, 4, 2],
    pars=pars,
    outdir="Paper/figs/sweep_f$(Int(100*fstar))"
)

MetawebDescriptive.plot_metaweb_spectrum(
    [pool_low, pool_mid, pool_high];
    names=["low","mid","high"],
    outdir="Paper/figs/metawebs"
)

