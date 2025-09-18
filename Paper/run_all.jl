# loading modules
include("../SetUp.jl");
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
include("src/plot_grids.jl")
include("src/metaweb_variety_sweep.jl"); using .MetawebVarietySweep
include("src/fig2_metrics.jl"); using .Fig2Metrics
mkpath("Paper/figs")

# ----------------------------
# 0) Settings
# ----------------------------
rng = MersenneTwister(123)
nx, ny = 120, 120
S = 200
basal_frac = 0.25

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

########################### REMOVE ME LATER ############################
# movement ON to make geometry relevant
pars = BSH.BAMParams(; τA=0.6, τB=0.6, τocc=0.2, γ=3.0, movement=:component, T=10)

# basal squeezed to low climate + consumers aligned to prey
A_fn_punch = (pool, grid; seed=1) -> BSH.abiotic_matrix_aligned(
    pool, grid; seed=seed, niche_basal=0.08, niche_cons=0.12, bias_basal=0.8, align=0.85)

# require at least 2 prey in a cell (low redundancy hurts)
relcurves = BSH.relative_loss_curves(rng, pool_low, grid_grad, pars;
    loss_fracs=0.2:0.1:0.8, seed_A=1, A_fn=A_fn_punch, agg=:kofn, kreq=2)

pfail     = BSH.pfail_curve(; rng, pool=pool_low, grid=grid_grad, pars,
    loss_fracs=0.2:0.1:0.8, seed_A=1, A_fn=A_fn_punch, agg=:kofn, kreq=2)
#################################################################################

# Parameters
# movement OFF
pars_off = BSH.BAMParams(; τA=0.5, τB=0.55, τocc=0.2, γ=3.0, movement=:off, T=6)
# movement ON (connected component gate with size T)
pars_on  = BSH.BAMParams(; τA=0.5, τB=0.5, τocc=0.2, γ=3.0, movement=:component, T=8)
# STRONG COMBINATION
pars_strong = BSH.BAMParams(; τA=0.55, τB=0.50, movement=:component, T=8)

pars = pars_on # TODO CHOOSE DEPENDING ON YOUR NEEDS

loss_fracs = 0.2:0.05:0.8
fstar      = 0.6
at_index   = findfirst(==(fstar), collect(loss_fracs))

# ------------------------------------
# 1) Fig.1 — core (grid=gradient, metaweb=mid)
# ------------------------------------
pool = pool_low
grid = grid_grad

relcurves = BSH.relative_loss_curves(rng, pool, grid, pars; loss_fracs, seed_A=1)
rel_strong = BSH.relative_loss_curves(rng, pool_low, grid_ridge, pars_strong; loss_fracs, seed_A=2)
placebo   = BSH.placebo_curves(; rng, pool, grid, pars, loss_fracs, seed_A=1)
pfail     = BSH.pfail_curve(; rng, pool, grid, pars, loss_fracs, seed_A=1)

fig1 = Figs.fig1_core(
    relcurves, pfail;
    placebo = placebo,
    with_placebo = false,
    title = "AM vs BAM (baseline equalized) — " *
            (pool === pool_high ? "high" : pool === pool_mid ? "mid" : "low") * " redundancy, " *
            (grid === grid_grad ? "gradient" : grid === grid_patch ? "patchy" :
             grid === grid_mosa ? "mosaic" : grid === grid_ridge ? "ridge" : "unknown") * " grid, and movement " *
            (pars.movement === :off ? "OFF" : "ON")
)

save("Paper/figs/Fig1_core.png", fig1); display(fig1)

# ------------------------------------
# 1.1) Fig.1.1 — FigsSweeps.sweep_fig1_3x3
# ------------------------------------
# choose which grid shows contrasts best (ridge often pops); swap to grid_grad/patch/mosa as desired
grid_for_sweep = "gradient"

# 3×3 per archetype (movement ON). Change tauB_list/T_list inside if you want.
FigSweeps.sweep_fig1_3x3(
    ; rng, archetype=:low, grids_vec=grids_vec, grid_type=grid_for_sweep,
    with_placebo=false,
    outpath="Paper/figs/Fig1Sweeps/Fig1_sweep_$(grid_for_sweep)_low.png"
)

FigSweeps.sweep_fig1_3x3(
    ; rng, archetype=:mid, grids_vec=grids_vec, grid_type=grid_for_sweep,
    with_placebo=false,
    outpath="Paper/figs/Fig1Sweeps/Fig1_sweep_$(grid_for_sweep)_mid.png"
)

FigSweeps.sweep_fig1_3x3(
    ; rng, archetype=:high, grids_vec=grids_vec, grid_type=grid_for_sweep,
    with_placebo=false,
    outpath="Paper/figs/Fig1Sweeps/Fig1_sweep_$(grid_for_sweep)_high.png"
)

# ------------------------------------
# 1.2) Fig.1 + P_fail — 3×3 per archetype, movement OFF and ON
# ------------------------------------
grid_for_sweep = "gradient"   # or "gradient" / "patchy" / "ridge" / "mosaic"

# movement OFF
FigSweeps.sweep_fig1_pair_3x3(
    ; rng, archetype=:low,  grids_vec=grids_vec, grid_type=grid_for_sweep,
      movement=:off, with_placebo=false,
      outpath="Paper/figs/Fig1Sweeps/Fig1PAIR_$(grid_for_sweep)_low_Moff.png")
FigSweeps.sweep_fig1_pair_3x3(
    ; rng, archetype=:mid,  grids_vec=grids_vec, grid_type=grid_for_sweep,
      movement=:off, with_placebo=false,
      outpath="Paper/figs/Fig1Sweeps/Fig1PAIR_$(grid_for_sweep)_mid_Moff.png")
FigSweeps.sweep_fig1_pair_3x3(
    ; rng, archetype=:high, grids_vec=grids_vec, grid_type=grid_for_sweep,
      movement=:off, with_placebo=false,
      outpath="Paper/figs/Fig1Sweeps/Fig1PAIR_$(grid_for_sweep)_high_Moff.png")

# movement ON (component gate)
FigSweeps.sweep_fig1_pair_3x3(
    ; rng, archetype=:low,  grids_vec=grids_vec, grid_type=grid_for_sweep,
      movement=:component, with_placebo=false,
      outpath="Paper/figs/Fig1Sweeps/Fig1PAIR_$(grid_for_sweep)_low_Mon.png")
FigSweeps.sweep_fig1_pair_3x3(
    ; rng, archetype=:mid,  grids_vec=grids_vec, grid_type=grid_for_sweep,
      movement=:component, with_placebo=false,
      outpath="Paper/figs/Fig1Sweeps/Fig1PAIR_$(grid_for_sweep)_mid_Mon.png")
FigSweeps.sweep_fig1_pair_3x3(
    ; rng, archetype=:high, grids_vec=grids_vec, grid_type=grid_for_sweep,
      movement=:component, with_placebo=false,
      outpath="Paper/figs/Fig1Sweeps/Fig1PAIR_$(grid_for_sweep)_high_Mon.png")

# ------------------------------------
# 2) Fig.2 — per-species CDFs at f*
# ------------------------------------
grid = grid_grad
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
fig2 = Figs.fig2_cdfs(cdfs; fstar, title="Distributional impact (per-species)_grid $(grid === grid_mosa ? "mosaic" : grid === grid_ridge ? "ridge" : grid === grid_grad ? "gradient" : "patchy")")
save("Paper/figs/Fig2_cdfs.png", fig2); display(fig2)

# ------------------------------------
# 2.1) Fig.2 — per-species CDFs WITH NUMBERS
# ------------------------------------
tailα = 0.80
metrix = Dict{Symbol,NamedTuple}()
fig2b = Figs.fig2_metric_numbers(metrix; fstar, tail=tailα)
save("figs/Fig2_metrics_grid_$(nameof(grid)).png", fig2b); display(fig2b)

# ------------------------------------
# 2.2) Fig.2 — per-species CDFs at f*, WITH METRICS
# ------------------------------------
# A) KS / tailΔ / GiniΔ table across geometries
Fig2Metrics.fig2_metrics_block(; rng, pool, grid, pars,
    fstar=fstar, tail_cut=0.80, seed_A=1,
    out="Paper/figs/fig2_metrics/Fig2_metrics.png")

# B) Optional: CDF of per-species biotic failure at f* (pick geometry)
Fig2Metrics.per_species_failure_cdf(; rng, pool, grid, pars,
    fstar=fstar, geometry=:front, seed_A=1,
    out="Paper/figs/Fig2_failure_front.png")

# C) Optional: 'who pays?' — penalty vs diet size (consumers only)
Fig2Metrics.penalty_vs_diet(; rng, pool, grid, pars,
    fstar=fstar, geometry=:front, seed_A=1,
    out="Paper/figs/penalty_vs_diet_front.png")

# ------------------------------------
# 3) Fig.3 — rank small multiples across archetypes × grids at f*
# ------------------------------------
combos = [("High-R, gradient", pool_high, grid_grad),
          ("Low-R, gradient",  pool_low,  grid_grad),
          ("Mid-R, patchy",    pool_mid,  grid_patch),
          ("Mid-R, ridge",     pool_mid,  grid_ridge)]
combos = [
    ("High-R, gradient", pool_high, grid_grad), ("Low-R, gradient",  pool_low,  grid_grad),
    ("Mid-R, patchy",    pool_mid,  grid_patch), ("Mid-R, ridge",     pool_mid,  grid_ridge),
    ("High-R, patchy",   pool_high, grid_patch), ("Low-R, ridge",     pool_low,  grid_ridge),
    ("High-R, mosaic",   pool_high, grid_mosa), ("Low-R, mosaic",    pool_low,  grid_mosa),
    ("Mid-R, gradient",  pool_mid,  grid_grad), ("High-R, ridge",    pool_high, grid_ridge),
    ("Low-R, patchy",    pool_low,  grid_patch), ("Mid-R, mosaic",    pool_mid,  grid_mosa),
]

ranks = NamedTuple[]
for (label, p, g) in combos
    rel = BSH.relative_loss_curves(rng, p, g, pars; loss_fracs)
    am, bm, flip = Metrics.rank_flip(rel; at_index)
    push!(ranks, (; label, am, bam=bm))
end

fig3 = Figs.fig3_rank_smallmultiples(ranks; title="Worst geometry at f* = $(fstar)", ncols=4)
save("Paper/figs/Fig3_rank.png", fig3); display(fig3)

# ------------------------------------
# 4) Fig.4 — regime maps over (D, R)
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

# === grid maps (already fixed by you) =========================================include("src/plot_grids.jl")
plot_all_grids(; nx=nx, ny=ny, out="Paper/figs/Grids.png")

# === metaweb sweep & descriptive =============================================
sweep_res60 = MetawebSweep.run_metaweb_sweep(
    ; grids=grids_vec,                                 # <- uses the pairs from GridsInit
    S=S, basal_frac=basal_frac, archetype=:mid,
    loss_fracs=loss_fracs, fstar=fstar,
    keep_probs=[1.0, 0.7, 0.5, 0.3], caps=[9999, 6, 4, 2],
    pars=pars,
    outdir="Paper/figs/sweep_f$(Int(100*fstar))"
)

# Descriptives for archetypes + some of the new metawebs (optional)
variety_pools = [
    pool_low, pool_mid, pool_high,
    Metawebs.build_metaweb_niche(rng; S=S, beta=1.8),
    Metawebs.build_metaweb_modular(rng; S=S, K=3, p_in=0.30, p_out=0.04),
    Metawebs.build_metaweb_powerlaw(rng; S=S, alpha=2.4, kmax=8)
]

variety_names = ["low","mid","high","niche_b1.8","mod_K3","pow_a2.4"]
MetawebDescriptive.plot_metaweb_spectrum(variety_pools; names=variety_names,
                                         outdir="Paper/figs/metawebs/")

# === Variety sweeps across more metawebs =====================================
# choose which grid you want to showcase; swap "patchy"/"mosaic"/"ridge"/"gradient"
MetawebVarietySweep.run_metaweb_variety_sweeps(
    ; rng, grids_vec=grids_vec, grid_type="patchy",
      S=S, tauB_list=[0.40,0.50,0.60], T_list=[6,8,12],
      movement=:component,    # try :off as well to compare
      τA=0.5, τocc=0.2, γ=3.0, loss_fracs=loss_fracs, seed_A=1,
      outdir="Paper/figs/metaweb_variety",
      display_figs=true
)
