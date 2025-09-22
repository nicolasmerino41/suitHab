include("../SetUp.jl");
include("src/grids.jl");      using .Grids
include("src/metawebs.jl");   using .Metawebs
include("src/hl.jl");         using .HL
include("src/bsh.jl");        using .BSH
include("src/metrics.jl");    using .Metrics
include("src/figs.jl");       using .Figs
include("src/grids_init.jl"); using .GridsInit
include("src/fig2_metrics.jl"); using .Fig2Metrics
include("src/metaweb_sweep.jl");       using .MetawebSweep
mkpath("Paper/AttemptOfFinalResults1")

# ---------- shared setup (once) ----------
rng = MersenneTwister(123)
nx, ny = 120, 120
S = 200
basal_frac = 0.25

# Grids (single source of truth)
G = GridsInit.build_default_grids(nx, ny; seeds=(11,12,13,14))
grid_grad  = G.grad
grid_patch = G.patch
grid_mosa  = G.mosaic
grid_ridge = G.ridge

# Metaweb archetypes
pool_low  = Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:low)
pool_mid  = Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:mid)
pool_high = Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:high)

# Parameters (movement ON; tweak τB/T if you want crisper contrasts)
pars = BSH.BAMParams(; τA=0.5, τB=0.5, τocc=0.2, γ=3.0, movement=:component, T=8)

loss_fracs = 0.2:0.1:0.8
pool = pool_mid
grid = grid_grad     # pick the grid you want to headline; gradient is interpretable

# ---------- R1: Absolute curves (fast) ----------
curves = Metrics.abs_bsh_vs_loss(; rng, pool, grid, pars,
    loss_fracs=loss_fracs, seed_A=1)

# Lightweight plot (Makie, inside Figs for aesthetics if you prefer)
begin
    fig = Figure(; size=(1100,360))
    for (j, geom) in enumerate((:front, :clustered, :random))
        ax = Axis(fig[1,j], title=String(geom),
                xlabel="area lost (fraction)",
                ylabel="BSH (mean over consumers / original area)")
        lines!(ax, curves[geom].loss, curves[geom].AM; color=(:gray,0.6), linestyle=:dash, label="AM")
        lines!(ax, curves[geom].loss, curves[geom].BAM; color=:black, label="BAM")
        if j==1; axislegend(ax, position=:rt); end
    end
    display(fig)
end
save("Paper/AttemptOfFinalResults1/R1_abs_BSH_vs_loss.png", fig)

# ---------- R2: per-species CDFs at data-driven f* ----------
best = Metrics.ks_best_fstar(; rng, pool, grid, pars,
                             loss_fracs=loss_fracs, geometry=:front, seed_A=1)
fstar = best.f   # KS-maximizing f*
@info "Fig R2 using f*=$(round(fstar,digits=2)) (KS=$(round(best.KS,digits=3)))"

function cdf_xy(v)
    s = sort(v); n = length(s)
    (s, (1:n)./n)
end

cdfs = Dict{Symbol,Any}()
for g in (:random,:clustered,:front)
    rAM, rBAM = BSH.per_species_relative_loss(rng, pool, grid, pars; fstar, geometry=g, seed_A=1)
    x1,y1 = cdf_xy(clamp.(-rAM, 0,1))
    x2,y2 = cdf_xy(clamp.(-rBAM, 0,1))
    cdfs[g] = (; xAM=x1, yAM=y1, xBAM=x2, yBAM=y2)
end

fig2 = Figs.fig2_cdfs(cdfs; fstar, title="Per-species losses at f*")
save("Paper/AttemptOfFinalResults1/R2_cdfs.png", fig2)

# Metrics panel (KS, tail, Gini) — single PNG
Fig2Metrics.fig2_metrics_block(; rng, pool, grid, pars,
    fstar=fstar, tail_cut=0.80, seed_A=1,
    out="Paper/AttemptOfFinalResults1/R2_metrics_block.png")

    # ---------- R3: filled ternary Δ(Random−Front) at f* ----------
# Define your mixtures once (reuse from earlier)
weights = [(1.0,0.0,0.0), (0.75,0.25,0.0), (0.5,0.5,0.0), (0.25,0.75,0.0),
           (0.0,1.0,0.0), (0.0,0.75,0.25), (0.0,0.5,0.5), (0.0,0.25,0.75),
           (0.0,0.0,1.0), (0.25,0.0,0.75), (0.5,0.0,0.5), (0.75,0.0,0.25),
           (0.33,0.34,0.33)]  # center

boot = Metrics.bootstrap_mixture_maps(; rng=MersenneTwister(11),
            pool, grid, pars, weights, fstar,
            A_fn=BSH.abiotic_matrix, agg=:mean, kreq=1,
            A_seeds=1:16)

maxabs = maximum(abs, boot.dRF_mean)
filled_ternary_map(
    weights, boot.dRF_mean;
    title="Δ (Random − Front) at f*=$(round(fstar,digits=2))",
    colormap=:balance, clim=(-maxabs, maxabs),
    res=24,
    fname="Paper/AttemptOfFinalResults1/R3_delta_RF_filled.png"
)

    # ---------- R4: DiD sweep (fast; single grid) ----------
mkpath("Paper/AttemptOfFinalResults1/sweep")

sweep = MetawebSweep.run_metaweb_sweep(
    ; grids=[("gradient",grid)],
      S=S, basal_frac=basal_frac, archetype=:mid,
      loss_fracs=loss_fracs, fstar=fstar,
      keep_probs=[1.0, 0.7, 0.5, 0.3],
      caps=[9999, 6, 4, 2],
      pars=pars,
      outdir="Paper/AttemptOfFinalResults1/sweep"
)
# Output image: Paper/AttemptOfFinalResults1/sweep/metaweb_sweep_gradient_f$(Int(round(100*fstar))).png
# ---------- R5 (optional): A-only realities ----------
realA = MetawebSweep.run_realities_Aonly(
    ; rng, grid, loss_fracs=loss_fracs,
      S=S, basal_frac=basal_frac, seed_A=1,
      Npools=10, geoms=(:front,:clustered,:random))

begin
    fig = Figure(; size=(1100,360), fontsize=13)
    for (j, g) in enumerate((:front,:clustered,:random))
        ax = Axis(fig[1,j], title=String(g),
                xlabel="area lost (fraction)",
                ylabel="suitable area (mean over consumers / original area)")
        lines!(ax, realA[g].loss, realA[g].ABM; label="ABM")
        lines!(ax, realA[g].loss, realA[g].MAB; label="MAB")
        lines!(ax, realA[g].loss, realA[g].BAM; label="BAM")
        if j==1; axislegend(ax, position=:lb); end
    end
    display(fig)
end
save("Paper/AttemptOfFinalResults1/R5_realities_Aonly.png", fig)
