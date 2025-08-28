include("Functions.jl")
include("Knee.jl")
## MAKE THE GRID AND SPECIES POOL ##
grid = make_grid(100, 100; seed=11)

pool_seeds = 1:5
mask_seeds = 1:30
metric = :fraction
sum_rand = sweep_ensemble(
    pool_seeds, mask_seeds; 
    kind=:random, grid=grid, 
    τ=0.55, S=200, 
    basal_frac=0.35, nseeds_cluster=1, metric=metric,
    sigma=0.22, density=0.12, pmax=0.70, R0_mean=10.0, R0_sd=0.25,
    niche_mode=:bimodal, mu_basal_centers=(0.5,0.75), mu_basal_sd=0.04,
    b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04
)
sum_clust = sweep_ensemble(
    pool_seeds, mask_seeds; 
    kind=:clustered, grid=grid, 
    τ=0.55, S=200, 
    basal_frac=0.35, nseeds_cluster=1, metric=metric,
    sigma=0.22, density=0.12, pmax=0.70, R0_mean=10.0, R0_sd=0.25,
    niche_mode=:bimodal, mu_basal_centers=(0.5,0.75), mu_basal_sd=0.04,
    b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04
)
# @show length(sum_rand.loss), length(sum_rand.ΔBSH_lo), length(sum_rand.ΔBSH_hi)
# @show length(sum_clust.loss), length(sum_clust.ΔBSH_lo), length(sum_clust.ΔBSH_hi)
# @show unique(sum_clust.loss)
# @show any(!isfinite.(sum_clust.ΔBSH_lo)), any(!isfinite.(sum_clust.ΔBSH_hi)), any(!isfinite.(sum_clust.ΔBSH_mean))
# @show length(sum_clust.loss), length(sum_clust.ΔBSH_lo), length(sum_clust.ΔBSH_hi), length(sum_clust.ΔBSH_mean)

# Plot with ribbons
begin
    fig = Figure(; size=(860,420))
    ax  = Axis(fig[1,1], xlabel="Area lost ($(string(metric)))", ylabel="ΔBSH (mean over consumers)")
    p = sortperm(sum_rand.loss)
    # band!(ax, sum_rand.loss[p], sum_rand.ΔBSH_lo[p], sum_rand.ΔBSH_hi[p]; transparency=true, label="Random (10–90%)")
    lines!(ax, sum_rand.loss[p], sum_rand.ΔBSH_mean[p])

    p = sortperm(sum_clust.loss)
    # band!(ax, sum_clust.loss[p], sum_clust.ΔBSH_lo[p], sum_clust.ΔBSH_hi[p]; transparency=true, label="Clustered (10–90%)")
    lines!(ax, sum_clust.loss[p], sum_clust.ΔBSH_mean[p])

    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    # axislegend(ax, position=:lb)
    display(fig)
end


