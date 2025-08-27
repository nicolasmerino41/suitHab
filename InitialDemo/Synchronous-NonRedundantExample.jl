sum_rand = sweep_ensemble(1:5, 1:30; kind=:random, grid=grid, τ=0.5, S=200, basal_frac=0.45)
sum_clust = sweep_ensemble(
    1:5, 1:30; 
    kind=:clustered, grid=grid, τ=0.55,
    S=200, basal_frac=0.35, nseeds_cluster=1, metric=:area,
    sigma=0.22, density=0.12, pmax=0.70, R0_mean=10.0, R0_sd=0.25,
    niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.04,
    b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04
)

