keep = 0.6
score = prey_vulnerability_score(pool, grid; τ=τ)
km_r = random_mask(grid.C, keep; seed=101)
km_h = hotspot_clustered_mask_bestfirst(grid, keep; score=score, nseeds=1, seed=303)

plot_TL_partition(pool, grid; τ=τ, keepmask=km_r,
    title_str="ΔBSH by TL — RANDOM (LowR–HighS)")
plot_TL_partition(pool, grid; τ=τ, keepmask=km_h,
    title_str="ΔBSH by TL — HOTSPOT (LowR–HighS)")
