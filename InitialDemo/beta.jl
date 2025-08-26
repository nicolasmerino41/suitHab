include("Functions.jl")
include("Knee.jl")
## MAKE THE GRID AND SPECIES POOL ##
grid = make_grid(60, 60; seed=11)

pool_seeds = 1:3
mask_seeds = 1:30
sum_rand  = sweep_ensemble(pool_seeds, mask_seeds; kind=:random, grid=grid, τ=0.7, S=300, basal_frac=0.2, metric=:fraction)
sum_clust = sweep_ensemble(pool_seeds, mask_seeds; kind=:clustered, grid=grid, τ=0.7, S=300, basal_frac=0.2, metric=:fraction)

@show length(sum_rand.loss), length(sum_rand.ΔBSH_lo), length(sum_rand.ΔBSH_hi)
@show length(sum_clust.loss), length(sum_clust.ΔBSH_lo), length(sum_clust.ΔBSH_hi)
@show unique(sum_clust.loss)
@show any(!isfinite.(sum_clust.ΔBSH_lo)), any(!isfinite.(sum_clust.ΔBSH_hi)), any(!isfinite.(sum_clust.ΔBSH_mean))
@show length(sum_clust.loss), length(sum_clust.ΔBSH_lo), length(sum_clust.ΔBSH_hi), length(sum_clust.ΔBSH_mean)

# Plot with ribbons
begin
    fig = Figure(; size=(860,420))
    ax  = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="ΔBSH (mean over consumers)")
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

# One loss level, decomposed and averaged across consumers
function decomp_at(pool, grid; τ=0.5, keepmask)
    Z0 = climate_pass(pool, grid; τ=τ)
    P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask)
    P1 = assemble(Z1, pool)

    cons = .!pool.basal
    A0 = climate_area_per_species(Z0)[cons]
    p0 = bsh1_cond_per_species(P0, Z0, pool)[cons]
    B0 = A0 .* p0

    A1 = climate_area_per_species(Z1)[cons]
    p1 = bsh1_cond_per_species(P1, Z1, pool)[cons]
    B1 = A1 .* p1

    d = map(decompose_delta, B0, A0, p0, B1, A1, p1)
    dB   = mean(getfield.(d, :dB))
    dAcl = mean(getfield.(d, :dA_part))
    dInt = mean(getfield.(d, :dp_part))
    dSyn = mean(getfield.(d, :synergy))
    (; dB, dAcl, dInt, dSyn)
end

# Example at 40% loss
keep = 0.6   # e.g., 40% loss
km_r = random_mask(grid.C, keep; seed=1)
km_c = clustered_mask(grid, keep; nseeds=1, seed=2)   # one big block
Z0 = climate_pass(pool, grid; τ=0.55)
basal_suit = vec(sum(Z0[pool.basal, :]; dims=1))  # basal richness by climate
top = sortperm(basal_suit, rev=true)[1:round(Int, 0.2*grid.C)]  # top 20%
seed_hot = rand(top)
km_hot = clustered_mask(grid, keep; nseeds=1, seed=seed_hot)     # use `seed` as starting index

pool_demo = build_pool(140; basal_frac=0.25, seed=3)
dr = decomp_at(pool_demo, grid; τ=0.6, keepmask=km_r)
dc = decomp_at(pool_demo, grid; τ=0.6, keepmask=km_c)

# Using your decomp_at() outputs dr (random) and dc (clustered)
begin
    using CairoMakie
    fig = Figure(size=(800,320))
    ax  = Axis(fig[1,1],
        title = "ΔBSH decomposition at 40% loss",
        ylabel = "Contribution",
        xticks = (1:3, ["Climate-only","Interaction-only","Synergy"])
    )

    x = 1:3
    off = 0.18
    bars_r = barplot!(ax, x .- off, [dr.dAcl, dr.dInt, dr.dSyn];
                      width=0.35, color=:dodgerblue, label="Random")
    bars_c = barplot!(ax, x .+ off, [dc.dAcl, dc.dInt, dc.dSyn];
                      width=0.35, color=:orange,    label="Clustered")

    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax, position=:rt)
    display(fig)
end
