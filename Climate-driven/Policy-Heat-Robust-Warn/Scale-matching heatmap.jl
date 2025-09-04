include("Helpers.jl")
begin
    grid   = make_grid(60,60; seed=11)
    τ      = 0.58
    keep   = 0.50
    σ_vals = [0.20, 0.25, 0.30, 0.35, 0.40]
    seeds  = [2, 4, 6, 9, 12, 15]   # fewer seeds → larger clusters
    inter_share = fill(NaN, length(σ_vals), length(seeds))

    for (i, σ) in enumerate(σ_vals)
        pool = build_pool(220; basal_frac=0.35, seed=1,
                          sigma=σ, density=0.12, pmax=0.70,
                          niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.03,
                          b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04)
        kmR = random_mask(grid.C, keep; seed=101)
        dR  = collect_damage_series([decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmR)])[1]

        for (j, nseeds) in enumerate(seeds)
            kmC = clustered_mask(grid, keep; nseeds=nseeds, seed=101)
            dC  = collect_damage_series([decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmC)])[1]

            R = collect_damage_series([decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmR)])
            C = collect_damage_series([decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmC)])
            ex = excess_from_damage(R, C)

            inter_share[i,j] = isnan(ex.share[1]) ? NaN : ex.share[1]
        end
    end

    fig = Figure(; size=(820,480))
    ax  = Axis(fig[1,1], xlabel="Cluster seeds (↓ ⇒ larger clusters)", ylabel="Diet width σ",
               title="Interaction share of excess (clustered @ 50% loss)")
    hm  = heatmap!(ax, 1:length(seeds), σ_vals, inter_share; nan_color=:gray90)
    Colorbar(fig[1,2], hm, label="interaction share (0–1)")
    ax.xticks = (1:length(seeds), string.(seeds))
    display(fig)
end

