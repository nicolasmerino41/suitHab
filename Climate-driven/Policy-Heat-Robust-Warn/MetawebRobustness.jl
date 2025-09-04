begin
    grid   = make_grid(60,60; seed=11)
    τ      = 0.58
    keep   = 0.50
    pool0  = build_pool(220; basal_frac=0.35, seed=1,
                        sigma=0.22, density=0.12, pmax=0.70,
                        niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.03,
                        b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04)
    Z0   = climate_pass(pool0, grid; τ=τ); P0 = assemble(Z0, pool0)
    Acel, Pcel = cell_scores(Z0, P0, pool0)

    fracs = 0.0:0.1:1.0
    shares = Float64[]; absEDi = Float64[]

    for f in fracs
        pool = rewire_metaweb_fraction(pool0; frac=f, seed=202)
        # rebuild scores fairly for each pool
        Z  = climate_pass(pool, grid; τ=τ); P = assemble(Z, pool)
        A, prey = cell_scores(Z, P, pool)

        kmR = random_mask(grid.C, keep; seed=101)
        kmH = a_matched_hotspot_mask(A, prey; keep_frac=keep, nbins=10, seed=101)

        dR = collect_damage_series([decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmR)])
        dH = collect_damage_series([decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmH)])
        ex = excess_from_damage(dR, dH)
        push!(shares, ex.share[1])
        push!(absEDi, ex.ED_int_abs[1])
    end

    fig = Figure(; size=(800,430))
    ax1 = Axis(fig[1,1], xlabel="Rewiring fraction (0 → 1)",
               ylabel="Interaction share of excess", title="Metaweb robustness")
    lines!(ax1, collect(fracs), shares)
    hlines!(ax1, [0.0, 1.0], color=(:gray,0.5), linestyle=:dash)

    ax2 = Axis(fig[2,1], xlabel="Rewiring fraction (0 → 1)",
               ylabel="|ED_int|", title="Magnitude of interaction excess")
    lines!(ax2, collect(fracs), absEDi)
    hlines!(ax2, [0.0], color=(:gray,0.5), linestyle=:dash)

    display(fig)
end
