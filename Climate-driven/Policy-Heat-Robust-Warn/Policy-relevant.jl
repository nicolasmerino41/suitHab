include("Helpers.jl")
# --- RUN & PLOT ---
begin
    grid   = make_grid(60,60; seed=11)
    pool   = build_pool(220; basal_frac=0.35, seed=1,
                        sigma=0.22, density=0.12, pmax=0.70,
                        niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.03,
                        b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04)
    τ      = 0.58
    losses = collect(0.15:0.05:0.60)
    Z0     = climate_pass(pool, grid; τ=τ)
    P0     = assemble(Z0, pool)
    A_cell, prey_cell = cell_scores(Z0, P0, pool)

    function keep_top(score; keep_frac)
        C=length(score); k=round(Int, keep_frac*C); idx=sortperm(score; rev=true)[1:k]
        keep=falses(C); keep[idx].=true; keep
    end

    D_rand = Float64[]; D_cli = Float64[]; D_int = Float64[]
    for loss in losses
        keep = 1 - loss
        kmR  = random_mask(grid.C, keep; seed=101)
        kmC  = keep_top(A_cell; keep_frac=keep)               # climate-only
        kmI  = keep_top(A_cell .* prey_cell; keep_frac=keep)  # interaction-aware

        dR = collect_damage_series([decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmR)])
        dC = collect_damage_series([decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmC)])
        dI = collect_damage_series([decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmI)])

        push!(D_rand, dR.D_tot[1]); push!(D_cli, dC.D_tot[1]); push!(D_int, dI.D_tot[1])
    end

    AD_cli = @. D_rand - D_cli    # avoided damage vs random (↑ better)
    AD_int = @. D_rand - D_int

    fig = Figure(; size=(1080,420))

    ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Damage D = −ΔBSH",
               title="Prioritization performance")
    lines!(ax1, losses, D_rand, label="Random")
    lines!(ax1, losses, D_cli,  label="Climate-only score")
    lines!(ax1, losses, D_int,  label="Interaction-aware score")
    hlines!(ax1, [0.0], color=(:gray,0.5), linestyle=:dash); axislegend(ax1, position=:lt)

    ax2 = Axis(fig[1,2], xlabel="Area lost (fraction)", ylabel="Avoided damage vs random (↑ better)")
    lines!(ax2, losses, AD_cli, label="Climate-only")
    lines!(ax2, losses, AD_int, label="Interaction-aware")
    hlines!(ax2, [0.0], color=(:gray,0.5), linestyle=:dash); axislegend(ax2, position=:lt)

    display(fig)
end

