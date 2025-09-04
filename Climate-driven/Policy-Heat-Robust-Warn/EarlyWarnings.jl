# Compute ΔA & Δp across consumers for a scenario at a given loss
function deltas_per_consumer(pool, grid; τ, keepmask)
    Z0 = climate_pass(pool, grid; τ=τ); P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask);       P1 = assemble(Z1, pool)

    cons = .!pool.basal
    A0 = climate_area_per_species(Z0)[cons]
    A1 = climate_area_per_species(Z1)[cons]
    p0 = bsh1_cond_per_species(P0, Z0, pool)[cons]
    p1 = bsh1_cond_per_species(P1, Z1, pool)[cons]
    nan_to_zero!(p0); nan_to_zero!(p1)
    return (dA = A1 .- A0, dp = p1 .- p0)
end

begin
    grid   = make_grid(60,60; seed=11)
    pool   = build_pool(220; basal_frac=0.35, seed=1,
                        sigma=0.22, density=0.12, pmax=0.70,
                        niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.03,
                        b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04)
    τ      = 0.58
    keep   = 0.50

    # A-matched hotspot
    Z0 = climate_pass(pool, grid; τ=τ); P0 = assemble(Z0, pool)
    A, prey = cell_scores(Z0, P0, pool)
    kmH = a_matched_hotspot_mask(A, prey; keep_frac=keep, nbins=10, seed=101)

    Δ = deltas_per_consumer(pool, grid; τ=τ, keepmask=kmH)

    # OLS slope & Spearman rho
    slope = cov(Δ.dA, Δ.dp) / (var(Δ.dA) + 1e-12)
    function spearman_rho(x,y)
        rx = sortperm(sortperm(x)) .+ 1
        ry = sortperm(sortperm(y)) .+ 1
        cov(rx,ry)/(std(rx)*std(ry)+1e-12)
    end
    rho = spearman_rho(Δ.dA, Δ.dp)

    fig = Figure(; size=(840,380))
    ax  = Axis(fig[1,1], xlabel="ΔA (change in climate-suitable fraction)",
               ylabel="Δp (change in prey support)",
               title="Δp vs ΔA @ 50% loss — slope=$(round(slope,digits=3)), ρ=$(round(rho,digits=3))")
    scatter!(ax, Δ.dA, Δ.dp)
    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    vlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    lines!(ax, sort(Δ.dA), slope .* sort(Δ.dA))
    display(fig)
end
