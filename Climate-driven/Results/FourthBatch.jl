# ---- NaN-safe fraction-basis decomposition of mean ΔBSH over consumers ----
struct Decomp{T}; dB::T; dA_only::T; dInt_only::T; synergy::T; end

decompose_delta(B0,A0,p0,B1,A1,p1) = begin
    dB = B1 - B0
    dA = A1 - A0
    dp = p1 - p0
    Abar = 0.5*(A0 + A1)
    pbar = 0.5*(p0 + p1)
    dA_only = pbar * dA
    dInt_only = Abar * dp
    synergy = dB - dA_only - dInt_only
    Decomp(dB, dA_only, dInt_only, synergy)
end

function decomp_at_mask_fraction_safe(pool::SpeciesPool, grid::Grid; τ::Float64, keepmask::BitVector)
    # build before/after
    Z0 = climate_pass(pool, grid; τ=τ);  P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask);       P1 = assemble(Z1, pool)

    cons = .!pool.basal

    A0 = climate_area_per_species(Z0)[cons]
    A1 = climate_area_per_species(Z1)[cons]

    p0 = bsh1_cond_per_species(P0, Z0, pool)[cons]
    p1 = bsh1_cond_per_species(P1, Z1, pool)[cons]

    # Treat undefined conditional probs as 0 (no suitable cells ⇒ no effect)
    replaceNaN!(v) = (for i in eachindex(v); if isnan(v[i]); v[i]=0.0; end; end; v)
    p0 = replaceNaN!(p0); p1 = replaceNaN!(p1)

    B0 = A0 .* p0
    B1 = A1 .* p1

    d = map(decompose_delta, B0, A0, p0, B1, A1, p1)

    # mean over consumers; also guard if any component still NaN (shouldn’t now)
    mean0(x) = mean(skipmissing((isnan(y) ? 0.0 : y) for y in x))
    (; clim = mean0(getfield.(d,:dA_only)),
       inter= mean0(getfield.(d,:dInt_only)),
       syn  = mean0(getfield.(d,:synergy)))
end

begin
    # --- regimes (edit as you like) ---
    regimes = [
        (; name="τ0.50–σ0.45–bB0.12", τ=0.50, sigma=0.45, density=0.35, pmax=0.90, niche_mode=:uniform,  b0_basal=0.12, bspread_basal=0.03, b0_cons=0.12, bspread_cons=0.04, mu_basal_sd=0.08),
        (; name="τ0.50–σ0.22–bB0.12", τ=0.50, sigma=0.22, density=0.12, pmax=0.70, niche_mode=:uniform,  b0_basal=0.12, bspread_basal=0.03, b0_cons=0.12, bspread_cons=0.04, mu_basal_sd=0.08),
        (; name="τ0.55–σ0.22–bimodal", τ=0.55, sigma=0.22, density=0.12, pmax=0.70, niche_mode=:bimodal,  b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04, mu_basal_sd=0.04),
        (; name="τ0.60–σ0.22–bimodal", τ=0.60, sigma=0.22, density=0.12, pmax=0.70, niche_mode=:bimodal,  b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04, mu_basal_sd=0.04),
        (; name="τ0.55–σ0.35–bB0.10", τ=0.55, sigma=0.35, density=0.25, pmax=0.80, niche_mode=:uniform,  b0_basal=0.10, bspread_basal=0.03, b0_cons=0.12, bspread_cons=0.04, mu_basal_sd=0.06),
    ]

    grid = make_grid(60,60; seed=11)
    S = 200; basal_frac=0.35
    loss = 0.50; keep = 1.0 - loss
    hotspot_power = 2.5
    nseeds_cluster = 1

    function scenario_decomp(pool, grid; τ, keep)
        kmR = random_mask(grid.C, keep; seed=101)
        kmC = clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=101)
        kmH = consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep, power=hotspot_power, what=:cons, seed=101)
        dR = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmR)
        dC = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmC)
        dH = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmH)
        return Dict(:Random=>dR, :Clustered=>dC, :Hotspot=>dH)
    end

    results = Dict{String,Dict{Symbol,NamedTuple}}()
    for reg in regimes
        pool = build_pool(S;
            basal_frac=basal_frac, seed=1,
            sigma=reg.sigma, density=reg.density, pmax=reg.pmax,
            niche_mode=reg.niche_mode, mu_basal_centers=(0.25,0.75), mu_basal_sd=reg.mu_basal_sd,
            b0_basal=reg.b0_basal, bspread_basal=reg.bspread_basal,
            b0_cons=reg.b0_cons,   bspread_cons=reg.bspread_cons
        )
        results[reg.name] = scenario_decomp(pool, grid; τ=reg.τ, keep=keep)
    end

    # ---- plot: stacked bars (polygon-based, version-agnostic) ----
    fig = Figure(; size=(1200, 420))
    scen_order = [:Random, :Clustered, :Hotspot]
    cols = [:dodgerblue, :orange, :forestgreen]

    for (j, scen) in enumerate(scen_order)
        ax = Axis(fig[1, j], title = String(scen),
                xlabel = "Regime",
                ylabel = (j == 1 ? "ΔBSH (50% loss)" : ""))

        x = collect(1:length(regimes))

        clim  = Float64.([results[r.name][scen].clim  for r in regimes])
        inter = Float64.([results[r.name][scen].inter for r in regimes])
        syn   = Float64.([results[r.name][scen].syn   for r in regimes])

        Y = hcat(clim, inter, syn) |> Array  # N×3

        stackedbars_poly!(ax, x, Y; width=0.8, colors=cols)

        ax.xticks = (x, [r.name for r in regimes])
        ax.xticklabelrotation[] = π/6
        hlines!(ax, [0.0], color = (:gray, 0.5), linestyle = :dash)

        if j == 3
            axislegend(ax,
                [PolyElement(color=cols[1]),
                PolyElement(color=cols[2]),
                PolyElement(color=cols[3])],
                ["Climate", "Interaction", "Synergy"]; position = :lt)
        end
    end

    display(fig)
end

begin
    # knobs (match your earlier ones)
    grid = make_grid(60,60; seed=11)
    S = 200; basal_frac = 0.35
    keep = 0.5                              # 50% kept => 50% loss
    τ_vals = 0.45:0.03:0.65
    σ_vals = 0.20:0.05:0.50
    hotspot_power = 2.5
    pool_seed = 1
    mask_seed = 101
    nseeds_cluster = 1

    # allocate: rows = σ index, cols = τ index
    inter_share_R = fill(NaN, length(σ_vals), length(τ_vals))
    inter_share_C = fill(NaN, length(σ_vals), length(τ_vals))
    inter_share_H = fill(NaN, length(σ_vals), length(τ_vals))

    for (i, σ) in enumerate(σ_vals), (j, τv) in enumerate(τ_vals)
        pool = build_pool(S;
            basal_frac=basal_frac, seed=pool_seed,
            sigma=σ, density=0.15, pmax=0.75,
            niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.04,
            b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04
        )

        # masks
        kmR = random_mask(grid.C, keep; seed=mask_seed)
        kmC = clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=mask_seed)
        kmH = consumer_hotspot_mask(grid, pool; τ=τv, keep_frac=keep,
                                    power=hotspot_power, what=:cons, seed=mask_seed)

        # fraction-basis decomposition (use your NaN-safe function)
        dR = decomp_at_mask_fraction_safe(pool, grid; τ=τv, keepmask=kmR)
        dC = decomp_at_mask_fraction_safe(pool, grid; τ=τv, keepmask=kmC)
        dH = decomp_at_mask_fraction_safe(pool, grid; τ=τv, keepmask=kmH)

        # interaction share ∈ [-1,1]
        denom(x) = abs(x.clim) + abs(x.inter) + abs(x.syn) + 1e-12
        inter_share_R[i,j] = dR.inter / denom(dR)
        inter_share_C[i,j] = dC.inter / denom(dC)
        inter_share_H[i,j] = dH.inter / denom(dH)
    end

    # plot the 3 heatmaps with same colorrange
    fig = Figure(; size=(1080, 380))
    τvec, σvec = collect(τ_vals), collect(σ_vals)
    rng = (-1.0, 1.0)

    axR = Axis(fig[1,1], title="Random",   xlabel="τ (climate threshold)", ylabel="σ (diet width)")
    hmR = heatmap!(axR, τvec, σvec, inter_share_R)#; colorrange=rng)

    axC = Axis(fig[1,2], title="Clustered", xlabel="τ (climate threshold)")
    hmC = heatmap!(axC, τvec, σvec, inter_share_C)#; colorrange=rng)

    axH = Axis(fig[1,3], title="Hotspot",   xlabel="τ (climate threshold)")
    hmH = heatmap!(axH, τvec, σvec, inter_share_H)#; colorrange=rng)

    Colorbar(fig[:, end+1], hmH, label="Interaction / Total")
    display(fig)
end
