"Make a copy of pool with interactions disabled for consumers."
function kill_interactions(pool::SpeciesPool)
    P2 = deepcopy(pool)
    for s in eachindex(P2.basal)
        if !P2.basal[s]
            P2.E[s] = Int[]       # no prey requirement
        end
    end
    return P2
end

"Decompose ΔBSH on a fraction basis for any keepmask (your NaN-safe version)."
decomp_at_mask_fraction_safe(pool, grid; τ, keepmask) = begin
    Z0 = climate_pass(pool, grid; τ=τ);  P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask);       P1 = assemble(Z1, pool)
    cons = .!pool.basal
    A0 = climate_area_per_species(Z0)[cons]
    A1 = climate_area_per_species(Z1)[cons]
    p0 = bsh1_cond_per_species(P0, Z0, pool)[cons]
    p1 = bsh1_cond_per_species(P1, Z1, pool)[cons]
    for v in (p0,p1); @inbounds for i in eachindex(v); if isnan(v[i]); v[i]=0.0; end; end; end
    B0, B1 = A0 .* p0, A1 .* p1
    dB = B1 .- B0; dA = A1 .- A0; dp = p1 .- p0
    Abar = 0.5 .* (A0 .+ A1); pbar = 0.5 .* (p0 .+ p1)
    dA_only  = pbar .* dA
    dInt_only= Abar .* dp
    syn      = dB .- dA_only .- dInt_only
    return (; clim=mean(dA_only), inter=mean(dInt_only), syn=mean(syn)) # (climate, interaction, synergy)
end

"Excess damage decomposition for a scenario vs Random at one loss fraction."
function excess_vs_random(pool, grid; τ, keep_frac, scenario::Symbol)
    kmR = random_mask(grid.C, keep_frac; seed=11)
    km  = scenario==:clustered ? clustered_mask(grid, keep_frac; nseeds=1, seed=11) :
          scenario==:hotspot   ? consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep_frac, power=2.5, what=:cons, seed=11) :
          error("scenario must be :clustered or :hotspot")
    climR, interR, synR = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmR)
    climS, interS, synS = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=km)
    return (; d_clim=climS-climR, d_inter=interS-interR, d_syn=synS-synR)
end

# >>> RUN
grid = make_grid(60,60; seed=11)
pool = build_pool(200; basal_frac=0.35, seed=1,
                  sigma=0.22, density=0.12, pmax=0.70,
                  niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.04,
                  b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04)
τ, keep = 0.55, 0.50

# Baseline world
exC = excess_vs_random(pool, grid; τ=τ, keep_frac=keep, scenario=:clustered)
exH = excess_vs_random(pool, grid; τ=τ, keep_frac=keep, scenario=:hotspot)

# No-interaction world
pool_C = kill_interactions(pool)
exC_noInt = excess_vs_random(pool_C, grid; τ=τ, keep_frac=keep, scenario=:clustered)
exH_noInt = excess_vs_random(pool_C, grid; τ=τ, keep_frac=keep, scenario=:hotspot)

@show exC exH exC_noInt exH_noInt

# --- helper: total ΔBSH from decomposition NamedTuple ---
totalΔ(d) = d.clim + d.inter + d.syn

function excess_curves_control1(pool, grid; τ, losses=0.05:0.05:0.8, nseeds_cluster=1, hotspot_power=2.5)
    exC = Float64[]; exC_no = Float64[]; exH = Float64[]; exH_no = Float64[]
    pool_no = kill_interactions(pool)
    for f in losses
        keep = 1.0 - f
        # masks
        kmR = random_mask(grid.C, keep; seed=11)
        kmC = clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=11)
        kmH = consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep, power=hotspot_power, what=:cons, seed=11)
        # with interactions
        dR  = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmR)
        dC  = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmC)
        dH  = decomp_at_mask_fraction_safe(pool, grid; τ=τ, keepmask=kmH)
        push!(exC, totalΔ(dC) - totalΔ(dR))
        push!(exH, totalΔ(dH) - totalΔ(dR))
        # no-interaction world
        kmH_no = consumer_hotspot_mask(grid, pool_no; τ=τ, keep_frac=keep, power=hotspot_power, what=:cons, seed=11)
        dRno = decomp_at_mask_fraction_safe(pool_no, grid; τ=τ, keepmask=kmR)
        dCno = decomp_at_mask_fraction_safe(pool_no, grid; τ=τ, keepmask=kmC)
        dHno = decomp_at_mask_fraction_safe(pool_no, grid; τ=τ, keepmask=kmH_no)
        push!(exC_no, totalΔ(dCno) - totalΔ(dRno))
        push!(exH_no, totalΔ(dHno) - totalΔ(dRno))
    end
    return (; losses=collect(losses), exC, exC_no, exH, exH_no)
end

# --- plot ---
begin
    grid = make_grid(60,60; seed=11)
    pool = build_pool(200; basal_frac=0.35, seed=1,
        sigma=0.2, density=0.1, pmax=0.6,
        niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.04,
        b0_basal=0.06, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04)

    τ = 0.65
    curves = excess_curves_control1(pool, grid; τ=τ, losses=0.05:0.05:0.8)

    fig = Figure(resolution=(900,330))
    for (j, lab, y1, y2) in (
            (1,"Clustered", curves.exC, curves.exC_no),
            (2,"Hotspot",   curves.exH, curves.exH_no)
        )
        ax = Axis(fig[1,j],
                title=lab,
                xlabel="Area lost (fraction)",
                ylabel=(j==1 ? "Excess ΔBSH (scenario − random)" : ""))

        l1 = lines!(ax, curves.losses, y1; linewidth=3, color=:dodgerblue)
        l2 = lines!(ax, curves.losses, y2; linestyle=:dash, linewidth=2, color=:gray)
        hlines!(ax, [0.0], color=(:black,0.3), linestyle=:dot)

        if j==2
            Legend(fig[1,3], [l1,l2], ["with interactions","no-interaction control"])
        end
    end

    display(fig)
end
