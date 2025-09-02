# --- keep-mask dispatcher (uses your existing functions + the new hotspot) ---
function keepmask_for(kind::Symbol, grid::Grid, pool::SpeciesPool;
                      τ::Float64, keep::Float64, nseeds_cluster::Int, seed::Int,
                      hotspot_power::Float64=2.5, hotspot_what::Symbol=:cons)
    if kind === :random
        return random_mask(grid.C, keep; seed=seed)
    elseif kind === :clustered
        return clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=seed)
    elseif kind === :hotspot
        return consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep,
                                     power=hotspot_power, what=hotspot_what, seed=seed)
    else
        error("unknown kind = $kind")
    end
end

# --- “blank out” columns instead of dropping (keeps the original denominator) ---
# returns a matrix the same width as Z, with removed cells set to false
function blank_mask(Z::BitMatrix, keep::BitVector)
    Zb = falses(size(Z))
    Zb[:, keep] .= Z[:, keep]
    return Zb
end

"""
sweep_ensemble(pool_seed_list, mask_seed_list;
               kind=:random, grid, τ=0.5, loss_fracs=0:0.05:0.9,
               S, basal_frac, nseeds_cluster=6,
               metric=:fraction,                      # :fraction or :area
               # metaweb / niche parameters (forwarded to build_pool)
               R0_mean=16.0, R0_sd=0.6, sigma=0.55, density=0.35, pmax=0.9,
               niche_mode=:uniform, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.06,
               b0_basal=0.10, bspread_basal=0.04, b0_cons=0.12, bspread_cons=0.04,
               # hotspot options
               hotspot_power=2.5, hotspot_what=:cons)
→ returns Summary(loss, ΔBSH_mean, ΔBSH_lo, ΔBSH_hi)
"""
function sweep_ensemble(pool_seed_list, mask_seed_list;
        kind::Symbol=:random,
        grid::Grid, τ::Float64=0.5,
        loss_fracs=0.0:0.05:0.9,
        S::Int, basal_frac::Float64, nseeds_cluster::Int=6,
        metric::Symbol=:fraction,
        R0_mean=16.0, R0_sd=0.6, sigma=0.55, density=0.35, pmax=0.9,
        niche_mode::Symbol=:uniform, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.06,
        b0_basal=0.10, bspread_basal=0.04, b0_cons=0.12, bspread_cons=0.04,
        hotspot_power=2.5, hotspot_what::Symbol=:cons)

    μ = Float64[]; lo = Float64[]; hi = Float64[]

    for f in loss_fracs
        keep = 1.0 - f
        vals = Float64[]

        for ps in pool_seed_list, ms in mask_seed_list
            # pool for this replicate
            pool = build_pool(S;
                basal_frac=basal_frac, seed=ps,
                R0_mean=R0_mean, R0_sd=R0_sd, sigma=sigma, density=density, pmax=pmax,
                niche_mode=niche_mode, mu_basal_centers=mu_basal_centers, mu_basal_sd=mu_basal_sd,
                b0_basal=b0_basal, bspread_basal=bspread_basal, b0_cons=b0_cons, bspread_cons=bspread_cons
            )

            # baseline on full grid
            Zfull = climate_pass(pool, grid; τ=τ)
            Pfull = assemble(Zfull, pool)
            cons  = .!pool.basal
            C0    = size(Zfull, 2)

            # choose keep mask for this replicate (may depend on pool)
            keepmask = keepmask_for(kind, grid, pool;
                                    τ=τ, keep=keep, nseeds_cluster=nseeds_cluster, seed=ms,
                                    hotspot_power=hotspot_power, hotspot_what=hotspot_what)

            if metric === :fraction
                # FRACTION basis: denominator = original C0
                # blank out removed cells so denominator stays C0
                Z1 = blank_mask(Zfull, keepmask)
                P1 = assemble(Z1, pool)
                base_BSH = mean(bsh1_per_species(Pfull, Zfull, pool)[cons])           # /C0
                BSH1     = mean(bsh1_per_species(P1,   Z1,   pool)[cons])             # /C0 (unchanged denom due to blanking)
                push!(vals, BSH1 - base_BSH)

            elseif metric === :area
                # AREA basis: denominator = kept cells; baseline evaluated on the kept set
                Z0k = Zfull[:, keepmask]; P0k = Pfull[:, keepmask]     # baseline on kept area
                Z1d = Z0k;                     # same set of kept columns…
                P1d = assemble(Z1d, pool)      # …but re-assembled (identical to P0k; kept explicit for clarity)
                base_BSH = mean(bsh1_per_species(P0k, Z0k, pool)[cons]) # /nkeep
                BSH1     = mean(bsh1_per_species(P1d, Z1d, pool)[cons]) # /nkeep
                push!(vals, BSH1 - base_BSH)

            else
                error("unknown metric = $metric (use :fraction or :area)")
            end
        end

        push!(μ,  mean(vals))
        push!(lo, quantile(vals, 0.10))
        push!(hi, quantile(vals, 0.90))
    end

    return Summary(collect(loss_fracs), μ, lo, hi)
end

# grid + parameter set
grid = make_grid(60, 60; seed=11)
pool_seeds = 1:5
mask_seeds = 1:30
S = 200; basal_frac = 0.35
τ = 0.55

# choose basis here
metric = :fraction          # or :area

sum_rand  = sweep_ensemble(pool_seeds, mask_seeds; kind=:random,   grid=grid, τ=τ,
                           S=S, basal_frac=basal_frac, metric=metric)
sum_clust = sweep_ensemble(pool_seeds, mask_seeds; kind=:clustered, grid=grid, τ=τ,
                           S=S, basal_frac=basal_frac, metric=metric, nseeds_cluster=6)
sum_hot   = sweep_ensemble(pool_seeds, mask_seeds; kind=:hotspot,  grid=grid, τ=τ,
                           S=S, basal_frac=basal_frac, metric=metric,
                           hotspot_power=2.5, hotspot_what=:cons)

# quick plot
begin
    fig = Figure(; size=(900,320))
    ax  = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="ΔBSH (mean over consumers)",
               title = "ΔBSH vs loss — $(uppercase(string(metric))) basis")

    for (s, col, lab) in ((sum_rand,  :steelblue, "Random"),
                          (sum_clust, :orange,    "Clustered"),
                          (sum_hot,   :seagreen,  "Hotspot"))
        p = sortperm(s.loss)
        band!(ax, s.loss[p], s.ΔBSH_lo[p], s.ΔBSH_hi[p]; color=(col, 0.15))
        lines!(ax, s.loss[p], s.ΔBSH_mean[p]; color=col, label=lab, linewidth=2)
    end
    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax, position=:lb)
    display(fig)
end

