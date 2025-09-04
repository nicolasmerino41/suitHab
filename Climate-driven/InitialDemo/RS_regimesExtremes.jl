"More extreme 2×2: low/high Redundancy × low/high Synchrony."
function RS_regimes_extreme()
    Dict(
        # Few prey per consumer (narrow & thin diets), basal niches diffuse
        :lowR_lowS  => (
            sigma=0.12, density=0.05, pmax=0.50,            # very low redundancy
            niche_mode=:uniform, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.15, # low synchrony (spread out)
            b0_basal=0.14, bspread_basal=0.05,              # broad basal breadths
            b0_cons=0.12, bspread_cons=0.04,
            R0_mean=10.0, R0_sd=0.15
        ),

        # Few prey per consumer, basal niches tightly clumped in two guilds
        :lowR_highS => (
            sigma=0.12, density=0.05, pmax=0.50,            # very low redundancy
            niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.01, # HIGH synchrony
            b0_basal=0.05, bspread_basal=0.01,              # narrow basal breadths
            b0_cons=0.10, bspread_cons=0.03,
            R0_mean=10.0, R0_sd=0.15
        ),

        # Many prey per consumer (wide & dense diets), basal niches diffuse
        :highR_lowS => (
            sigma=0.70, density=0.60, pmax=1.00,            # very high redundancy
            niche_mode=:uniform, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.15, # low synchrony
            b0_basal=0.16, bspread_basal=0.06,              # broad basal breadths
            b0_cons=0.18, bspread_cons=0.06,
            R0_mean=10.0, R0_sd=0.60
        ),

        # Many prey per consumer, basal niches tightly clumped
        :highR_highS=> (
            sigma=0.70, density=0.60, pmax=1.00,            # very high redundancy
            niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.01, # HIGH synchrony
            b0_basal=0.05, bspread_basal=0.01,              # narrow basal breadths
            b0_cons=0.18, bspread_cons=0.06,
            R0_mean=10.0, R0_sd=0.60
        )
    )
end

res = run_RS_once(; grid=grid, keep_frac=0.5, τ=0.65, S=500, basal_frac=0.35,
                  nseeds_cluster=1, regimes=RS_regimes_extreme())
plot_RS_bars(res)  # from earlier
