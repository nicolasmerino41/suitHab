"Return a Dict of 4 regimes (low/high redundancy × low/high synchrony)."
function RS_regimes()
    Dict(
        :lowR_lowS  => (sigma=0.22, density=0.12, pmax=0.70, niche_mode=:uniform, mu_basal_sd=0.08),
        :lowR_highS => (sigma=0.22, density=0.12, pmax=0.70, niche_mode=:bimodal, mu_basal_sd=0.03),
        :highR_lowS => (sigma=0.45, density=0.35, pmax=0.90, niche_mode=:uniform, mu_basal_sd=0.08),
        :highR_highS=> (sigma=0.45, density=0.35, pmax=0.90, niche_mode=:bimodal, mu_basal_sd=0.03)
    )
end

"""
run_RS_once(grid; keep_frac, τ, S, basal_frac, nseeds_cluster=1, regimes=RS_regimes())

For each regime, builds one pool and computes area-basis decomposition
means over consumers for RANDOM and CLUSTERED masks at the same loss.
Returns a dict of NamedTuples.
"""
function run_RS_once(; grid::Grid, keep_frac::Float64=0.6, τ::Float64=0.62,
                      S::Int=500, basal_frac::Float64=0.35, nseeds_cluster::Int=1,
                      seed_pool::Int=1, regimes=RS_regimes())
    km_r = random_mask(grid.C, keep_frac; seed=101)
    km_c = clustered_mask(grid, keep_frac; nseeds=nseeds_cluster, seed=202)

    results = Dict{Symbol,Any}()
    for (name, kw) in regimes
        pool = build_pool(S; basal_frac=basal_frac, seed=seed_pool,
                          sigma=kw[:sigma], density=kw[:density], pmax=kw[:pmax],
                          niche_mode=kw[:niche_mode], mu_basal_sd=kw[:mu_basal_sd])

        dr = decomp_at_mask(pool, grid; τ=τ, keepmask=km_r)
        dc = decomp_at_mask(pool, grid; τ=τ, keepmask=km_c)

        results[name] = (
            random   = dr.agg,    # (dB, dAcl, dInt, dSyn)
            clustered= dc.agg,
            pool     = pool
        )
    end
    return results
end

"Plot climate / interaction / synergy bars for the 4 regimes (random vs clustered)."
function plot_RS_bars(results::Dict{Symbol,Any})
    order = [:lowR_lowS, :lowR_highS, :highR_lowS, :highR_highS]
    titles = Dict(
        :lowR_lowS=>"Low R, Low S", :lowR_highS=>"Low R, High S",
        :highR_lowS=>"High R, Low S", :highR_highS=>"High R, High S"
    )

    begin
        fig = Figure(resolution=(980,620))
        i = 1
        for key in order
            r = results[key]
            ax = Axis(fig[div(i-1,2)+1, mod(i-1,2)+1],
                      title=titles[key], xticks=(1:3, ["Climate","Interact","Synergy"]),
                      ylabel = (i in (1,3)) ? "Contribution" : "")
            x = 1:3; off=0.18
            barplot!(ax, x .- off, [r.random.dAcl, r.random.dInt, r.random.dSyn];
                     width=0.35, color=:dodgerblue, label="Random")
            barplot!(ax, x .+ off, [r.clustered.dAcl, r.clustered.dInt, r.clustered.dSyn];
                     width=0.35, color=:orange,     label="Clustered")
            hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
            if i==2; axislegend(ax, position=:rt); end
            i += 1
        end
        display(fig)
    end
end

# choose loss & τ once (same for all comparisons)
keep = 0.6
τ    = 0.45

# 2×2 redundancy×synchrony bar comparison (consumers only)
res = run_RS_once(; grid=grid, keep_frac=keep, τ=τ, S=200, basal_frac=0.35, nseeds_cluster=1)
plot_RS_bars(res)

# Pick one regime (e.g., low redundancy + high synchrony) and partition by TL
pool_LRHS = res[:lowR_highS].pool
km_c = clustered_mask(grid, keep; nseeds=1, seed=202)
plot_TL_partition(
    pool, grid; τ=τ, keepmask=km_c,
    title_str="ΔBSH by TL (clustered, LowR–HighS)"
    )
pool = build_pool(300;
    basal_frac=0.35, seed=1,
    # deeper hierarchy
    R0_mean=3.0, R0_sd=0.10,
    sigma=0.14, density=0.08, pmax=0.65,
    # keep synchrony as you like (e.g., bimodal & tight)
    niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.02,
    # ensure consumers have some climate area
    b0_cons=0.14, bspread_cons=0.05
)
τ = 0.52
