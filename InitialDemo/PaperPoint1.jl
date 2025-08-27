# Hotspot clustering + decomposition vs loss
# ---------- helpers you can paste as-is ----------

# 4-neighborhood on the grid (re-include if you don't already have it)
# --- Frontier helpers (no deps) ---
neighbors4(i, nx, ny) = begin
    x = (i - 1) % nx + 1
    y = (i - 1) ÷ nx + 1
    nb = Int[]
    if x > 1;   push!(nb, i - 1)    end
    if x < nx;  push!(nb, i + 1)    end
    if y > 1;   push!(nb, i - nx)   end
    if y < ny;  push!(nb, i + nx)   end
    nb
end

"""
prey_support_score(pool, grid; τ)
Return a vector S(c) counting, for each cell c, how many consumers have
(i) climate-suitable and (ii) ≥1 prey present in PRE-LOSS assembly.
"""
# Pre-loss prey-support score (count of consumers with ≥1 prey in c)
function prey_support_score(pool::SpeciesPool, grid::Grid; τ::Float64=0.64)
    Z0 = climate_pass(pool, grid; τ=τ)
    P0 = assemble(Z0, pool)
    cons = .!pool.basal
    S = zeros(Int, grid.C)
    @inbounds for c in 1:grid.C
        cnt = 0
        for s in findall(cons)
            if Z0[s,c]
                has_pre = false
                for q in pool.E[s]; if P0[q,c]; has_pre = true; break; end; end
                cnt += has_pre ? 1 : 0
            end
        end
        S[c] = cnt
    end
    S
end

# Optional: "vulnerability" score (emphasize low redundancy)
function prey_vulnerability_score(pool::SpeciesPool, grid::Grid; τ::Float64=0.64)
    Z0 = climate_pass(pool, grid; τ=τ)
    P0 = assemble(Z0, pool)
    cons = .!pool.basal
    S = zeros(Float64, grid.C)
    @inbounds for c in 1:grid.C
        tot = 0.0
        for s in findall(cons)
            if Z0[s,c]
                k = 0
                for q in pool.E[s]; k += (P0[q,c] ? 1 : 0); end
                if k > 0
                    tot += 1.0 / k  # heavier weight when fewer prey present
                end
            end
        end
        S[c] = tot
    end
    S
end

"""
hotspot_clustered_mask(grid, keep_frac; score, nseeds=1, seed=0)

Grow a single (or few) BFS cluster starting at the top 'score' cells
until the target number of cells is removed.
Returns a Bool keep-mask (true = keep).
"""
function hotspot_clustered_mask(grid::Grid, keep_frac::Float64;
                                score::Vector{<:Real},
                                nseeds::Int=1, seed::Int=0)
    Random.seed!(seed)
    C  = grid.C
    tgt_remove = C - round(Int, keep_frac*C)
    removed = falses(C)

    # choose seed(s) as the highest-scoring cells (break ties randomly)
    order = sortperm(score; rev=true, by=identity)
    seeds = order[1:nseeds]

    Q = collect(seeds)  # BFS queue
    ptr = 1; removed_ct = 0

    while removed_ct < tgt_remove
        if ptr > length(Q)
            # fall back: push another high-score not yet removed or queued
            for i in order
                if !removed[i] && !(i in Q)  # cheap check is fine
                    push!(Q, i); break
                end
            end
        end
        v = Q[ptr]; ptr += 1
        removed[v] && continue
        removed[v] = true
        removed_ct += 1
        for nb in neighbors4(v, grid.nx, grid.ny)
            if !removed[nb]
                push!(Q, nb)
            end
        end
    end
    return .!removed
end

"""
Hotspot cluster with best-first growth (keeps hugging high-score cells).
- score: Vector{<:Real} (e.g., prey_support_score or prey_vulnerability_score)
- nseeds: number of initial hotspots (use 1 for maximal coherence)
"""
function hotspot_clustered_mask_bestfirst(grid::Grid, keep_frac::Float64;
                                          score::AbstractVector,
                                          nseeds::Int=1, seed::Int=0)
    Random.seed!(seed)
    C  = grid.C
    target_remove = C - round(Int, keep_frac*C)
    removed = falses(C)
    queued  = falses(C)

    order = sortperm(score; rev=true)       # highest score first
    seeds  = order[1:nseeds]

    frontier = Int[]
    for s in seeds
        push!(frontier, s); queued[s] = true
    end

    removed_ct = 0
    while removed_ct < target_remove
        if isempty(frontier)
            # if frontier dries up, re-seed with next best cell not used yet
            for i in order
                if !removed[i] && !queued[i]
                    push!(frontier, i); queued[i] = true; break
                end
            end
        end
        # pick the frontier cell with max score
        imax = findmax(score[frontier])[2]
        v = frontier[imax]
        frontier[imax] = frontier[end]; pop!(frontier)

        if removed[v]; continue; end
        removed[v] = true
        removed_ct += 1

        # add neighbors into frontier
        for nb in neighbors4(v, grid.nx, grid.ny)
            if !removed[nb] && !queued[nb]
                push!(frontier, nb); queued[nb] = true
            end
        end
    end
    return .!removed   # keep mask
end

# ---------- decomposition vs loss for 3 strategies ----------

"""
decomp_series_3strat(pool, grid; τ, losses, nseeds_cluster=1)

For each loss fraction, compute area-basis mean over consumers of:
  climate-only (dAcl), interaction-only (dInt), synergy (dSyn)
under three strategies: random, clustered, hotspot-clustered.
Returns a NamedTuple with fields .losses and .rand/.clust/.hotspot,
each containing (dAcl, dInt, dSyn) vectors.
"""
function decomp_series_3strat(pool::SpeciesPool, grid::Grid;
                              τ::Float64=0.64,
                              losses=collect(0.0:0.05:0.8),
                              nseeds_cluster::Int=1,
                              score_type::Symbol=:support)
    score = score_type === :support ?
        prey_support_score(pool, grid; τ=τ) :
        prey_vulnerability_score(pool, grid; τ=τ)

    dAcl_r = Float64[]; dInt_r = Float64[]; dSyn_r = Float64[]
    dAcl_c = Float64[]; dInt_c = Float64[]; dSyn_c = Float64[]
    dAcl_h = Float64[]; dInt_h = Float64[]; dSyn_h = Float64[]

    for f in losses
        keep = 1.0 - f
        km_r = random_mask(grid.C, keep; seed=101)
        km_c = clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=202)
        km_h = hotspot_clustered_mask_bestfirst(grid, keep; score=score,
                                                nseeds=nseeds_cluster, seed=303)
        for (km, A, I, S) in ((km_r, dAcl_r, dInt_r, dSyn_r),
                              (km_c, dAcl_c, dInt_c, dSyn_c),
                              (km_h, dAcl_h, dInt_h, dSyn_h))
            r = decomp_at_mask(pool, grid; τ=τ, keepmask=km)
            push!(A, r.agg.dAcl); push!(I, r.agg.dInt); push!(S, r.agg.dSyn)
        end
    end
    return (losses=losses,
            rand=(dAcl=dAcl_r, dInt=dInt_r, dSyn=dSyn_r),
            clust=(dAcl=dAcl_c, dInt=dInt_c, dSyn=dSyn_c),
            hotspot=(dAcl=dAcl_h, dInt=dInt_h, dSyn=dSyn_h))
end

"Plot three subpanels: Random / Clustered / Hotspot; each shows 3 lines."
function plot_decomp_three(series; title_str="ΔBSH decomposition vs loss (area basis)")
    L = series.losses
    begin
        fig = Figure(; size=(1080,360))
        for (j, key, name) in zip(1:3, (:rand,:clust,:hotspot), ["Random","Clustered","Hotspot"])
            ax = Axis(fig[1,j], title=name, xlabel="Area lost (fraction)",
                      ylabel = j==1 ? "ΔBSH decomposition" : "")
            X = L
            Y = getfield(series, key)
            lines!(ax, X, Y.dAcl, label="Climate-only")
            lines!(ax, X, Y.dInt, label="Interaction-only")
            lines!(ax, X, Y.dSyn, label="Synergy")
            hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
            if j==3; axislegend(ax, position=:lb); end
        end
        display(fig)
    end
end

pool = build_pool(220;
    basal_frac=0.35, seed=1,
    sigma=0.12, density=0.05, pmax=0.50,              # low redundancy
    niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.01,  # high synchrony
    b0_basal=0.05, bspread_basal=0.01,
    b0_cons=0.14, bspread_cons=0.05,
    R0_mean=3.0, R0_sd=0.10)

τ = 0.66
losses = collect(0.0:0.05:0.8)

series = decomp_series_3strat(pool, grid; τ=τ, losses=losses,
                              nseeds_cluster=1, score_type=:vulnerable)

plot_decomp_three(series; title_str="ΔBSH decomposition vs loss")

