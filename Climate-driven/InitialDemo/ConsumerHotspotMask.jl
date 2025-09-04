# ---- 4-neighbourhood helper (only if you don't already have it) ----
if !@isdefined(neighbors4)
function neighbors4(idx::Int, nx::Int, ny::Int)
    i = ((idx-1) % nx) + 1
    j = ((idx-1) ÷ nx) + 1
    nb = Int[]
    if i>1  push!(nb, idx-1)   end
    if i<nx push!(nb, idx+1)   end
    if j>1  push!(nb, idx-nx)  end
    if j<ny push!(nb, idx+nx)  end
    return nb
end
end

# ---- consumer "relevance" score per cell, from the baseline assembly ----
"""
consumer_hotspot_scores(pool, grid; τ=0.5, θ=0.0)

Returns a length-C vector `score` where `score[c]` counts how many **consumers**
would treat cell `c` as contributing to their viable habitat in the **baseline**:
- The cell must be climate-suitable for that consumer, and
- It must contain at least one prey (θ=0) or a fraction ≥ θ of the consumer's potential prey.

Use a higher θ (e.g. 0.5–0.7) to make "relevance" stricter.
"""
function consumer_hotspot_scores(pool::SpeciesPool, grid::Grid; τ::Float64=0.5, θ::Float64=0.0)
    Z0 = climate_pass(pool, grid; τ=τ)
    P0 = assemble(Z0, pool)
    S, C = size(P0)
    cons = .!pool.basal
    score = zeros(Float64, C)

    @inbounds for c in 1:C
        cnt = 0
        for s in 1:S
            cons[s] || continue
            if Z0[s,c]
                prey = pool.E[s]
                isempty(prey) && continue
                present = 0
                for q in prey
                    present += (P0[q,c] ? 1 : 0)
                    if θ == 0.0 && present>0
                        break
                    end
                end
                ok = (θ == 0.0) ? (present>0) : ((present/length(prey)) >= θ)
                cnt += ok ? 1 : 0
            end
        end
        score[c] = cnt
    end
    return score
end

# ---- weighted choice without external deps ----
@inline function _weighted_pick(idxs::Vector{Int}, w::Vector{Float64})
    s = 0.0
    tot = 0.0
    @inbounds for i in idxs
        tot += w[i]
    end
    u = rand() * (tot + eps())
    @inbounds for i in idxs
        s += w[i]
        if s >= u
            return i
        end
    end
    return last(idxs)
end

# ---- hotspot-clustered mask (contiguous growth, biased to high-score cells) ----
"""
hotspot_clustered_mask(grid, keep_frac; scores, nseeds=1, seed=0, power=2.0)

Removes a contiguous blob until the target number of cells is removed.
Growth is **biased** toward cells with larger `scores` (from `consumer_hotspot_scores`),
using weights `scores.^power`.

Returns a BitVector `keep` (true = kept).
"""
function hotspot_clustered_mask(grid::Grid, keep_frac::Float64;
                                scores::Vector{Float64},
                                nseeds::Int=1, seed::Int=0, power::Float64=2.0)
    Random.seed!(seed)
    C  = grid.C
    nremove = C - round(Int, keep_frac*C)
    removed = falses(C)

    # positive weights for bias
    w = (scores .+ 1e-9) .^ power

    # pick initial seeds (weighted by w among not-removed)
    seeds = Int[]
    avail = collect(1:C)
    for _ in 1:nseeds
        s = _weighted_pick(avail, w)
        push!(seeds, s)
        filter!(x -> x != s, avail)
    end

    frontier = copy(seeds)
    ptr = 1
    removed_ct = 0

    while removed_ct < nremove
        if ptr > length(frontier)
            # pick a new weighted seed in any not-removed location
            avail = [i for i in 1:C if !removed[i]]
            isempty(avail) && break
            newseed = _weighted_pick(avail, w)
            push!(frontier, newseed)
        end

        # pick next frontier cell **with weight bias** among frontier
        active_front = [i for i in frontier[ptr:end] if !removed[i]]
        if isempty(active_front)
            ptr = length(frontier) + 1
            continue
        end
        v = _weighted_pick(active_front, w)
        ptr = findfirst(==(v), frontier) + 1  # move pointer past chosen element

        if removed[v]; continue; end
        removed[v] = true
        removed_ct += 1

        # push 4-neighbors
        for nb in neighbors4(v, grid.nx, grid.ny)
            if !removed[nb]
                push!(frontier, nb)
            end
        end
    end

    return .!removed
end

# ---- small helper: one-shot decomposition at a given keep mask and basis ----
"""
decomp_at_mask(pool, grid; τ, keepmask, basis=:area)

Returns a NamedTuple with the **mean over consumers** of:
- total ΔB (mean dB),
- climate-only part (mean p̄ ΔA),
- interaction-only part (mean Ā Δp),
- synergy (mean ΔA Δp).

`basis=:fraction` -> A is relative to **original** area.
`basis=:area`     -> both baseline and scenario are evaluated on the **kept** area.
"""
function decomp_at_mask(pool::SpeciesPool, grid::Grid;
                        τ::Float64=0.5, keepmask::BitVector,
                        basis::Symbol=:area)
    # baseline
    Z0_full = climate_pass(pool, grid; τ=τ)
    P0_full = assemble(Z0_full, pool)

    # basis choice: evaluate baseline on kept area for :area basis
    if basis == :area
        Z0 = apply_mask(Z0_full, keepmask)
        P0 = assemble(Z0, pool)
    else
        Z0 = Z0_full
        P0 = P0_full
    end

    # scenario
    Z1 = apply_mask(Z0_full, keepmask)
    P1 = assemble(Z1, pool)

    cons = .!pool.basal
    A0 = climate_area_per_species(Z0)[cons]
    p0 = bsh1_cond_per_species(P0, Z0, pool)[cons]
    B0 = A0 .* p0

    A1 = climate_area_per_species(Z1)[cons]
    p1 = bsh1_cond_per_species(P1, Z1, pool)[cons]
    B1 = A1 .* p1

    dA = A1 .- A0
    dp = p1 .- p0
    dB = B1 .- B0
    Abar = 0.5 .* (A0 .+ A1)
    pbar = 0.5 .* (p0 .+ p1)

    climate = mean(pbar .* dA)
    inter   = mean(Abar .* dp)
    synergy = mean(dA .* dp)
    total   = mean(dB)

    return (total=total, climate=climate, interaction=inter, synergy=synergy)
end

# ---- run curves across loss for three strategies; return BOTH bases ----
"""
run_decomp_curves(; grid, τ, S, basal_frac, losses, pool_seeds, mask_seeds,
                    nseeds_cluster=1, hotspot_power=3.0, hotspot_θ=0.0,
                    pool_kwargs...)

Builds a pool for each `pool_seed`, computes consumer-hotspot **scores** once,
then for each loss `f` and each `mask_seed` computes RANDOM, CLUSTERED,
and HOTSPOT masks and returns mean curves for **both** bases.

Returns:
    curves_area::Dict{Symbol,NamedTuple}
    curves_frac::Dict{Symbol,NamedTuple}
Each NamedTuple has fields (:loss, :climate, :interaction, :synergy).
"""
function run_decomp_curves(; grid::Grid, τ::Float64=0.5, S::Int=200, basal_frac::Float64=0.35,
                           losses = 0.0:0.05:0.8,
                           pool_seeds = 1:3, mask_seeds = 1:10,
                           nseeds_cluster::Int=1, hotspot_power::Float64=3.0,
                           hotspot_θ::Float64=0.0,
                           pool_kwargs...)
    kinds = (:random, :clustered, :hotspot)

    data_area  = Dict(k => (climate=Float64[], interaction=Float64[], synergy=Float64[], loss=Float64[]) for k in kinds)
    data_frac  = Dict(k => (climate=Float64[], interaction=Float64[], synergy=Float64[], loss=Float64[]) for k in kinds)

    for f in losses
        keep = 1.0 - f

        # accumulators per kind & basis
        acc_area  = Dict(k => (Float64[], Float64[], Float64[]) for k in kinds)   # (clim, inter, syn)
        acc_frac  = Dict(k => (Float64[], Float64[], Float64[]) for k in kinds)

        for ps in pool_seeds
            pool = build_pool(S; basal_frac=basal_frac, seed=ps, pool_kwargs...)
            # hotspot scores from the **baseline** for this pool
            scores = consumer_hotspot_scores(pool, grid; τ=τ, θ=hotspot_θ)

            for ms in mask_seeds
                km_random   = random_mask(grid.C, keep; seed=ms)
                km_cluster  = clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=ms)
                km_hotspot  = hotspot_clustered_mask(grid, keep;
                                scores=scores, nseeds=nseeds_cluster, seed=ms, power=hotspot_power)

                for (k, km) in zip(kinds, (km_random, km_cluster, km_hotspot))
                    # area basis
                    dA = decomp_at_mask(pool, grid; τ=τ, keepmask=km, basis=:area)
                    push!(acc_area[k][1], dA.climate)
                    push!(acc_area[k][2], dA.interaction)
                    push!(acc_area[k][3], dA.synergy)
                    # fraction basis
                    dF = decomp_at_mask(pool, grid; τ=τ, keepmask=km, basis=:fraction)
                    push!(acc_frac[k][1], dF.climate)
                    push!(acc_frac[k][2], dF.interaction)
                    push!(acc_frac[k][3], dF.synergy)
                end
            end
        end

        # write means into the curve containers
        for k in kinds
            climA = mean(acc_area[k][1]);  interA = mean(acc_area[k][2]);  synA = mean(acc_area[k][3])
            push!(data_area[k].climate, climA);  push!(data_area[k].interaction, interA);  push!(data_area[k].synergy, synA); push!(data_area[k].loss, f)
            climF = mean(acc_frac[k][1]);  interF = mean(acc_frac[k][2]);  synF = mean(acc_frac[k][3])
            push!(data_frac[k].climate, climF);  push!(data_frac[k].interaction, interF);  push!(data_frac[k].synergy, synF); push!(data_frac[k].loss, f)
        end
    end

    # pack as NamedTuples
    curves_area = Dict(k => (; loss = data_area[k].loss,
                               climate = data_area[k].climate,
                               interaction = data_area[k].interaction,
                               synergy = data_area[k].synergy) for k in kinds)
    curves_frac = Dict(k => (; loss = data_frac[k].loss,
                               climate = data_frac[k].climate,
                               interaction = data_frac[k].interaction,
                               synergy = data_frac[k].synergy) for k in kinds)
    return curves_area, curves_frac
end

# --- set your grid & pool knobs (these can be anything; below pushes contrast) ---
grid = make_grid(60, 60; seed=11)

curves_area, curves_frac = run_decomp_curves(
    grid = grid, τ = 0.55,
    S = 200, basal_frac = 0.35,
    losses = 0.0:0.05:0.8,
    pool_seeds = 1:4, mask_seeds = 1:10,
    nseeds_cluster = 1,
    hotspot_power = 4.0,     # ↑ makes the hotspot blob home in on top consumer cells
    hotspot_θ = 0.4;         # ↑ stricter "relevance" for consumers (optional)
    # you can pass diet/niche params of build_pool here, e.g.:
    # sigma=0.25, density=0.15, pmax=0.7, niche_mode=:bimodal, mu_basal_sd=0.04
)

# ---- plotting helper ----
function plot_decomp_grid(curves; title_left::String)
    kinds = (:random, :clustered, :hotspot)
    titles = Dict(:random=>"Random", :clustered=>"Clustered", :hotspot=>"Hotspot")
    fig = Figure(resolution=(1150, 380))
    for (j,k) in enumerate(kinds)
        ax = Axis(fig[1,j], title=titles[k], xlabel="Area lost (fraction)", ylabel = (j==1 ? "ΔBSH decomposition" : ""))
        lines!(ax, curves[k].loss, curves[k].climate,     label="Climate-only")
        lines!(ax, curves[k].loss, curves[k].interaction, label="Interaction-only")
        lines!(ax, curves[k].loss, curves[k].synergy,     label="Synergy")
        hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
        j==3 && axislegend(ax, position=:lt)
    end
    fig[0,1:3] = Label(fig, title_left, halign = :left, fontsize=16, padding=(0,0,6,0))
    display(fig)
    # return fig
end

figA = plot_decomp_grid(curves_area; title_left="AREA basis (evaluate baseline on kept cells)")
figF = plot_decomp_grid(curves_frac; title_left="FRACTION basis (vs original area)")

"""
consumer_hotspot_mask(grid, pool; τ, keep_frac, power=2.0, what=:cons, seed=0)

Returns a Bool vector of length grid.C where `true` means **keep** the cell.
Cells are ranked on the *baseline* community (full grid) using a hotspot score:

  what = :cons    -> number of consumer species present in the cell
  what = :basal   -> number of basal species that are climate-suitable in the cell
  what = :hybrid  -> 0.5*consumers_present + 0.5*basal_suitable

`power` (>1) sharpens the ranking (emphasizes peaks). `keep_frac` is the
fraction of cells to keep. Ties are broken with tiny random jitter (seeded).
"""
function consumer_hotspot_mask(grid::Grid, pool::SpeciesPool;
        τ::Float64, keep_frac::Float64, power::Float64=2.0,
        what::Symbol=:cons, seed::Int=0)

    # Baseline (full grid)
    Z0 = climate_pass(pool, grid; τ=τ)     # climate pass (S×C BitMatrix)
    P0 = assemble(Z0, pool)                # assembled presence (S×C BitMatrix)

    basal = pool.basal

    # Scores per cell
    cons_present = Float64.(vec(sum(@view P0[.!basal, :]; dims=1)))
    basal_ok     = Float64.(vec(sum(@view Z0[basal,  :]; dims=1)))

    score = if what === :cons
        cons_present
    elseif what === :basal
        basal_ok
    elseif what === :hybrid
        0.5 .* cons_present .+ 0.5 .* basal_ok
    else
        error("unknown `what` = $what (use :cons, :basal, or :hybrid)")
    end

    # Sharpen & rank
    score .= (score .+ 1e-12) .^ power
    Random.seed!(seed)
    jitter = 1e-9 .* rand(length(score))
    ord = sortperm(score .+ jitter; rev=true)

    # Keep top fraction
    C = grid.C
    nkeep = clamp(round(Int, keep_frac * C), 0, C)
    keep = falses(C)
    keep[ord[1:nkeep]] .= true
    return keep
end
