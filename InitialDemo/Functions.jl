# --- GRID + CLIMATE ----------------------------------------------------------
struct Grid
    nx::Int                  # number of columns
    ny::Int                  # number of rows
    C::Int                   # total cells = nx*ny
    xy::Matrix{Float64}      # 2 x C matrix of (x,y) coords in [0,1]
    climate::Vector{Float64} # length-C vector in [0,1]
end

"""
make_grid(nx, ny; seed=42, texture=0.10)

Creates a regular nx-by-ny grid with a smooth climate field in [0,1].
`texture` controls small-scale variation on top of a diagonal gradient.
"""
function make_grid(nx::Int, ny::Int; seed::Int=42, texture::Float64=0.10)
    Random.seed!(seed)
    xs = range(0, 1; length=nx)
    ys = range(0, 1; length=ny)

    C = nx * ny
    xy = Matrix{Float64}(undef, 2, C)
    clim = Vector{Float64}(undef, C)

    k = 1
    for j in 1:ny, i in 1:nx
        x = xs[i]; y = ys[j]
        xy[:, k] .= (x, y)
        base    = 0.7*x + 0.3*y                     # diagonal gradient
        wiggle  = texture * sin(6π*x) * sin(6π*y)   # optional texture
        clim[k] = base + wiggle
        k += 1
    end

    # normalize climate to [0,1]
    clim .-= minimum(clim)
    clim ./= (maximum(clim) + eps())

    return Grid(nx, ny, C, xy, clim)
end

# Neighbor indices (4-neighborhood) for clustered masks
function neighbors4(ix::Int, nx::Int, ny::Int)
    i = ((ix - 1) % nx) + 1
    j = ((ix - 1) ÷ nx) + 1
    out = Int[]
    if i > 1;       push!(out, ix - 1);     end
    if i < nx;      push!(out, ix + 1);     end
    if j > 1;       push!(out, ix - nx);    end
    if j < ny;      push!(out, ix + nx);    end
    out
end

# --- SPECIES POOL + METAWEB --------------------------------------------------

struct SpeciesPool
    S::Int                      # number of species
    masses::Vector{Float64}     # body mass per species
    basal::BitVector            # basal flag per species (true = basal)
    mu::Vector{Float64}         # climatic niche center in [0,1]
    b::Vector{Float64}          # climatic niche breadth
    E::Vector{Vector{Int}}      # adjacency: E[s] = list of prey indices for predator s
end

"""
build_pool(S; kwargs...) -> SpeciesPool

Create a species pool and metaweb with **body-mass ratio** feeding rules.

Key knobs you’ll likely tune later:

  Redundancy (diet breadth / link density)
    - R0_mean::Float64 = 12.0     # preferred predator:prey mass ratio (≈ log-normal mean)
    - R0_sd::Float64   = 0.50     # SD (on log-scale) of preferred ratio among predators
    - sigma::Float64   = 0.50     # width of the log-ratio kernel (within predator)
    - density::Float64 = 0.30     # global thinning multiplier
    - pmax::Float64    = 0.90     # max acceptance probability at kernel peak

  Synchrony (climatic niches; affects spatial co-occurrence of basal prey)
    - niche_mode::Symbol = :uniform    # :uniform or :bimodal (basal only)
    - mu_basal_centers::Tuple = (0.25, 0.75)  # centers if bimodal
    - mu_basal_sd::Float64   = 0.05           # SD around each center (bimodal)
    - b0_basal::Float64 = 0.12    # baseline niche breadth for basal
    - b0_cons::Float64  = 0.12    # baseline niche breadth for consumers
    - bspread_basal::Float64 = 0.05
    - bspread_cons::Float64  = 0.05

Other:
    - basal_frac::Float64 = 0.45  # lightest fraction of species set as basal
    - seed::Int = 1               # RNG seed (thread-safe, local RNG)

Returns: SpeciesPool(S, masses, basal, mu, b, E)
"""
function build_pool(
    S::Int;
    basal_frac::Float64 = 0.45,
    seed::Int = 1,
    # diet / redundancy
    R0_mean::Float64 = 12.0,
    R0_sd::Float64   = 0.50,
    sigma::Float64   = 0.50,
    density::Float64 = 0.30,
    pmax::Float64    = 0.90,
    # climate / synchrony
    niche_mode::Symbol = :uniform,              # :uniform | :bimodal (basal only)
    mu_basal_centers::Tuple{Float64,Float64} = (0.25, 0.75),
    mu_basal_sd::Float64 = 0.05,
    b0_basal::Float64 = 0.12, bspread_basal::Float64 = 0.05,
    b0_cons::Float64  = 0.12, bspread_cons::Float64  = 0.05
)
    rng = MersenneTwister(seed)

    # body masses (log-spread)
    logm   = collect(range(log(1e-2), log(10.0); length=S))
    shuffle!(rng, logm)
    masses = exp.(logm)

    # basal = LIGHTEST fraction
    order = sortperm(masses)               # increasing mass
    nB    = round(Int, basal_frac * S)
    basal = falses(S); basal[order[1:nB]] .= true

    # climatic niches
    mu = similar(masses)
    b  = similar(masses)

    # consumers first (uniform)
    cons_ids = findall(!, basal)
    mu[cons_ids] .= rand(rng, length(cons_ids))
    b[cons_ids]  .= b0_cons .+ bspread_cons .* rand(rng, length(cons_ids))

    # basal (uniform or bimodal)
    bas_ids = findall(basal)
    if niche_mode === :bimodal
        nb  = length(bas_ids)
        nb1 = nb ÷ 2; nb2 = nb - nb1
        c1, c2 = mu_basal_centers
        mu[bas_ids[1:nb1]]     .= clamp.(c1 .+ mu_basal_sd .* randn(rng, nb1), 0, 1)
        mu[bas_ids[nb1+1:end]] .= clamp.(c2 .+ mu_basal_sd .* randn(rng, nb2), 0, 1)
    else
        mu[bas_ids] .= rand(rng, length(bas_ids))
    end
    b[bas_ids] .= b0_basal .+ bspread_basal .* rand(rng, length(bas_ids))

    # predator-specific preferred ratio (log-normal)
    R0 = exp.(log(R0_mean) .+ R0_sd .* randn(rng, S))

    # metaweb: probabilistic acceptance on log mass ratio
    E = [Int[] for _ in 1:S]
    for (ii, s) in pairs(order)
        basal[s] && continue
        for jj in 1:ii-1
            q = order[jj]
            r = masses[s] / masses[q]                    # ratio > 1
            z = (log(r) - log(R0[s])) / sigma
            p = pmax * exp(-0.5 * z^2) * density        # acceptance prob
            if rand(rng) < p
                push!(E[s], q)
            end
        end
        # ensure every consumer has ≥1 potential prey (fallback = nearest ratio)
        if isempty(E[s]) && ii > 1
            cand   = order[1:ii-1]
            target = log(masses[s]) - log(R0[s])
            qstar  = cand[argmin(abs.(log.(masses[cand]) .- target))]
            push!(E[s], qstar)
        end
    end

    SpeciesPool(S, masses, basal, mu, b, E)
end

# ---------- Climate pass (Z) ----------
# Suitability (Gaussian in climate); pass if ≥ τ
function climate_pass(pool::SpeciesPool, grid::Grid; τ=0.5)
    S, C = pool.S, grid.C
    Z = falses(S, C)
    @inbounds for s in 1:S, c in 1:C
        su = exp(-((grid.climate[c] - pool.mu[s])^2) / (2 * pool.b[s]^2 + 1e-12))
        if su ≥ τ
            Z[s,c] = true
        end
    end
    Z
end

# ---------- Assembly: prune unsupported consumers until fixed point ----------
function assemble(Z::BitMatrix, pool::SpeciesPool)
    S, C = size(Z)
    P = copy(Z) # start with climate pass
    changed = true
    while changed
        changed = false
        @inbounds for c in 1:C, s in 1:S
            if P[s,c] && !pool.basal[s]
                prey = pool.E[s]
                has_pre = false
                @inbounds for q in prey
                    if P[q,c]; has_pre = true; break; end
                end
                if !has_pre
                    P[s,c] = false
                    changed = true
                end
            end
        end
    end
    P
end

# ---------- Metrics ----------
# Biotically Safe Climate Habitat (BSH, r=1) per species (fraction of ALL cells)
function bsh1_per_species(P::BitMatrix, Z::BitMatrix, pool::SpeciesPool)
    S, C = size(P)
    if C == 0
        return zeros(Float64, S)   # or fill(0.0, S) — define your convention
    end
    out = zeros(Float64, S)
    @inbounds for s in 1:S
        if isempty(pool.E[s])              # basal
            out[s] = count(@view Z[s, :]) / C
        else
            ok = 0
            for c in 1:C
                if Z[s,c]
                    for q in pool.E[s]
                        if P[q,c]; ok += 1; break; end
                    end
                end
            end
            out[s] = ok / C
        end
    end
    out
end

# Count of BSH cells for each species (no normalization)
function bsh1_count(P::BitMatrix, Z::BitMatrix, pool::SpeciesPool)
    S, C = size(P)
    out = zeros(Int, S)
    @inbounds for s in 1:S
        if pool.basal[s]               # << use basal flag (not isempty(E))
            out[s] = count(@view Z[s, :])
        else
            cnt = 0
            for c in 1:C
                if Z[s,c]
                    for q in pool.E[s]
                        if P[q,c]; cnt += 1; break; end
                    end
                end
            end
            out[s] = cnt
        end
    end
    return out
end

# Fraction of ORIGINAL area that is BSH after masking
# (divide the BSH counts computed on the masked grid by the ORIGINAL C_full)
function bsh1_area_fraction(P::BitMatrix, Z::BitMatrix, pool::SpeciesPool, C_full::Int)
    counts = bsh1_count(P, Z, pool)
    return counts ./ C_full
end

# ILRF per species: share of climate-suitable cells where prey are absent
function ilrf_per_species(P::BitMatrix, Z::BitMatrix, pool::SpeciesPool)
    S, C = size(P)
    out = fill(Float64(NaN), S)
    @inbounds for s in 1:S
        suitable = 0
        no_prey  = 0
        if pool.basal[s]
            out[s] = 0.0
        else
            for c in 1:C
                if Z[s,c]
                    suitable += 1
                    prey_ct = 0
                    for q in pool.E[s]
                        prey_ct += P[q,c] ? 1 : 0
                        if prey_ct > 0; break; end
                    end
                    if prey_ct == 0; no_prey += 1; end
                end
            end
            out[s] = (suitable == 0) ? NaN : no_prey / suitable
        end
    end
    out
end

# Mean over consumers only
consumer_mask(pool::SpeciesPool) = .!pool.basal

# ---------- Habitat-loss masks ----------
# Apply a keep-mask (true = keep cell; false = removed)
apply_mask(Z::BitMatrix, keep::BitVector) = Z[:, keep]

# Random keep mask for given keep_frac
function random_mask(C::Int, keep_frac::Float64; seed::Int=0)
    Random.seed!(seed)
    keep = falses(C)
    nkeep = round(Int, keep_frac * C)
    keep[randperm(C)[1:nkeep]] .= true
    keep
end

# Clustered removal via BFS growth from multiple seeds until target reached
function clustered_mask(grid::Grid, keep_frac::Float64; nseeds::Int=6, seed::Int=0)
    Random.seed!(seed)
    C = grid.C
    target_remove = C - round(Int, keep_frac*C)
    removed = falses(C)
    q = Int[]
    # pick initial seeds
    seeds = randperm(C)[1:nseeds]
    append!(q, seeds)
    removed_ct = 0
    ptr = 1
    while removed_ct < target_remove
        if ptr > length(q)
            # start a new seed in a not-yet-removed cell
            newseed = first(filter(i -> !removed[i], randperm(C)))
            push!(q, newseed)
        end
        v = q[ptr]; ptr += 1
        if removed[v]; continue; end
        removed[v] = true
        removed_ct += 1
        # enqueue neighbors
        for nb in neighbors4(v, grid.nx, grid.ny)
            if !removed[nb]
                push!(q, nb)
            end
        end
    end
    .!removed
end

# BSH conditional on climate-suitable cells (interaction lens)
function bsh1_cond_per_species(P::BitMatrix, Z::BitMatrix, pool::SpeciesPool)
    S, C = size(P)
    out = fill(Float64(NaN), S)
    @inbounds for s in 1:S
        suitable = findall(@view Z[s, :])
        denom = length(suitable)
        if denom == 0; continue; end
        if isempty(pool.E[s]) # basal
            out[s] = 1.0  # within suitable climate, basal are always 'effective'
        else
            ok = 0
            for c in suitable
                prey_ct = 0
                for q in pool.E[s]
                    prey_ct += P[q,c] ? 1 : 0
                    if prey_ct > 0; break; end
                end
                ok += (prey_ct > 0) ? 1 : 0
            end
            out[s] = ok / denom
        end
    end
    out
end

# Climate-suitable fraction per species (A_s)
function climate_area_per_species(Z::BitMatrix)
    S, C = size(Z)
    [count(@view Z[s, :]) / C for s in 1:S]
end

# Midpoint (Shapley-style) decomposition for each species:
# B = A * p  (A = climate area; p = conditional BSH)
# ΔB = p̄ * ΔA + Ā * Δp + ΔA * Δp, where bars are midpoints
struct Decomp{T}
    dB::T; dA_only::T; dInt_only::T; synergy::T
end

function decompose_delta(B0, A0, p0, B1, A1, p1)
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

struct Summary
    loss::Vector{Float64}
    ΔBSH_mean::Vector{Float64}
    ΔBSH_lo::Vector{Float64}
    ΔBSH_hi::Vector{Float64}
end

"""
sweep_ensemble(pool_seed_list, mask_seeds_list; kwargs...) -> Summary

Run a loss sweep and return ΔBSH envelopes (mean, 10%, 90%) across pool × mask
replicates. This function is *agnostic* to redundancy/synchrony; pass those via
`build_pool`’s keyword arguments.

Arguments (most important):
  - kind::Symbol = :random        # :random | :clustered
  - grid::Grid                    # your grid
  - τ::Float64 = 0.5              # climate threshold
  - loss_fracs = 0.0:0.05:0.9     # fraction of area removed
  - S::Int, basal_frac::Float64   # pool size and basal share
  - nseeds_cluster::Int = 6       # number of clusters for :clustered
  - metric::Symbol = :fraction    # :fraction (on kept cells) or :area (vs original area)
  - pool_kwargs...                # **forwarded to build_pool** (e.g., sigma, density, niche_mode,…)

Returns: Summary(loss, ΔBSH_mean, ΔBSH_lo, ΔBSH_hi)
"""
function sweep_ensemble(pool_seed_list, mask_seeds_list;
        kind::Symbol = :random,
        grid,
        τ::Float64 = 0.5,
        loss_fracs = 0.0:0.05:0.9,
        S::Int,
        basal_frac::Float64,
        nseeds_cluster::Int = 6,
        metric::Symbol = :fraction,
        pool_kwargs...
)
    n    = length(loss_fracs)
    μ    = fill(NaN, n)
    lo   = fill(NaN, n)
    hi   = fill(NaN, n)
    Cfull = grid.C

    Threads.@threads for k in 1:n
        f    = loss_fracs[k]
        keep = 1.0 - f
        vals = Float64[]

        for ps in pool_seed_list, ms in mask_seeds_list
            # Build pool with whatever parametrisation you pass via pool_kwargs
            pool  = build_pool(S; basal_frac=basal_frac, seed=ps, pool_kwargs...)
            Zfull = climate_pass(pool, grid; τ=τ)
            baseP = assemble(Zfull, pool)

            base = metric === :area ?
                mean(bsh1_area_fraction(baseP, Zfull, pool, Cfull)[.!pool.basal]) :
                mean(bsh1_per_species(baseP, Zfull, pool)[.!pool.basal])

            # keepmask = kind === :random ?
            #     random_mask(Cfull, keep; seed=ms) :
            #     clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=ms)
            # # inside your sweep_ensemble inner loop, where you choose a keepmask:

            keepmask = begin
                if kind === :random
                    random_mask(grid.C, keep; seed=ms)
                elseif kind === :clustered
                    clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=ms)
                elseif kind === :hotspot_cons
                    # NOTE: depends on the pool, so build pool first in the loop
                    consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep, power=hotspot_power, seed=ms)
                elseif kind === :hotspot_prey
                    prey_hotspot_mask(grid, pool; τ=τ, keep_frac=keep, power=hotspot_power, seed=ms)
                else
                    error("Unknown kind = $kind")
                end
            end

            Z = apply_mask(Zfull, keepmask)
            if size(Z, 2) == 0
                push!(vals, 0.0); continue
            end
            P = assemble(Z, pool)

            b = metric === :area ?
                mean(bsh1_area_fraction(P, Z, pool, Cfull)[.!pool.basal]) :
                mean(bsh1_per_species(P, Z, pool)[.!pool.basal])

            push!(vals, b - base)
        end

        μ[k]  = mean(vals)
        lo[k] = quantile(vals, 0.10)
        hi[k] = quantile(vals, 0.90)
    end

    Summary(collect(loss_fracs), μ, lo, hi)
end

# Redundancy: mean #prey present per suitable cell (baseline)
function redundancy_per_consumer(P, Z, pool)
    S, C = size(P)
    out = fill(NaN, S)
    for s in 1:S
        pool.basal[s] && continue
        idx = findall(@view Z[s, :])
        if !isempty(idx)
            cnt = 0; denom = 0
            for c in idx
                denom += 1
                for q in pool.E[s]
                    if P[q,c]; cnt += 1; end
                end
            end
            out[s] = cnt/denom
        end
    end
    out
end

# Synchrony: average pairwise correlation of prey presences across cells
function prey_synchrony_per_consumer(P, pool)
    S, C = size(P)
    out = fill(NaN, S)
    for s in 1:S
        pool.basal[s] && continue
        prey = pool.E[s]
        np = length(prey)
        if np ≥ 2
            M = BitMatrix(P[prey, :])
            v = Vector{Float64}()
            for i in 1:np-1, j in i+1:np
                c = cor(Float64.(M[i, :]), Float64.(M[j, :]))
                if isfinite(c); push!(v, c); end
            end
            if !isempty(v); out[s] = mean(v); end
        end
    end
    out
end
