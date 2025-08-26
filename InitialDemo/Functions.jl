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
build_pool(S; basal_frac=0.45, seed=1,
           Rmin=3.0, Rmax=300.0,   # mass-ratio diet window
           b0=0.12, bspread=0.08)  # niche breadths

Creates S species with:
- log-spread body masses,
- the LIGHTEST ~basal_frac fraction set as basal,
- random climatic niches (mu, breadth b),
- a metaweb where a predator s eats prey q if mass ratio m_s/m_q ∈ [Rmin, Rmax].
Guarantees every consumer has ≥1 potential prey (fallback to immediate smaller).
"""
# Replace your build_pool with this variant
function build_pool(S::Int=140; basal_frac=0.5, seed=1,
                    R0_mean=16.0, R0_sd=0.6,   # predator-specific preferred ratio
                    sigma=0.35,                # spread of log-ratio kernel
                    density=0.05, pmax=0.7,    # overall thinning
                    b0=0.12, bspread=0.08)

    Random.seed!(seed)

    # masses (log-spread)
    logm   = collect(range(log(1e-2), log(10.0); length=S))
    shuffle!(logm)
    masses = exp.(logm)

    # basal = LIGHTEST species
    order = sortperm(masses)     # increasing mass
    nB    = round(Int, basal_frac*S)
    basal = falses(S); basal[order[1:nB]] .= true

    # climate niches
    mu = rand(S)
    b  = b0 .+ bspread .* rand(S)

    # predator-specific preferred ratio
    R0 = exp.(log(R0_mean) .+ R0_sd .* randn(S))

    # build metaweb with probabilistic kernel on log ratio
    E = [Int[] for _ in 1:S]
    for (ii, s) in pairs(order)
        basal[s] && continue
        for jj in 1:ii-1
            q = order[jj]
            r = masses[s] / masses[q]              # ratio > 1
            z = (log(r) - log(R0[s])) / sigma
            p = pmax * exp(-0.5*z^2) * density     # acceptance prob
            if rand() < p
                push!(E[s], q)
            end
        end
        # guarantee ≥1 potential prey (pick nearest in log-ratio)
        if isempty(E[s]) && ii > 1
            cand = order[1:ii-1]
            target = log(masses[s]) - log(R0[s])
            qstar = cand[argmin(abs.(log.(masses[cand]) .- target))]
            push!(E[s], qstar)
        end
    end
    return SpeciesPool(S, masses, basal, mu, b, E)
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
    dB::T; dA_part::T; dp_part::T; synergy::T
end

function decompose_delta(B0, A0, p0, B1, A1, p1)
    dB = B1 - B0
    dA = A1 - A0
    dp = p1 - p0
    Abar = 0.5*(A0 + A1)
    pbar = 0.5*(p0 + p1)
    dA_part = pbar * dA
    dp_part = Abar * dp
    synergy = dB - dA_part - dp_part
    Decomp(dB, dA_part, dp_part, synergy)
end

struct Summary
    loss::Vector{Float64}
    ΔBSH_mean::Vector{Float64}
    ΔBSH_lo::Vector{Float64}
    ΔBSH_hi::Vector{Float64}
end

function sweep_ensemble(pool_seed_list, mask_seeds_list;
        kind::Symbol = :random,
        grid,
        τ::Float64 = 0.5,
        loss_fracs = 0.0:0.05:0.9,
        S::Int,
        basal_frac::Float64,
        nseeds_cluster::Int = 6,
        metric::Symbol = :area)   # :area or :fraction

    n = length(loss_fracs)
    μ  = fill(NaN, n); lo = fill(NaN, n); hi = fill(NaN, n)
    Cfull = grid.C

    Threads.@threads for k in 1:n
        f    = loss_fracs[k]
        keep = 1.0 - f
        vals = Float64[]

        for ps in pool_seed_list, ms in mask_seeds_list
            pool = build_pool(
                200;
                basal_frac = 0.35,        # not too low (need prey biomass), not too high (consumers still picky)
                # diet kernel: narrow + focused → few prey per consumer
                R0_mean = 10.0,
                R0_sd   = 0.25,           # predators prefer similar ratios → less flexibility
                sigma   = 0.22,           # narrow diet window → low redundancy
                density = 0.12,           # thin metaweb
                pmax    = 0.70,
                # climate niches: basal narrower & clustered → prey synchrony
                b0      = 0.10,           # this is the *consumer* base; we’ll override basal below
                bspread = 0.04
            )
            # prey synchrony: split basal into two tight climate guilds
            basal_ids = findall(pool.basal)
            nb = length(basal_ids)
            nb1 = nb ÷ 2
            nb2 = nb - nb1

            pool.mu[basal_ids[1:nb1]] .= clamp.(0.25 .+ 0.03 .* randn(nb1), 0, 1)
            pool.mu[basal_ids[nb1+1:end]] .= clamp.(0.75 .+ 0.03 .* randn(nb2), 0, 1)
            pool.b[basal_ids] .= 0.06 .+ 0.01 .* rand(nb)    # basal narrower
            # consumers a bit broader (leave pool.b for non-basal as created)
                        Zfull = climate_pass(pool, grid; τ=τ)
            baseP = assemble(Zfull, pool)

            base = if metric === :area
                mean(bsh1_area_fraction(baseP, Zfull, pool, Cfull)[.!pool.basal])
            else
                mean(bsh1_per_species(baseP, Zfull, pool)[.!pool.basal])  # fraction on current cells
            end

            keepmask = kind === :random ?
                random_mask(Cfull, keep; seed=ms) :
                clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=ms)

            Z = apply_mask(Zfull, keepmask)
            P = assemble(Z, pool)

            b = if metric === :area
                mean(bsh1_area_fraction(P, Z, pool, Cfull)[.!pool.basal])
            else
                mean(bsh1_per_species(P, Z, pool)[.!pool.basal])
            end

            push!(vals, b - base)
            Z = climate_pass(pool, grid; τ=0.5)

            th = trophic_height_per_cell(P, pool)
            @show maximum(th), count(>(1), th)/length(th)   # want max ≥ 3 and many cells >1
        end

        μ[k]  = mean(vals)
        lo[k] = quantile(vals, 0.10)
        hi[k] = quantile(vals, 0.90)
    end

    return Summary(collect(loss_fracs), μ, lo, hi)
end

