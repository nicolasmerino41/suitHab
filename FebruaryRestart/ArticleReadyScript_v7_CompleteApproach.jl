#!/usr/bin/env julia
# Julia 1.11
# Run:  julia --threads auto biotic_divergence_pipeline.jl
#
# PURPOSE (self-contained, no external data)
# ------------------------------------------------------------
# A vs AB baseline divergence (no habitat loss):
# - A: climatic suitability only
# - AB: climatic suitability + one-prey trophic rule (fixed point)
# - Optional spatial viability filter: species persists only if LCC ≥ Emin (applied AFTER A/AB maps are built)
#
# Sweeps:
# - Network family (4): Random (ER-like), Modular (SBM), Heavy-tailed out-degree, Cascade hierarchy
# - Connectance C (x-axis)
# - Mechanistic niche correlation r (y-axis): consumer μ correlated with mean prey μ
# - Niche breadth regime (4): Narrow/Broad × LowVar/HighVar
# - Environmental space (2): Random vs Autocorrelated field
#
# Outputs:
# - Heatmaps (per env type) for 3 divergence metrics:
#   (i) Relative richness loss: 1 - S_AB / S_A   (consumers-only)
#   (ii) Mean per-species Jaccard mismatch: mean_i [1 - J(A_i, AB_i)]  (consumers-only, |A_i|>0)
#   (iii) Fraction affected: frac_i [A_i != AB_i]  (consumers-only, |A_i|>0)
# - Plus an extra heatmap figure for realized overlap (diagnostic/mediator).
#
# Notes:
# - "One-prey rule" is hard: consumer present in a cell iff it is climatically suitable AND at least one prey is present there.
# - Basal species have no prey and are not trophically restricted (AB = A).
# - Connectance here follows the previous convention: Ltarget ≈ round(C * S^2),
#   while links are only from consumers -> prey, so actual realized connectance is reported.

using Random
using Statistics
using LinearAlgebra
using CairoMakie
using Printf
using Dates

# ============================================================
# 0) Global parameters (tune here)
# ============================================================
# Spatial grid
const NX = 60
const NY = 60
const NCELLS = NX * NY

# Species pool
const S = 250
const BASAL_FRAC = 0.30  # basal species fraction

# Spatial viability filter (movement/connectivity proxy)
const USE_CONNECTIVITY_FILTER = true
const Emin_patch = 60  # LCC threshold for persistence (set smaller if you want less stringent)

# Environmental field domain (e.g., temperature)
const E_MIN = 0.0
const E_MAX = 100.0

# Niche suitability: Gaussian with threshold
# suitability = exp(-0.5 * ((E-μ)/σ)^2 ) >= SUIT_THRESH
const SUIT_THRESH = 0.25

# Environmental autocorrelation
const AUTOCORR_ITERS = 18
const AUTOCORR_ALPHA = 0.55  # 0..1 (higher = smoother)

# Sweep axes
const CONNECTANCE_RANGE = (0.01, 0.1)
const CORR_RANGE       = (0.0, 0.9)
const N_CONNECT = 12
const N_CORR    = 12

# Replicates per heatmap cell (increase for final)
const NREP = 8

# Network-family knobs
const N_MODULES = 6
const MODULAR_IN_BIAS = 6.0        # >1 increases within-module links vs between
const HEAVYTAIL_GAMMA = 2.2        # out-degree heaviness
const HEAVYTAIL_KMAX_FRAC = 0.35   # kmax = round(frac*(S-1))
const CASCADE_LAMBDA = 2.5         # 0 = uniform among lower ranks; higher = interval-like

# Mechanistic niche-correlation builder knobs
const RELAX_ITERS = 30
const MU_NOISE_SD = 1.8            # consumer μ jitter during relaxation
const TARGET_R_TOL = 0.03

# Thread-safe seeds
const BASE_SEED = 20260202

# Output directory
# ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
OUTDIR = joinpath(pwd(), "output_biotic_divergence_")
isdir(OUTDIR) || mkpath(OUTDIR)

# ============================================================
# 1) Index helpers + neighbors (4-neighbour)
# ============================================================

@inline linidx(x::Int, y::Int) = (y - 1) * NX + x
@inline x_of(i::Int) = ((i - 1) % NX) + 1
@inline y_of(i::Int) = ((i - 1) ÷ NX) + 1

function build_neighbors_4()
    neigh = Vector{Vector{Int}}(undef, NCELLS)
    for i in 1:NCELLS
        x = x_of(i); y = y_of(i)
        nb = Int[]
        if x > 1;  push!(nb, linidx(x-1,y)); end
        if x < NX; push!(nb, linidx(x+1,y)); end
        if y > 1;  push!(nb, linidx(x,y-1)); end
        if y < NY; push!(nb, linidx(x,y+1)); end
        neigh[i] = nb
    end
    return neigh
end
const NEIGH_4 = build_neighbors_4()

# ============================================================
# 2) Connectivity utilities (LCC and largest component mask)
# ============================================================
mutable struct CCWorkspace
    seen::Vector{Int32}
    stamp::Int32
    queue::Vector{Int}
end

function make_workspaces()
    nt = Threads.nthreads()
    wss = Vector{CCWorkspace}(undef, nt)
    for t in 1:nt
        wss[t] = CCWorkspace(fill(Int32(0), NCELLS), Int32(0), Int[])
    end
    return wss
end
const WSS = make_workspaces()

function lcc_size(ws::CCWorkspace, mask::BitVector)
    count(mask) == 0 && return 0
    ws.stamp += 1
    stamp = ws.stamp
    seen = ws.seen
    empty!(ws.queue)

    best = 0
    @inbounds for i in 1:NCELLS
        if mask[i] && seen[i] != stamp
            seen[i] = stamp
            push!(ws.queue, i)
            qpos = 1
            compsize = 0
            while qpos <= length(ws.queue)
                v = ws.queue[qpos]; qpos += 1
                compsize += 1
                for nb in NEIGH_4[v]
                    if mask[nb] && seen[nb] != stamp
                        seen[nb] = stamp
                        push!(ws.queue, nb)
                    end
                end
            end
            best = max(best, compsize)
            empty!(ws.queue)
        end
    end
    return best
end

function largest_component_mask(ws::CCWorkspace, mask::BitVector)
    # returns largest component as mask; if none, returns all-false
    count(mask) == 0 && return BitVector(falses(NCELLS))
    ws.stamp += 1
    stamp = ws.stamp
    seen = ws.seen
    empty!(ws.queue)

    best = 0
    best_nodes = Int[]
    tmp_nodes = Int[]

    @inbounds for i in 1:NCELLS
        if mask[i] && seen[i] != stamp
            seen[i] = stamp
            push!(ws.queue, i)
            qpos = 1
            empty!(tmp_nodes)
            while qpos <= length(ws.queue)
                v = ws.queue[qpos]; qpos += 1
                push!(tmp_nodes, v)
                for nb in NEIGH_4[v]
                    if mask[nb] && seen[nb] != stamp
                        seen[nb] = stamp
                        push!(ws.queue, nb)
                    end
                end
            end
            if length(tmp_nodes) > best
                best = length(tmp_nodes)
                best_nodes = copy(tmp_nodes)
            end
            empty!(ws.queue)
        end
    end

    out = BitVector(falses(NCELLS))
    @inbounds for v in best_nodes
        out[v] = true
    end
    return out
end

# Override the old behavior:
# - if USE_CONNECTIVITY_FILTER = false -> return mask unchanged
# - else -> keep all components >= Emin (and drop the rest)
# ============================================================
# Connectivity filter (PATCH): keep ALL components ≥ Emin
# ============================================================

# Keeps every connected component (4-neigh) whose size is >= Emin
# Drops smaller components. If nothing survives -> all-false.
function keep_components_ge_Emin(ws::CCWorkspace, mask::BitVector, Emin::Int)
    count(mask) == 0 && return BitVector(falses(NCELLS))

    ws.stamp += 1
    stamp = ws.stamp
    seen = ws.seen
    empty!(ws.queue)

    out = BitVector(falses(NCELLS))
    comp_nodes = Int[]

    @inbounds for i in 1:NCELLS
        if mask[i] && seen[i] != stamp
            empty!(comp_nodes)

            # BFS collect this component
            seen[i] = stamp
            push!(ws.queue, i)
            qpos = 1

            while qpos <= length(ws.queue)
                v = ws.queue[qpos]; qpos += 1
                push!(comp_nodes, v)
                for nb in NEIGH_4[v]
                    if mask[nb] && seen[nb] != stamp
                        seen[nb] = stamp
                        push!(ws.queue, nb)
                    end
                end
            end

            # keep it if large enough
            if length(comp_nodes) >= Emin
                for v in comp_nodes
                    out[v] = true
                end
            end

            empty!(ws.queue)
        end
    end

    return out
end

function apply_connectivity_filter(ws::CCWorkspace, mask::BitVector, Emin::Int)
    if !USE_CONNECTIVITY_FILTER
        return mask
    end

    # cheap early-out: total area < Emin can't contain any component >= Emin
    if count(mask) < Emin
        return BitVector(falses(NCELLS))
    end

    kept = keep_components_ge_Emin(ws, mask, Emin)
    return (count(kept) == 0) ? BitVector(falses(NCELLS)) : kept
end

# ============================================================
# 3) Environmental fields
# ============================================================

function rescale_to_range!(v::Vector{Float64}, lo::Float64, hi::Float64)
    mn = minimum(v); mx = maximum(v)
    if mx == mn
        fill!(v, 0.5*(lo+hi))
        return v
    end
    @inbounds for i in eachindex(v)
        v[i] = lo + (hi-lo) * (v[i] - mn) / (mx - mn)
    end
    return v
end

function smooth_field_once!(E::Vector{Float64}, tmp::Vector{Float64}, α::Float64)
    @inbounds for i in 1:NCELLS
        s = E[i]
        nb = NEIGH_4[i]
        m = s
        for j in nb
            m += E[j]
        end
        m /= (1 + length(nb))
        tmp[i] = (1-α)*E[i] + α*m
    end
    copyto!(E, tmp)
    return E
end

function make_environment(rng::AbstractRNG, kind::Symbol)
    E = randn(rng, NCELLS)
    if kind == :autocorr
        tmp = similar(E)
        for _ in 1:AUTOCORR_ITERS
            smooth_field_once!(E, tmp, AUTOCORR_ALPHA)
        end
    elseif kind != :random
        error("Unknown env kind: $kind")
    end
    rescale_to_range!(E, E_MIN, E_MAX)
    return E
end

# ============================================================
# 4) Niches: suitability masks (1D Gaussian threshold)
# ============================================================

@inline function suitability_mask_1d(E::Vector{Float64}, μ::Float64, σ::Float64, thresh::Float64)
    lim = sqrt(-2.0 * log(thresh))  # |(E-μ)/σ| <= lim
    invσ = 1.0 / max(σ, 1e-6)
    m = BitVector(undef, NCELLS)
    @inbounds for i in 1:NCELLS
        z = (E[i] - μ) * invσ
        m[i] = abs(z) <= lim
    end
    return m
end

# Breadth regimes
struct BreadthRegime
    name::String
    meanσ::Float64
    logsd::Float64
end

const regimes = [
    BreadthRegime("Narrow + LowVar",  7.5, 0.20),
    BreadthRegime("Narrow + HighVar", 7.5, 0.55),
    BreadthRegime("Broad + LowVar",   16.0, 0.20),
    BreadthRegime("Broad + HighVar",  16.0, 0.55),
]

function draw_sigmas(rng::AbstractRNG, regime::BreadthRegime)
    σ = Vector{Float64}(undef, S)
    # lognormal around meanσ with logsd; clamp to reasonable bounds
    for i in 1:S
        val = exp(log(regime.meanσ) + regime.logsd * randn(rng))
        σ[i] = clamp(val, 1.5, 45.0)
    end
    return σ
end

# ============================================================
# 5) Network builders (4 families)
# Consumers are species nb+1:S; basal 1:nb have no prey.
# Each consumer is guaranteed ≥1 prey.
# ============================================================

function consumers_and_basal()
    nb = round(Int, BASAL_FRAC * S)
    basal_mask = BitVector(falses(S))
    basal_mask[1:nb] .= true
    consumers = collect((nb+1):S)
    return nb, basal_mask, consumers
end

function realized_connectance(prey::Vector{Vector{Int}}, basal_mask::BitVector)
    # realized L / S^2 (to match your earlier convention)
    L = 0
    for i in 1:S
        basal_mask[i] && continue
        L += length(prey[i])
    end
    return L / (S^2)
end

function ensure_min1_prey!(rng::AbstractRNG, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    # Guarantee each consumer has at least one prey.
    for i in 1:S
        basal_mask[i] && continue
        if isempty(prey[i])
            # prefer basal prey
            candidates = findall(basal_mask)
            if isempty(candidates)
                explained = [j for j in 1:S if j != i]
                push!(prey[i], explained[rand(rng, 1:length(explained))])
            else
                push!(prey[i], candidates[rand(rng, 1:length(candidates))])
            end
        end
    end
end

function build_metaweb_random(rng::AbstractRNG, C::Float64, basal_mask::BitVector)
    # Random consumer->prey edges; exact Ltarget achieved (after ensuring 1 prey each).
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal_mask)
    Ltarget = round(Int, C * S^2)
    # Step 1: 1 prey each
    for i in consumers
        cand = [j for j in 1:S if j != i]
        push!(prey[i], cand[rand(rng, 1:length(cand))])
    end
    # Step 2: fill remaining edges uniformly without duplicates
    L = length(consumers)
    while L < Ltarget
        i = consumers[rand(rng, 1:length(consumers))]
        cand = rand(rng, 1:S)
        if cand == i; continue; end
        if cand ∉ prey[i]
            push!(prey[i], cand)
            L += 1
        end
    end
    return prey
end

function build_metaweb_modular(rng::AbstractRNG, C::Float64, basal_mask::BitVector)
    # Modular SBM-like: higher probability within module.
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal_mask)
    # module assignment for all species
    MODULE = Vector{Int}(undef, S)
    for i in 1:S
        MODULE[i] = 1 + (i - 1) % N_MODULES
    end
    Ltarget = round(Int, C * S^2)

    # helper: sample prey with within-module bias
    function sample_prey(i::Int)
        inmod = Int[]
        outmod = Int[]
        mi = MODULE[i]
        for j in 1:S
            if j == i; continue; end
            if MODULE[j] == mi
                push!(inmod, j)
            else
                push!(outmod, j)
            end
        end
        # weight within vs outside
        if isempty(inmod)
            return outmod[rand(rng, 1:length(outmod))]
        elseif isempty(outmod)
            return inmod[rand(rng, 1:length(inmod))]
        else
            if rand(rng) < MODULAR_IN_BIAS / (MODULAR_IN_BIAS + 1.0)
                return inmod[rand(rng, 1:length(inmod))]
            else
                return outmod[rand(rng, 1:length(outmod))]
            end
        end
    end

    # Step 1: 1 prey each
    for i in consumers
        push!(prey[i], sample_prey(i))
    end
    # Step 2: fill remaining
    L = length(consumers)
    while L < Ltarget
        i = consumers[rand(rng, 1:length(consumers))]
        j = sample_prey(i)
        if j ∉ prey[i]
            push!(prey[i], j)
            L += 1
        end
    end
    return prey, MODULE
end

function build_metaweb_heavytail(rng::AbstractRNG, C::Float64, basal_mask::BitVector)
    # Heavy-tailed consumer out-degree, overall Ltarget.
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal_mask)
    nC = length(consumers)
    Ltarget = round(Int, C * S^2)

    # weights ~ Pareto-like
    w = zeros(Float64, nC)
    for k in 1:nC
        # Pareto with exponent gamma: use inverse-CDF with u^( -1/(gamma-1) )
        u = rand(rng)
        w[k] = u^(-1/(HEAVYTAIL_GAMMA-1))
    end
    w ./= sum(w)

    # allocate degrees: start with 1 each
    deg = ones(Int, nC)
    remaining = max(0, Ltarget - nC)
    # multinomial allocation of remaining edges
    for _ in 1:remaining
        r = rand(rng)
        acc = 0.0
        idx = 1
        for k in 1:nC
            acc += w[k]
            if r <= acc
                idx = k
                break
            end
        end
        deg[idx] += 1
    end

    # cap degrees and reassign overflow to keep total ~ Ltarget
    kmax = max(2, round(Int, HEAVYTAIL_KMAX_FRAC * (S-1)))
    overflow = 0
    for k in 1:nC
        if deg[k] > kmax
            overflow += deg[k] - kmax
            deg[k] = kmax
        end
    end
    # redistribute overflow
    for _ in 1:overflow
        k = rand(rng, 1:nC)
        if deg[k] < kmax
            deg[k] += 1
        end
    end

    # choose prey sets
    for (kk, i) in enumerate(consumers)
        candidates = [j for j in 1:S if j != i]
        shuffle!(rng, candidates)
        d = min(deg[kk], length(candidates))
        prey[i] = candidates[1:d]
    end
    ensure_min1_prey!(rng, prey, basal_mask)
    return prey
end

function build_metaweb_cascade(rng::AbstractRNG, C::Float64, basal_mask::BitVector)
    # Cascade hierarchy: consumers feed on lower-ranked species; optional interval bias via CASCADE_LAMBDA.
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal_mask)
    nC = length(consumers)
    Ltarget = round(Int, C * S^2)

    # ranks for all species; enforce every consumer has at least one lower-ranked prey candidate
    ranks = zeros(Float64, S)
    for attempt in 1:200
        for i in 1:S
            ranks[i] = rand(rng)
        end
        ok = true
        for i in consumers
            lower = findall(j -> ranks[j] < ranks[i] && j != i, 1:S)
            if isempty(lower)
                ok = false
                break
            end
        end
        ok && break
        attempt == 200 && error("Failed to sample cascade ranks with valid lower-prey candidates for all consumers")
    end

    # helper: sample prey among lower ranks with exponential bias toward nearby ranks
    function sample_lower_prey(i::Int)
        lower = Int[]
        w = Float64[]
        ri = ranks[i]
        for j in 1:S
            if j == i; continue; end
            if ranks[j] < ri
                push!(lower, j)
                # weight: exp(-λ * (ri-rj))
                push!(w, exp(-CASCADE_LAMBDA * (ri - ranks[j])))
            end
        end
        # normalize and sample
        sw = sum(w)
        r = rand(rng) * sw
        acc = 0.0
        for k in 1:length(lower)
            acc += w[k]
            if r <= acc
                return lower[k]
            end
        end
        return lower[end]
    end

    # Step 1: 1 prey each
    for i in consumers
        push!(prey[i], sample_lower_prey(i))
    end
    # Step 2: fill remaining edges
    L = nC
    while L < Ltarget
        i = consumers[rand(rng, 1:nC)]
        j = sample_lower_prey(i)
        if j ∉ prey[i]
            push!(prey[i], j)
            L += 1
        end
    end

    return prey, ranks
end

# ============================================================
# 6) Mechanistic niche correlation builder
# Target correlation across consumers between μ_cons and mean μ_prey(cons)
# ============================================================

function pearson_r(a::Vector{Float64}, b::Vector{Float64})
    ma = mean(a); mb = mean(b)
    sa = std(a); sb = std(b)
    (sa == 0 || sb == 0) && return 0.0
    return mean((a .- ma) .* (b .- mb)) / (sa * sb)
end

function prey_means(mu::Vector{Float64}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    pm = fill(NaN, S)
    for i in 1:S
        if basal_mask[i] || isempty(prey[i])
            continue
        end
        pm[i] = mean(mu[prey[i]])
    end
    return pm
end

function mechanistic_corr(mu::Vector{Float64}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    pm = prey_means(mu, prey, basal_mask)
    idx = findall(i -> !basal_mask[i] && !isnan(pm[i]), 1:S)
    length(idx) < 6 && return 0.0
    return pearson_r(mu[idx], pm[idx])
end

function assign_mus_with_target_corr!(
    rng::AbstractRNG,
    mu::Vector{Float64},
    prey::Vector{Vector{Int}},
    basal_mask::BitVector,
    target_r::Float64;
    relax_iters::Int = RELAX_ITERS
)
    # Basal μ already set; this updates consumer μ using relaxation with α.
    consumers = findall(!, basal_mask)

    function relax(alpha::Float64)
        m = copy(mu)
        for _ in 1:relax_iters
            pm = prey_means(m, prey, basal_mask)
            @inbounds for i in consumers
                if isnan(pm[i]); continue; end
                m[i] = clamp((1-alpha)*m[i] + alpha*pm[i] + MU_NOISE_SD*randn(rng), E_MIN, E_MAX)
            end
        end
        return mechanistic_corr(m, prey, basal_mask), m
    end

    lo, hi = -0.98, 0.98
    best_err = Inf
    best_mu = copy(mu)
    best_r = 0.0

    for _ in 1:26
        mid = 0.5*(lo + hi)
        rmid, mmid = relax(mid)
        err = abs(rmid - target_r)
        if err < best_err
            best_err = err
            best_mu = mmid
            best_r = rmid
        end
        if rmid < target_r
            lo = mid
        else
            hi = mid
        end
    end

    copyto!(mu, best_mu)
    return best_r
end

# ============================================================
# 7) AB fixed point with one-prey rule
# ============================================================

function fixed_point_AB(A_pres::Vector{BitVector}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    pres = [copy(A_pres[i]) for i in 1:S]
    newp = [BitVector(falses(NCELLS)) for _ in 1:S]
    changed = true
    iter = 0
    while changed
        iter += 1
        changed = false
        for i in 1:S
            if basal_mask[i] || isempty(prey[i])
                newp[i] = pres[i]
            else
                u = BitVector(falses(NCELLS))
                @inbounds for j in prey[i]
                    u .|= pres[j]
                end
                newp[i] = A_pres[i] .& u
            end
        end
        for i in 1:S
            if newp[i] != pres[i]
                pres[i] = newp[i]
                changed = true
            end
        end
        iter > 80 && break
    end
    return pres
end

# ============================================================
# 8) Metrics
# ============================================================

function gamma_richness_cons(pres::Vector{BitVector}, basal_mask::BitVector)
    c = 0
    for i in 1:S
        basal_mask[i] && continue
        c += (count(pres[i]) > 0) ? 1 : 0
    end
    return c
end

function mean_jaccard_mismatch(A::Vector{BitVector}, AB::Vector{BitVector}, basal_mask::BitVector)
    vals = Float64[]
    for i in 1:S
        basal_mask[i] && continue
        Ai = A[i]
        count(Ai) == 0 && continue   # exclude climatically absent consumers

        ABi = AB[i]
        inter = count(Ai .& ABi)
        uni   = count(Ai .| ABi)
        J = uni == 0 ? 1.0 : (inter / uni)
        push!(vals, 1 - J)
    end
    return isempty(vals) ? NaN : mean(vals)
end

function frac_affected(A::Vector{BitVector}, AB::Vector{BitVector}, basal_mask::BitVector)
    num = 0
    den = 0
    for i in 1:S
        basal_mask[i] && continue
        Ai = A[i]
        count(Ai) == 0 && continue
        den += 1
        num += (Ai != AB[i]) ? 1 : 0
    end
    return den == 0 ? NaN : num / den
end

function realized_overlap(A_raw::Vector{BitVector}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    # O_i = avg_j |A_i ∩ A_j| / |A_i|
    vals = Float64[]
    for i in 1:S
        basal_mask[i] && continue
        Ai = A_raw[i]
        a = count(Ai)
        (a == 0 || isempty(prey[i])) && continue
        s = 0.0
        for j in prey[i]
            s += count(Ai .& A_raw[j]) / a
        end
        push!(vals, s / length(prey[i]))
    end
    return isempty(vals) ? NaN : mean(vals)
end

# ============================================================
# 9) One replicate at (envkind, networkfamily, regime, C, target_r)
# ============================================================

@inline function count_links(prey::Vector{Vector{Int}}, basal_mask::BitVector)
    L = 0
    for i in 1:S
        basal_mask[i] && continue
        L += length(prey[i])
    end
    return L
end

function simulate_one!(
    rng::AbstractRNG,
    ws::CCWorkspace,
    envkind::Symbol,
    netfamily::Symbol,
    regime::BreadthRegime,
    C::Float64,
    target_r::Float64
)
    nb, basal_mask, consumers = consumers_and_basal()

    # Environment
    E = make_environment(rng, envkind)

    # Network
    MODULE = nothing
    ranks = nothing
    prey = nothing

    if netfamily == :random
        prey = build_metaweb_random(rng, C, basal_mask)
    elseif netfamily == :modular
        prey, MODULE = build_metaweb_modular(rng, C, basal_mask)
    elseif netfamily == :heavytail
        prey = build_metaweb_heavytail(rng, C, basal_mask)
    elseif netfamily == :cascade
        prey, ranks = build_metaweb_cascade(rng, C, basal_mask)
    else
        error("Unknown netfamily: $netfamily")
    end

    # Niche parameters: μ and σ
    σ = draw_sigmas(rng, regime)

    μ = Vector{Float64}(undef, S)
    # basal μ uniform
    for i in 1:nb
        μ[i] = rand(rng) * (E_MAX - E_MIN) + E_MIN
    end
    # initial consumer μ uniform
    for i in (nb+1):S
        μ[i] = rand(rng) * (E_MAX - E_MIN) + E_MIN
    end

    achieved_r = assign_mus_with_target_corr!(rng, μ, prey, basal_mask, target_r)

    # A (raw climatic)
    A_raw = Vector{BitVector}(undef, S)
    for i in 1:S
        A_raw[i] = suitability_mask_1d(E, μ[i], σ[i], SUIT_THRESH)
    end

    # AB fixed point (raw)
    AB_raw = fixed_point_AB(A_raw, prey, basal_mask)

    # Apply connectivity filter post-hoc (extinction rule)
    A = Vector{BitVector}(undef, S)
    AB = Vector{BitVector}(undef, S)
    for i in 1:S
        A[i]  = apply_connectivity_filter(ws, A_raw[i], Emin_patch)
        AB[i] = apply_connectivity_filter(ws, AB_raw[i], Emin_patch)
    end

    # Metrics (consumers-only)
    SA  = gamma_richness_cons(A, basal_mask)
    SAB = gamma_richness_cons(AB, basal_mask)
    dSrel = (SA == 0) ? NaN : (1.0 - (SAB / SA))

    mjm = mean_jaccard_mismatch(A, AB, basal_mask)
    fa  = frac_affected(A, AB, basal_mask)
    ov  = realized_overlap(A_raw, prey, basal_mask)

    # network diagnostics
    L = count_links(prey, basal_mask)
    Creal = L / (S^2)

    return (
        dSrel=dSrel,
        mean_jaccard_mismatch=mjm,
        frac_affected=fa,
        realized_overlap=ov,
        achieved_r=achieved_r,
        Creal=Creal
    )
end

# ============================================================
# 10) Sweep runner
# ============================================================

const ENVKINDS = [:random, :autocorr]
const NETFAMS  = [:random, :modular, :heavytail, :cascade]
const NETNAMES = Dict(
    :random=>"Random",
    :modular=>"Modular",
    :heavytail=>"Heavy-tail",
    :cascade=>"Cascade"
)

function regime_name(reg::BreadthRegime)
    return reg.name
end

function sweep_all()
    Cvals = collect(range(0.01, 0.1, length=N_CONNECT))
    Rvals = collect(range(0.0, 0.9, length=N_CORR))

    # Store: Dict keyed by (env, netfam, regime_index, metric) => Matrix(N_CORR, N_CONNECT)
    metrics = [:dSrel, :mean_jaccard_mismatch, :frac_affected, :realized_overlap, :achieved_r, :Creal]
    store = Dict{Tuple{Symbol,Symbol,Int,Symbol}, Matrix{Float64}}()

    for env in ENVKINDS, net in NETFAMS, (ri, reg) in enumerate(regimes), m in metrics
        store[(env, net, ri, m)] = fill(NaN, N_CORR, N_CONNECT)
    end

    println("Threads: ", Threads.nthreads())
    println("Grid: $(NX)x$(NY) cells=$(NCELLS), S=$(S), basal_frac=$(BASAL_FRAC), Emin=$(Emin_patch), connectivity_filter=$(USE_CONNECTIVITY_FILTER)")
    println("OUTDIR: ", OUTDIR)
    println("Sweep: env=$(length(ENVKINDS)) × net=$(length(NETFAMS)) × regimes=$(length(regimes)) × (C=$(N_CONNECT), r=$(N_CORR)) × reps=$(NREP)\n")

    total_cells = length(ENVKINDS)*length(NETFAMS)*length(regimes)*N_CONNECT*N_CORR
    cell_counter = 0

    for env in ENVKINDS
        for net in NETFAMS
            for (ri, reg) in enumerate(regimes)

                # parallelize over grid cells
                Threads.@threads for idx in 1:(N_CONNECT * N_CORR)
                    tid = Threads.threadid()
                    ws = WSS[tid]
                    rng = MersenneTwister(BASE_SEED + 100_000*tid + 17)

                    cind = (idx - 1) % N_CONNECT + 1
                    rind = (idx - 1) ÷ N_CONNECT + 1
                    C = Cvals[cind]
                    r = Rvals[rind]

                    dS = Float64[]
                    jm = Float64[]
                    fa = Float64[]
                    ov = Float64[]
                    rg = Float64[]
                    cr = Float64[]

                    for rep in 1:NREP
                        seed = BASE_SEED +
                               10_000_000 * (findfirst(==(env), ENVKINDS)) +
                               1_000_000  * (findfirst(==(net), NETFAMS)) +
                               100_000    * ri +
                               10_000     * rind +
                               1_000      * cind +
                               rep
                        Random.seed!(rng, seed)
                        out = simulate_one!(rng, ws, env, net, reg, C, r)
                        push!(dS, out.dSrel)
                        push!(jm, out.mean_jaccard_mismatch)
                        push!(fa, out.frac_affected)
                        push!(ov, out.realized_overlap)
                        push!(rg, out.achieved_r)
                        push!(cr, out.Creal)
                    end

                    store[(env, net, ri, :dSrel)][rind, cind] = mean(dS)
                    store[(env, net, ri, :mean_jaccard_mismatch)][rind, cind] = mean(jm)
                    store[(env, net, ri, :frac_affected)][rind, cind] = mean(fa)
                    store[(env, net, ri, :realized_overlap)][rind, cind] = mean(ov)
                    store[(env, net, ri, :achieved_r)][rind, cind] = mean(rg)
                    store[(env, net, ri, :Creal)][rind, cind] = mean(cr)
                end

                cell_counter += N_CONNECT*N_CORR
                println("Done: env=$(env), net=$(net), regime=$(reg.name)  ($(cell_counter)/$(total_cells) cells)")
            end
        end
    end

    return store, Cvals, Rvals
end

# ============================================================
# 11) Plotting: facet grid (net families × regimes) of heatmaps
# ============================================================

function global_minmax(mats::Vector{Matrix{Float64}})
    vals = Float64[]
    for M in mats
        for x in M
            isfinite(x) && push!(vals, x)
        end
    end
    isempty(vals) && return (0.0, 1.0)
    return (minimum(vals), maximum(vals))
end

function facet_heatmaps(store, Cvals, Rvals, env::Symbol, metric::Symbol;
                        title::String="", outfile::Union{Nothing,String}=nothing,
                        fixed_colorbar = false)

    netfams = (:random, :modular, :heavytail, :cascade)
    reg_titles = ("Narrow + LowVar", "Narrow + HighVar",
                  "Broad + LowVar", "Broad + HighVar")
    net_titles = Dict(
        :random    => "Random",
        :modular   => "Modular",
        :heavytail => "Heavy-tail",
        :cascade   => "Cascade"
    )

    nC = length(Cvals)
    nR = length(Rvals)

    # --- Keep your storage orientation safety (but we will NOT rely on heatmap!(ax, M))
    # We enforce the semantic convention: M_RC[r, c] => (Rvals[r], Cvals[c])
    function orient_to_RC(M::AbstractMatrix)
        nr, nc = size(M)
        if nr == nR && nc == nC
            return M                      # already (R,C)
        elseif nr == nC && nc == nR
            return permutedims(M)         # stored as (C,R) -> convert to (R,C)
        else
            error("Matrix size $(size(M)) incompatible with (R=$nR, C=$nC)")
        end
    end

    f = Figure(size=(1600, 900))

    # ---- explicit grid (NO resizing calls)
    gl = f[1, 1] = GridLayout()
    rowgap!(gl, 14)
    colgap!(gl, 14)

    # ---- title
    if !isempty(title)
        Label(gl[1, 1:6], title, fontsize=11)
    end

    # ---- column headers
    for ci in 1:4
        Label(gl[2, ci+1], reg_titles[ci], fontsize=10, halign=:center)
    end

    xtlabs = [@sprintf("%.2f", x) for x in Cvals]
    ytlabs = [@sprintf("%.2f", y) for y in Rvals]

    hobj = nothing

    # IMPORTANT FIX:
    # Makie heatmap uses x-length == size(Z,1) and y-length == size(Z,2).
    # We want x = C (length nC) and y = R (length nR),
    # so Z MUST be sized (nC, nR) at plotting time.
    # If our semantic matrix is M_RC sized (nR, nC), we plot permutedims(M_RC) (=> nC x nR).
    xcoords = 1:nC
    ycoords = 1:nR

    # ---- panels
    for (net_i, net) in enumerate(netfams)

        Label(gl[net_i+2, 1], net_titles[net],
              rotation=pi/2, fontsize=10, valign=:center)

        for ri in 1:4
            ax = Axis(
                gl[net_i+2, ri+1],
                xticks = (xcoords, xtlabs),
                yticks = (ycoords, ytlabs),
                xlabel = "Connectance",
                ylabel = "Niche corr",
                xticklabelrotation = pi/4,
                xticklabelsize = 8,
                yticklabelsize = 8
            )

            Mraw = store[(env, net, ri, metric)]
            M_RC = orient_to_RC(Mraw)                 # (nR, nC) with semantic meaning
            Z = permutedims(M_RC)                     # (nC, nR) for Makie heatmap x/y mapping

            if !fixed_colorbar
                h = heatmap!(ax, xcoords, ycoords, Z)     # AXIS-SAFE
            else
                h = heatmap!(ax, xcoords, ycoords, Z, colorrange = (0.0, 1.0))
            end

            hobj === nothing && (hobj = h)

            # reduce clutter (unchanged)
            net_i < length(netfams) && (ax.xticklabelsvisible = false)
            ri > 1 && (ax.yticklabelsvisible = false)
        end
    end

    # ---- colorbar
    if hobj !== nothing
        Colorbar(gl[3:6, 6], hobj, label=string(metric))
    end

    if outfile !== nothing
        save(joinpath(OUTDIR, outfile), f)
    end

    return f
end

# ============================================================
# 12) MAIN
# ============================================================
store, Cvals, Rvals = sweep_all()
using Serialization

# cache_path = joinpath(OUTDIR, "sweep_cache_smallerConn_001_smallerRrange_60x60cells.jls")
# serialize(cache_path, (store=store, Cvals=Cvals, Rvals=Rvals))
# println("Saved sweep cache to: ", cache_path)

cache_path = joinpath(OUTDIR, "sweep_cache_smallerConn_001_smallerRrange_60x60cells.jls")
data = deserialize(cache_path)

store = data.store
Cvals = data.Cvals
Rvals = data.Rvals

println("Loaded sweep cache from: ", cache_path)

# Save numeric matrices as TSV (simple, no external packages)
function save_matrix_tsv(path::String, M::Matrix{Float64})
    open(path, "w") do io
        nr, nc = size(M)
        for r in 1:nr
            println(io, join([@sprintf("%.6f", M[r,c]) for c in 1:nc], '\t'))
        end
    end
end

println("\nExploding outputs...")
for env in ENVKINDS, net in NETFAMS, (ri, reg) in enumerate(regimes)
    for metric in [:dSrel, :mean_jaccard_mismatch, :frac_affected, :realized_overlap, :achieved_r, :Creal]
        M = store[(env, net, ri, metric)]
        fname = "mat_$(env)_$(net)_reg$(ri)_$(metric).tsv"
        save_matrix_tsv(joinpath(OUTDIR, fname), M)
    end
end

# Plot heatmaps for the 3 main divergence metrics + realized overlap (mediator diagnostic)
for env in ENVKINDS
    envname = env == :random ? "Random environment" : "Autocorrelated environment"

    f1 = facet_heatmaps(
        store, Cvals, Rvals, env, :dSrel;
        title = "Relative richness loss (consumers-only): 1 - S_AB / S_A — $(envname)",
        outfile = "heatmaps_$(env)_metric_dSrel_smallerConn_smallerR_60x60.png",
        fixed_colorbar = true
    )
    display(f1)

    f2 = facet_heatmaps(
        store, Cvals, Rvals, env, :mean_jaccard_mismatch;
        title = "Mean per-species Jaccard mismatch (consumers-only): mean(1 - J(A_i,AB_i)) — $(envname)",
        outfile = "heatmaps_$(env)_metric_mean_jaccard_mismatch_smallerConn_smallerR_60x60.png",
        # fixed_colorbar = true
    )
    display(f2)

    f3 = facet_heatmaps(
        store, Cvals, Rvals, env, :frac_affected;
        title = "Fraction affected (consumers-only): frac(A_i != AB_i) — $(envname)",
        outfile = "heatmaps_$(env)_metric_frac_affected_smallerConn_smallerR_60x60.png",
        fixed_colorbar = true
    )
    display(f3)

    f4 = facet_heatmaps(
        store, Cvals, Rvals, env, :realized_overlap;
        title = "Realized prey-support overlap: mean_i avg_j |A_i∩A_j|/|A_i| — $(envname)",
        outfile = "heatmaps_$(env)_diagnostic_realized_overlap_smallerConn_smallerR_60x60.png",
        fixed_colorbar = true
    )
    display(f4)

    f5 = facet_heatmaps(
        store, Cvals, Rvals, env, :achieved_r;
        title = "Achieved mechanistic niche correlation — $(envname)",
        outfile = "heatmaps_$(env)_diagnostic_achieved_r_smallerConn_smallerR_60x60.png",
        fixed_colorbar = false
    )
    display(f5)

    f6 = facet_heatmaps(
        store, Cvals, Rvals, env, :Creal;
        title = "Realized connectance: L/S^2 — $(envname)",
        outfile = "heatmaps_$(env)_diagnostic_Creaal_smallerConn_smallerR_60x60.png",
        fixed_colorbar = false
    )
    display(f6)
end

println("\nDone. Output directory:")
println(OUTDIR)