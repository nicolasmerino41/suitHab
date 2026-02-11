#!/usr/bin/env julia
# Julia 1.11
# Run:  julia --threads auto biotic_divergence_jaccard_tail_pipeline.jl
#
# PURPOSE (self-contained, no external data)
# ------------------------------------------------------------
# A vs AB baseline divergence (no habitat loss), focusing on PER-SPECIES Jaccard mismatch:
#   mismatch_i = 1 - J(A_i, AB_i)  for consumer species with |A_i|>0
#
# Outputs:
# 1) Heatmaps (faceted: net families × niche-breadth regimes) for:
#    - mismatch_q90: 90th percentile of mismatch across consumers
#    - mismatch_frac_gt: fraction of consumers with mismatch > TAIL_THRESH
#    (optionally: achieved mechanistic r and realized connectance diagnostics)
# 2) A few distribution figures for selected parameter cells (hist + ECDF)
#
# Notes:
# - A: climatic suitability only
# - AB: climatic suitability + one-prey rule fixed point
# - Optional spatial connectivity filter: keep only largest component if LCC ≥ Emin (applied after A/AB built)
# - Networks: Random, Modular, Heavy-tail out-degree, Cascade hierarchy
# - Environmental space: Random vs Autocorrelated field
# - Mechanistic niche correlation: correlation across consumers between μ_cons and mean μ_prey(cons)

using Random
using Statistics
using LinearAlgebra
using CairoMakie
using Printf
using Dates
using Serialization

# ============================================================
# 0) Global parameters (tune here)
# ============================================================

# Spatial grid
const NX = 60
const NY = 60
const NCELLS = NX * NY

# Species pool
const S = 250
const BASAL_FRAC = 0.30

# Spatial viability filter (movement/connectivity proxy)
const USE_CONNECTIVITY_FILTER = true
const Emin_patch = 60

# Environmental field domain (e.g., temperature)
const E_MIN = 0.0
const E_MAX = 100.0

# Niche suitability: Gaussian with threshold
const SUIT_THRESH = 0.25

# Environmental autocorrelation
const AUTOCORR_ITERS = 18
const AUTOCORR_ALPHA = 0.55

# Sweep axes
const CONNECTANCE_RANGE = (0.005, 0.1)
const CORR_RANGE       = (0.0, 0.9)
const N_CONNECT = 15
const N_CORR    = 15

# Replicates per heatmap cell
const NREP = 8

# Tail threshold for "fraction strongly affected"
const TAIL_THRESH = 0.8

# Network-family knobs
const N_MODULES = 6
const MODULAR_IN_BIAS = 6.0
const HEAVYTAIL_GAMMA = 2.2
const HEAVYTAIL_KMAX_FRAC = 0.35
const CASCADE_LAMBDA = 2.5

# Mechanistic niche-correlation builder knobs
const RELAX_ITERS = 30
const MU_NOISE_SD = 1.8
const TARGET_R_TOL = 0.03

# Thread-safe seeds
const BASE_SEED = 20260202

# Output directory
# ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
OUTDIR = joinpath(pwd(), "output_jaccard_tail_60x60_")
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

function apply_connectivity_filter(ws::CCWorkspace, mask::BitVector, Emin::Int)
    if !USE_CONNECTIVITY_FILTER
        return mask
    end
    if count(mask) < Emin
        return BitVector(falses(NCELLS))
    end
    lcc = lcc_size(ws, mask)
    if lcc < Emin
        return BitVector(falses(NCELLS))
    end
    return largest_component_mask(ws, mask)
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
        m = E[i]
        nb = NEIGH_4[i]
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
    for i in 1:S
        val = exp(log(regime.meanσ) + regime.logsd * randn(rng))
        σ[i] = clamp(val, 1.5, 45.0)
    end
    return σ
end

# ============================================================
# 5) Network builders (4 families)
# ============================================================

function consumers_and_basal()
    nb = round(Int, BASAL_FRAC * S)
    basal_mask = BitVector(falses(S))
    basal_mask[1:nb] .= true
    consumers = collect((nb+1):S)
    return nb, basal_mask, consumers
end

function ensure_min1_prey!(rng::AbstractRNG, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    for i in 1:S
        basal_mask[i] && continue
        if isempty(prey[i])
            candidates = findall(basal_mask)
            if isempty(candidates)
                cand = [j for j in 1:S if j != i]
                push!(prey[i], cand[rand(rng, 1:length(cand))])
            else
                push!(prey[i], candidates[rand(rng, 1:length(candidates))])
            end
        end
    end
end

function build_metaweb_random(rng::AbstractRNG, C::Float64, basal_mask::BitVector)
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal_mask)
    Ltarget = round(Int, C * S^2)

    for i in consumers
        cand = [j for j in 1:S if j != i]
        push!(prey[i], cand[rand(rng, 1:length(cand))])
    end

    L = length(consumers)
    while L < Ltarget
        i = consumers[rand(rng, 1:length(consumers))]
        j = rand(rng, 1:S)
        (j == i || j ∈ prey[i]) && continue
        push!(prey[i], j)
        L += 1
    end
    return prey
end

function build_metaweb_modular(rng::AbstractRNG, C::Float64, basal_mask::BitVector)
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal_mask)

    MODULE = Vector{Int}(undef, S)
    for i in 1:S
        MODULE[i] = 1 + (i - 1) % N_MODULES
    end
    Ltarget = round(Int, C * S^2)

    function sample_prey(i::Int)
        inmod = Int[]
        outmod = Int[]
        mi = MODULE[i]
        for j in 1:S
            (j == i) && continue
            if MODULE[j] == mi
                push!(inmod, j)
            else
                push!(outmod, j)
            end
        end
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

    for i in consumers
        push!(prey[i], sample_prey(i))
    end

    L = length(consumers)
    while L < Ltarget
        i = consumers[rand(rng, 1:length(consumers))]
        j = sample_prey(i)
        (j ∈ prey[i]) && continue
        push!(prey[i], j)
        L += 1
    end
    return prey, MODULE
end

function build_metaweb_heavytail(rng::AbstractRNG, C::Float64, basal_mask::BitVector)
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal_mask)
    nC = length(consumers)
    Ltarget = round(Int, C * S^2)

    w = zeros(Float64, nC)
    for k in 1:nC
        u = rand(rng)
        w[k] = u^(-1/(HEAVYTAIL_GAMMA-1))
    end
    w ./= sum(w)

    deg = ones(Int, nC)
    remaining = max(0, Ltarget - nC)
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

    kmax = max(2, round(Int, HEAVYTAIL_KMAX_FRAC * (S-1)))
    overflow = 0
    for k in 1:nC
        if deg[k] > kmax
            overflow += deg[k] - kmax
            deg[k] = kmax
        end
    end
    for _ in 1:overflow
        k = rand(rng, 1:nC)
        if deg[k] < kmax
            deg[k] += 1
        end
    end

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
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal_mask)
    nC = length(consumers)
    Ltarget = round(Int, C * S^2)

    ranks = zeros(Float64, S)
    for attempt in 1:200
        for i in 1:S
            ranks[i] = rand(rng)
        end
        ok = true
        for i in consumers
            lower = findall(j -> ranks[j] < ranks[i] && j != i, 1:S)
            isempty(lower) && (ok = false; break)
        end
        ok && break
        attempt == 200 && error("Failed to sample cascade ranks")
    end

    function sample_lower_prey(i::Int)
        lower = Int[]
        w = Float64[]
        ri = ranks[i]
        for j in 1:S
            (j == i) && continue
            if ranks[j] < ri
                push!(lower, j)
                push!(w, exp(-CASCADE_LAMBDA * (ri - ranks[j])))
            end
        end
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

    for i in consumers
        push!(prey[i], sample_lower_prey(i))
    end

    L = nC
    while L < Ltarget
        i = consumers[rand(rng, 1:nC)]
        j = sample_lower_prey(i)
        (j ∈ prey[i]) && continue
        push!(prey[i], j)
        L += 1
    end

    return prey, ranks
end

# ============================================================
# 6) Mechanistic niche correlation builder
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
    consumers = findall(!, basal_mask)

    function relax(alpha::Float64)
        m = copy(mu)
        for _ in 1:relax_iters
            pm = prey_means(m, prey, basal_mask)
            @inbounds for i in consumers
                isnan(pm[i]) && continue
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
# 8) Per-species mismatch vector + summaries
# ============================================================

function jaccard_mismatch_vec(A::Vector{BitVector}, AB::Vector{BitVector}, basal_mask::BitVector)
    vals = Float64[]
    for i in 1:S
        basal_mask[i] && continue
        Ai = A[i]
        count(Ai) == 0 && continue   # exclude climatically absent consumers
        ABi = AB[i]
        inter = count(Ai .& ABi)
        uni   = count(Ai .| ABi)
        J = (uni == 0) ? 1.0 : (inter / uni)
        push!(vals, 1.0 - J)
    end
    return vals
end

function q_from_sorted(v::Vector{Float64}, p::Float64)
    n = length(v)
    n == 0 && return NaN
    idx = clamp(round(Int, 1 + (n-1)*p), 1, n)
    return v[idx]
end

function mismatch_q90(vals::Vector{Float64})
    isempty(vals) && return NaN
    v = sort(vals)
    return q_from_sorted(v, 0.90)
end

function mismatch_frac_gt(vals::Vector{Float64}, thr::Float64)
    isempty(vals) && return NaN
    return mean(vals .> thr)
end

# Diagnostics
@inline function count_links(prey::Vector{Vector{Int}}, basal_mask::BitVector)
    L = 0
    for i in 1:S
        basal_mask[i] && continue
        L += length(prey[i])
    end
    return L
end

# ============================================================
# 9) One replicate at (envkind, networkfamily, regime, C, target_r)
# ============================================================

function simulate_one!(
    rng::AbstractRNG,
    ws::CCWorkspace,
    envkind::Symbol,
    netfamily::Symbol,
    regime::BreadthRegime,
    C::Float64,
    target_r::Float64;
    return_mvec::Bool=false
)
    nb, basal_mask, consumers = consumers_and_basal()

    E = make_environment(rng, envkind)

    prey = nothing
    if netfamily == :random
        prey = build_metaweb_random(rng, C, basal_mask)
    elseif netfamily == :modular
        prey, _ = build_metaweb_modular(rng, C, basal_mask)
    elseif netfamily == :heavytail
        prey = build_metaweb_heavytail(rng, C, basal_mask)
    elseif netfamily == :cascade
        prey, _ = build_metaweb_cascade(rng, C, basal_mask)
    else
        error("Unknown netfamily: $netfamily")
    end

    σ = draw_sigmas(rng, regime)

    μ = Vector{Float64}(undef, S)
    for i in 1:S
        μ[i] = rand(rng) * (E_MAX - E_MIN) + E_MIN
    end

    achieved_r = assign_mus_with_target_corr!(rng, μ, prey, basal_mask, target_r)

    A_raw = Vector{BitVector}(undef, S)
    for i in 1:S
        A_raw[i] = suitability_mask_1d(E, μ[i], σ[i], SUIT_THRESH)
    end

    AB_raw = fixed_point_AB(A_raw, prey, basal_mask)

    # Apply connectivity filter post-hoc
    A = Vector{BitVector}(undef, S)
    AB = Vector{BitVector}(undef, S)
    for i in 1:S
        A[i]  = apply_connectivity_filter(ws, A_raw[i], Emin_patch)
        AB[i] = apply_connectivity_filter(ws, AB_raw[i], Emin_patch)
    end

    mvec = jaccard_mismatch_vec(A, AB, basal_mask)
    q90 = mismatch_q90(mvec)
    frac = mismatch_frac_gt(mvec, TAIL_THRESH)

    L = count_links(prey, basal_mask)
    Creal = L / (S^2)

    return (
        mismatch_q90 = q90,
        mismatch_frac_gt = frac,
        achieved_r = achieved_r,
        Creal = Creal,
        mvec = return_mvec ? mvec : Float64[]
    )
end

# ============================================================
# 10) Sweep runner (heatmaps)
# ============================================================

const ENVKINDS = [:random, :autocorr]
const NETFAMS  = [:random, :modular, :heavytail, :cascade]

function sweep_all()
    Cvals = collect(range(CONNECTANCE_RANGE[1], CONNECTANCE_RANGE[2], length=N_CONNECT))
    Rvals = collect(range(0.0, 0.9, length=N_CORR))

    metrics = [:mismatch_q90, :mismatch_frac_gt, :achieved_r, :Creal]
    store = Dict{Tuple{Symbol,Symbol,Int,Symbol}, Matrix{Float64}}()

    for env in ENVKINDS, net in NETFAMS, (ri, _) in enumerate(regimes), m in metrics
        store[(env, net, ri, m)] = fill(NaN, N_CORR, N_CONNECT) # semantic: (R,C)
    end

    println("Threads: ", Threads.nthreads())
    println("Grid: $(NX)x$(NY) cells=$(NCELLS), S=$(S), basal_frac=$(BASAL_FRAC), Emin=$(Emin_patch), connectivity_filter=$(USE_CONNECTIVITY_FILTER)")
    println("OUTDIR: ", OUTDIR)
    println("Sweep: env=$(length(ENVKINDS)) × net=$(length(NETFAMS)) × regimes=$(length(regimes)) × (C=$(N_CONNECT), r=$(N_CORR)) × reps=$(NREP)\n")

    total_cells = length(ENVKINDS)*length(NETFAMS)*length(regimes)*N_CONNECT*N_CORR
    done_cells = 0

    for env in ENVKINDS
        for net in NETFAMS
            for (ri, reg) in enumerate(regimes)

                Threads.@threads for idx in 1:(N_CONNECT * N_CORR)
                    tid = Threads.threadid()
                    ws = WSS[tid]
                    rng = MersenneTwister(BASE_SEED + 100_000*tid + 17)

                    cind = (idx - 1) % N_CONNECT + 1
                    rind = (idx - 1) ÷ N_CONNECT + 1
                    C = Cvals[cind]
                    r = Rvals[rind]

                    q90s = Float64[]
                    fracs = Float64[]
                    rgs = Float64[]
                    crs = Float64[]

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
                        push!(q90s, out.mismatch_q90)
                        push!(fracs, out.mismatch_frac_gt)
                        push!(rgs, out.achieved_r)
                        push!(crs, out.Creal)
                    end

                    store[(env, net, ri, :mismatch_q90)][rind, cind] = mean(q90s)
                    store[(env, net, ri, :mismatch_frac_gt)][rind, cind] = mean(fracs)
                    store[(env, net, ri, :achieved_r)][rind, cind] = mean(rgs)
                    store[(env, net, ri, :Creal)][rind, cind] = mean(crs)
                end

                done_cells += N_CONNECT*N_CORR
                println("Done: env=$(env), net=$(net), regime=$(reg.name)  ($(done_cells)/$(total_cells) cells)")
            end
        end
    end

    return store, Cvals, Rvals
end

# ============================================================
# 11) Plotting: facet grid (net families × regimes) of heatmaps
#     AXIS-SAFE orientation (semantic M[R,C] plotted with x=C, y=R)
# ============================================================

function facet_heatmaps(store, Cvals, Rvals, env::Symbol, metric::Symbol;
                        title::String="", outfile::Union{Nothing,String}=nothing)

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

    function orient_to_RC(M::AbstractMatrix)
        nr, nc = size(M)
        if nr == nR && nc == nC
            return M
        elseif nr == nC && nc == nR
            return permutedims(M)
        else
            error("Matrix size $(size(M)) incompatible with (R=$nR, C=$nC)")
        end
    end

    f = Figure(size=(1600, 900))
    gl = f[1, 1] = GridLayout()
    rowgap!(gl, 14)
    colgap!(gl, 14)

    if !isempty(title)
        Label(gl[1, 1:6], title, fontsize=11)
    end

    for ci in 1:4
        Label(gl[2, ci+1], reg_titles[ci], fontsize=10, halign=:center)
    end

    xtlabs = [@sprintf("%.3f", x) for x in Cvals]
    ytlabs = [@sprintf("%.2f", y) for y in Rvals]

    hobj = nothing
    xcoords = 1:nC
    ycoords = 1:nR

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
            M_RC = orient_to_RC(Mraw)     # semantic (R,C)
            Z = permutedims(M_RC)         # Makie wants Z as (x,y) => (C,R)

            h = heatmap!(ax, xcoords, ycoords, Z)
            hobj === nothing && (hobj = h)

            net_i < length(netfams) && (ax.xticklabelsvisible = false)
            ri > 1 && (ax.yticklabelsvisible = false)
        end
    end

    if hobj !== nothing
        Colorbar(gl[3:6, 6], hobj, label=string(metric))
    end

    if outfile !== nothing
        save(joinpath(OUTDIR, outfile), f)
    end

    return f
end

# ============================================================
# 12) Distribution figures for selected cells (hist + ECDF)
# ============================================================

function ecdf_xy(vals::Vector{Float64})
    isempty(vals) && return (Float64[], Float64[])
    v = sort(vals)
    n = length(v)
    y = range(1/n, 1.0, length=n)
    return v, collect(y)
end

function plot_distribution_case(mvec_all::Vector{Float64};
                                title::String="", outfile::Union{Nothing,String}=nothing)

    fig = Figure(size=(1200, 500))

    # Left: histogram
    ax1 = Axis(fig[1,1], title=title, xlabel="Per-species mismatch 1 - J(A,AB)", ylabel="Count")
    hist!(ax1, mvec_all, bins=30)
    vlines!(ax1, [TAIL_THRESH], linestyle=:dash, linewidth=3)

    # Right: ECDF
    ax2 = Axis(fig[1,2], title="ECDF", xlabel="Per-species mismatch", ylabel="F(x)")
    x, y = ecdf_xy(mvec_all)
    lines!(ax2, x, y, linewidth=2)
    vlines!(ax2, [TAIL_THRESH], linestyle=:dash, linewidth=3)

    # Footer with stats
    q90 = mismatch_q90(mvec_all)
    frac = mismatch_frac_gt(mvec_all, TAIL_THRESH)
    n = length(mvec_all)
    Label(fig[2,1:2],
        "n=$(n) | q90=$(round(q90, digits=3)) | frac(mismatch > $(TAIL_THRESH))=$(round(frac, digits=3))",
        fontsize=12
    )

    if outfile !== nothing
        save(joinpath(OUTDIR, outfile), fig)
    end
    return fig
end

# A small, sensible default set of distribution probes:
# - use autocorr env
# - show both ends of r and C + a mid-point
# - for each net family
# - for two regimes (one narrow-highvar, one broad-highvar)
const DIST_ENV = :autocorr
const DIST_REGIMES = [2, 4]  # indices in regimes: 2=Narrow+HighVar, 4=Broad+HighVar
const NREP_DIST = 25

function run_distribution_probes(store, Cvals, Rvals)
    # choose indices: low, mid, high
    cinds = [1, Int(cld(length(Cvals),2)), length(Cvals)]
    rinds = [1, Int(cld(length(Rvals),2)), length(Rvals)]

    for net in NETFAMS
        for ri in DIST_REGIMES
            reg = regimes[ri]
            for (cind, rind) in ((cinds[1], rinds[1]), (cinds[1], rinds[end]),
                                 (cinds[end], rinds[1]), (cinds[end], rinds[end]),
                                 (cinds[2], rinds[2]))

                C = Cvals[cind]
                r = Rvals[rind]

                # pool mismatches across replicates (and species)
                pooled = Float64[]
                for rep in 1:NREP_DIST
                    tid = Threads.threadid()
                    ws = WSS[tid]
                    rng = MersenneTwister(BASE_SEED + 900_000 + 31*rep + 7*tid)

                    seed = BASE_SEED +
                           77_000_000 * (findfirst(==(DIST_ENV), ENVKINDS)) +
                           9_000_000  * (findfirst(==(net), NETFAMS)) +
                           700_000    * ri +
                           70_000     * rind +
                           7_000      * cind +
                           rep
                    Random.seed!(rng, seed)

                    out = simulate_one!(rng, ws, DIST_ENV, net, reg, C, r; return_mvec=true)
                    append!(pooled, out.mvec)
                end

                ttl = "Dist probe — env=$(DIST_ENV), net=$(net), regime=$(reg.name)\nC=$(round(C,digits=4)), target_r=$(round(r,digits=2)), reps=$(NREP_DIST)"
                fname = "dist_env_$(DIST_ENV)_net_$(net)_reg$(ri)_C$(cind)_R$(rind).png"
                fig = plot_distribution_case(pooled; title=ttl, outfile=fname)
                display(fig)
            end
        end
    end
end

# ============================================================
# 13) MAIN
# ============================================================
println("OUTDIR: ", OUTDIR)

# --- Run sweep
store, Cvals, Rvals = sweep_all()

# --- Cache
cache_path = joinpath(OUTDIR, "sweep_cache_jaccard_tail.jls")
serialize(cache_path, (store=store, Cvals=Cvals, Rvals=Rvals))
println("Saved sweep cache to: ", cache_path)

using Serialization

data = deserialize(cache_path)
store = data.store
Cvals = data.Cvals
Rvals = data.Rvals
# --- Heatmaps
for env in ENVKINDS
    envname = env == :random ? "Random environment" : "Autocorrelated environment"

    f1 = facet_heatmaps(
        store, Cvals, Rvals, env, :mismatch_q90;
        title = "Jaccard mismatch tail: q90 of (1 - J(A,AB)) — $(envname)",
        outfile = "heatmaps_$(env)_mismatch_q90.png"
    )
    display(f1)

    f2 = facet_heatmaps(
        store, Cvals, Rvals, env, :mismatch_frac_gt;
        title = "Jaccard mismatch tail: frac(1 - J(A,AB) > $(TAIL_THRESH)) — $(envname)",
        outfile = "heatmaps_$(env)_mismatch_frac_gt.png"
    )
    display(f2)

    # Optional diagnostics (often useful to sanity-check)
    f3 = facet_heatmaps(
        store, Cvals, Rvals, env, :achieved_r;
        title = "Diagnostic: achieved mechanistic niche correlation — $(envname)",
        outfile = "heatmaps_$(env)_diagnostic_achieved_r.png"
    )
    display(f3)

    f4 = facet_heatmaps(
        store, Cvals, Rvals, env, :Creal;
        title = "Diagnostic: realized connectance L/S^2 — $(envname)",
        outfile = "heatmaps_$(env)_diagnostic_Creal.png"
    )
    display(f4)
end

# --- Distribution probes (hist + ECDF)
println("\nRunning distribution probes (hist + ECDF) ...")
run_distribution_probes(store, Cvals, Rvals)

println("\nDone. Output directory:")
println(OUTDIR)
