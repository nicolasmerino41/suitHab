#!/usr/bin/env julia
# Julia 1.11
# Run:  julia --threads auto fig1_relative_and_shared_10scenarios.jl
#
# PURPOSE
# Produce Figure-1-style panels for 10 scenarios under 3 habitat-loss geometries, using two comparison options:
#
#   OPTION A (relative-to-baseline): plot S(h)/S(0) for each line (A-only, AB, SAR-A, SAR-AB).
#   OPTION B (shared-cohort): restrict to species extant under BOTH A and AB at h=0 (given LCC≥Emin),
#                             then plot cohort richness under A and AB (plus cohort-anchored SAR-A and SAR-AB).
#
# OUTPUT
# - 10 figures total:
#     * 5 figures for OPTION A, each a 2×3 panel (two scenarios per figure × three geometries).
#     * 5 figures for OPTION B, each a 2×3 panel.
#
# NOTES
# - No external data needed.
# - “One-prey rule” is kept: in AB, a consumer is present in a cell if climate-suitable AND
#   at least one of its prey is present in that cell (OR across its prey set).
# - SAR-A is fit from baseline A-only occupancy; SAR-AB from baseline AB occupancy (both area-only).
# - In shared-cohort mode, SAR curves are fit/anchored using ONLY the cohort’s baseline occupancy.

using Random
using Statistics
using LinearAlgebra
using CairoMakie
using Printf
using Dates

# ============================================================
# 0) Global parameters (tune for speed/quality)
# ============================================================

const NX = 45
const NY = 45
const NCELLS = NX * NY

const S = 250
const BASAL_FRAC = 0.33
const Emin_patch = 80

# Habitat loss levels
const HL = collect(0.0:0.05:1.0)# Replicates (increase for smoother means)
const NREP = 10

# Climate
const CLIM_MIN = 0.0
const CLIM_MAX = 100.0
const CLIM_NOISE_SD = 7.0

# Niche model: 2D Gaussian suitability threshold
const SUIT_THRESH = 0.25

# Habitat-loss geometry params
const CLUSTER_SEED_PER_DESTROYED = 150
const FRONT_JITTER_MAX = 5
const FRONT_SMOOTH = :strong

# SAR fitting
const QUADRAT_SIZES = [3, 5, 9, 15, 45]  # divisors of 45

# Baseline viability control (prevents "already extinct" via LCC<Emin under A-only at h=0)
# Set to false if you want to allow baseline extinctions (cohort will shrink accordingly).
const ENFORCE_BASELINE_VIABILITY = true
const MAX_SIGMA_RESAMPLE = 30
const MAX_OCC_TRIES = 25

# Output
ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
OUTDIR = joinpath(pwd(), "fig1_compare_rel_vs_shared_" * ts)
isdir(OUTDIR) || mkpath(OUTDIR)

# Global seed
const BASE_SEED = 20260202

println("Threads: ", Threads.nthreads())
println("Grid: $(NX)x$(NY) cells=$(NCELLS)  Emin=$(Emin_patch)")
println("Saving to: ", OUTDIR)

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
# 2) Thread-safe RNG + workspace for BFS
# ============================================================

mutable struct CCWorkspace
    seen::Vector{Int32}
    stamp::Int32
    queue::Vector{Int}
end

function make_thread_rngs_and_workspaces(base_seed::Int)
    nt = Threads.nthreads()
    rngs = Vector{MersenneTwister}(undef, nt)
    wss  = Vector{CCWorkspace}(undef, nt)
    for t in 1:nt
        rngs[t] = MersenneTwister(base_seed + 10_000 * t + 137)
        wss[t]  = CCWorkspace(fill(Int32(0), NCELLS), Int32(0), Int[])
    end
    return rngs, wss
end
RNGS, WSS = make_thread_rngs_and_workspaces(BASE_SEED)

# ============================================================
# 3) Climate generation
# ============================================================

function make_climate(rng::AbstractRNG)
    env1 = Vector{Float64}(undef, NCELLS)
    env2 = Vector{Float64}(undef, NCELLS)
    for y in 1:NY, x in 1:NX
        i = linidx(x,y)
        g1 = (x - 1) / (NX - 1)
        g2 = (y - 1) / (NY - 1)
        env1[i] = CLIM_MIN + (CLIM_MAX - CLIM_MIN) * g1 + randn(rng)*CLIM_NOISE_SD
        env2[i] = CLIM_MIN + (CLIM_MAX - CLIM_MIN) * g2 + randn(rng)*CLIM_NOISE_SD
    end
    @inbounds for i in 1:NCELLS
        env1[i] = clamp(env1[i], CLIM_MIN, CLIM_MAX)
        env2[i] = clamp(env2[i], CLIM_MIN, CLIM_MAX)
    end
    return env1, env2
end

# ============================================================
# 4) Habitat loss orders
# ============================================================

function order_random(rng::AbstractRNG)
    randperm(rng, NCELLS)
end

function order_clustered(rng::AbstractRNG)
    nseed_total = max(1, round(Int, NCELLS / CLUSTER_SEED_PER_DESTROYED))
    seeds = randperm(rng, NCELLS)[1:nseed_total]
    destroyed = BitVector(falses(NCELLS))
    frontier = Int[]
    order = Int[]
    for s in seeds
        destroyed[s] = true
        push!(frontier, s)
        push!(order, s)
    end
    while length(order) < NCELLS
        cand = Int[]
        for v in frontier
            for nb in NEIGH_4[v]
                if !destroyed[nb]
                    push!(cand, nb)
                end
            end
        end
        if isempty(cand)
            for i in 1:NCELLS
                if !destroyed[i]
                    push!(cand, i)
                end
            end
        end
        shuffle!(rng, cand)
        take_n = min(length(cand), 40)
        newfront = Int[]
        for k in 1:take_n
            c = cand[k]
            if !destroyed[c]
                destroyed[c] = true
                push!(order, c)
                push!(newfront, c)
                length(order) == NCELLS && break
            end
        end
        frontier = newfront
    end
    return order
end

function smooth1d!(v::Vector{Float64}, strength::Symbol)
    strength == :none && return v
    w = strength == :mild ? 3 : 7
    tmp = copy(v)
    for i in 1:length(v)
        lo = max(1, i-w)
        hi = min(length(v), i+w)
        v[i] = mean(@view tmp[lo:hi])
    end
    return v
end

function order_front(rng::AbstractRNG)
    jit = [rand(rng)*2*FRONT_JITTER_MAX - FRONT_JITTER_MAX for _ in 1:NY]
    smooth1d!(jit, FRONT_SMOOTH)
    scores = Vector{Tuple{Float64,Int}}(undef, NCELLS)
    k = 1
    for y in 1:NY, x in 1:NX
        i = linidx(x,y)
        s = x + jit[y] + 0.03*randn(rng)
        scores[k] = (s, i)
        k += 1
    end
    sort!(scores, by=first)
    return [scores[i][2] for i in 1:NCELLS]
end

function habitat_mask_from_order(order::Vector{Int}, k_destroy::Int)
    m = BitVector(trues(NCELLS))
    k_destroy = clamp(k_destroy, 0, NCELLS)
    @inbounds for i in 1:k_destroy
        m[order[i]] = false
    end
    return m
end

# ============================================================
# 5) Connectivity: LCC≥Emin
# ============================================================

function lcc_size_ge!(ws::CCWorkspace, mask::BitVector, Emin::Int)
    if count(mask) < Emin
        return 0, false
    end
    ws.stamp += 1
    stamp = ws.stamp
    seen = ws.seen
    q = ws.queue
    empty!(q)

    best = 0
    @inbounds for i in 1:NCELLS
        if mask[i] && seen[i] != stamp
            seen[i] = stamp
            push!(q, i)
            qpos = 1
            compsize = 0
            while qpos <= length(q)
                v = q[qpos]; qpos += 1
                compsize += 1
                for nb in NEIGH_4[v]
                    if mask[nb] && seen[nb] != stamp
                        seen[nb] = stamp
                        push!(q, nb)
                    end
                end
            end
            best = max(best, compsize)
            best >= Emin && return best, true
            empty!(q)
        end
    end
    return best, best >= Emin
end

# ============================================================
# 6) Metaweb + niche-correlation control
# ============================================================

# Build a directed metaweb with target connectance C := L / S^2
# Basal species have no prey; consumers can have many prey (still OR-rule in AB).
function make_metaweb(rng::AbstractRNG, basal_mask::BitVector, C::Float64)
    Ltarget = round(Int, C * S^2)
    prey = [Int[] for _ in 1:S]

    candidates = Tuple{Int,Int}[]
    sizehint!(candidates, S*S)
    for i in 1:S
        if basal_mask[i]; continue; end
        for j in 1:S
            i == j && continue
            push!(candidates, (i,j))
        end
    end
    L = min(Ltarget, length(candidates))
    perm = randperm(rng, length(candidates))
    for k in 1:L
        (i,j) = candidates[perm[k]]
        push!(prey[i], j)
    end
    return prey
end

function pearson_r(a::Vector{Float64}, b::Vector{Float64})
    ma = mean(a); mb = mean(b)
    sa = std(a);  sb = std(b)
    (sa == 0 || sb == 0) && return 0.0
    return mean((a .- ma) .* (b .- mb)) / (sa * sb)
end

function compute_prey_means(mu1::Vector{Float64}, mu2::Vector{Float64}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    pm1 = similar(mu1)
    pm2 = similar(mu2)
    for i in 1:S
        if basal_mask[i] || isempty(prey[i])
            pm1[i] = NaN
            pm2[i] = NaN
        else
            ps = prey[i]
            pm1[i] = mean(mu1[ps])
            pm2[i] = mean(mu2[ps])
        end
    end
    return pm1, pm2
end

function centroid_corr(mu1::Vector{Float64}, mu2::Vector{Float64}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    pm1, pm2 = compute_prey_means(mu1, mu2, prey, basal_mask)
    idx = findall(i -> !(basal_mask[i] || isempty(prey[i])) && isfinite(pm1[i]) && isfinite(pm2[i]), 1:S)
    length(idx) < 5 && return 0.0
    r1 = pearson_r(mu1[idx], pm1[idx])
    r2 = pearson_r(mu2[idx], pm2[idx])
    return 0.5*(r1 + r2)
end

function assign_centroids_with_target_corr(rng::AbstractRNG, prey::Vector{Vector{Int}}, basal_mask::BitVector, target_r::Float64;
                                          relax_iters::Int=30)
    mu1 = rand(rng, S) .* (CLIM_MAX - CLIM_MIN) .+ CLIM_MIN
    mu2 = rand(rng, S) .* (CLIM_MAX - CLIM_MIN) .+ CLIM_MIN

    function relax(alpha::Float64)
        m1 = copy(mu1)
        m2 = copy(mu2)
        for _ in 1:relax_iters
            pm1, pm2 = compute_prey_means(m1, m2, prey, basal_mask)
            @inbounds for i in 1:S
                if basal_mask[i] || isempty(prey[i]); continue; end
                m1[i] = clamp((1-alpha)*m1[i] + alpha*pm1[i] + randn(rng)*1.0, CLIM_MIN, CLIM_MAX)
                m2[i] = clamp((1-alpha)*m2[i] + alpha*pm2[i] + randn(rng)*1.0, CLIM_MIN, CLIM_MAX)
            end
        end
        return centroid_corr(m1, m2, prey, basal_mask), m1, m2
    end

    lo, hi = -0.98, 0.98
    best_err = Inf
    best = (0.0, copy(mu1), copy(mu2))
    for _ in 1:24
        mid = 0.5*(lo + hi)
        rmid, m1, m2 = relax(mid)
        err = abs(rmid - target_r)
        if err < best_err
            best_err = err
            best = (rmid, m1, m2)
        end
        if rmid < target_r
            lo = mid
        else
            hi = mid
        end
    end
    rgot, m1, m2 = best
    return m1, m2, rgot
end

# ============================================================
# 7) Niche model + occupancy-targeted sigma assignment
# ============================================================

@inline function suitability_mask(env1, env2, mu1, mu2, s1, s2, thresh)
    lim = -2.0 * log(thresh)
    m = BitVector(undef, NCELLS)
    invs1 = 1.0/s1
    invs2 = 1.0/s2
    @inbounds for i in 1:NCELLS
        d1 = (env1[i] - mu1) * invs1
        d2 = (env2[i] - mu2) * invs2
        m[i] = (d1*d1 + d2*d2) <= lim
    end
    return m
end

# Gamma sampler (for Beta via ratio-of-gammas)
function randgamma(rng::AbstractRNG, k::Float64)
    if k < 1
        return randgamma(rng, k+1) * rand(rng)^(1/k)
    end
    d = k - 1/3
    c = 1 / sqrt(9d)
    while true
        x = randn(rng)
        v = (1 + c*x)^3
        if v <= 0; continue; end
        u = rand(rng)
        if u < 1 - 0.0331*(x^4)
            return d*v
        end
        if log(u) < 0.5*x^2 + d*(1 - v + log(v))
            return d*v
        end
    end
end

struct BetaLike
    a::Float64
    b::Float64
end
import Base: rand
function rand(rng::AbstractRNG, d::BetaLike)
    x = randgamma(rng, d.a)
    y = randgamma(rng, d.b)
    return x / (x + y)
end

# Given a target occupancy fraction p_target and anisotropy factor aniso, solve for (s1,s2)
function solve_sigmas_for_occ(env1, env2, mu1, mu2, aniso, p_target; iters=16)
    p_target = clamp(p_target, 0.002, 0.995)
    s_lo, s_hi = 1.2, 80.0
    best_s, best_err = 12.0, Inf
    for _ in 1:iters
        s_mid = 0.5*(s_lo + s_hi)
        s1 = max(0.6, s_mid * aniso)
        s2 = max(0.6, s_mid / aniso)
        mask = suitability_mask(env1, env2, mu1, mu2, s1, s2, SUIT_THRESH)
        p = count(mask) / NCELLS
        err = abs(p - p_target)
        if err < best_err
            best_err = err
            best_s = s_mid
        end
        if p < p_target
            s_lo = s_mid
        else
            s_hi = s_mid
        end
    end
    s1 = max(0.6, best_s * aniso)
    s2 = max(0.6, best_s / aniso)
    return s1, s2
end

# ============================================================
# 8) AB fixed point (one-prey rule retained)
# ============================================================

function fixed_point_AB(A_presence::Vector{BitVector}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    pres = [copy(A_presence[i]) for i in 1:S]
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
                newp[i] = A_presence[i] .& u
            end
        end
        for i in 1:S
            if newp[i] != pres[i]
                pres[i] = newp[i]
                changed = true
            end
        end
        iter > 60 && break
    end
    return pres
end

# ============================================================
# 9) Extant mask + richness for a given set (optionally restricted to a cohort)
# ============================================================
function extant_mask!(
    ws::CCWorkspace,
    pres::Vector{BitVector},
    basal_mask::BitVector;
    which::Symbol = :consumers,
    cohort::Union{BitVector,Nothing} = nothing
)
    ext = BitVector(falses(S))
    @inbounds for sp in 1:S
        if which == :consumers && basal_mask[sp]
            continue
        end
        if cohort !== nothing && !cohort[sp]
            continue
        end
        _, ok = lcc_size_ge!(ws, pres[sp], Emin_patch)
        ext[sp] = ok
    end
    return ext
end

function gamma_richness_from_ext(
    ext::BitVector,
    basal_mask::BitVector;
    which::Symbol = :consumers,
    cohort::Union{BitVector,Nothing} = nothing
)
    c = 0
    @inbounds for sp in 1:S
        if which == :consumers && basal_mask[sp]
            continue
        end
        if cohort !== nothing && !cohort[sp]
            continue
        end
        c += ext[sp] ? 1 : 0
    end
    return c
end

# ============================================================
# 10) SAR: fit z from non-overlapping quadrats for a specified species subset
# ============================================================
function sar_points_gamma(
    pres::Vector{BitVector},
    basal_mask::BitVector;
    which::Symbol = :consumers,
    cohort::Union{BitVector,Nothing} = nothing
)
    A = Float64[]
    Y = Float64[]

    for q in QUADRAT_SIZES
        nx = NX ÷ q
        ny = NY ÷ q
        vals = Float64[]

        for by in 0:(ny-1), bx in 0:(nx-1)
            inds = Int[]
            for yy in 1:q, xx in 1:q
                push!(inds, linidx(bx*q + xx, by*q + yy))
            end

            c = 0
            for sp in 1:S
                if which == :consumers && basal_mask[sp]
                    continue
                end
                if cohort !== nothing && !cohort[sp]
                    continue
                end

                m = pres[sp]
                present = false
                for i in inds
                    if m[i]
                        present = true
                        break
                    end
                end
                c += present ? 1 : 0
            end

            push!(vals, c)
        end

        push!(A, q^2)
        push!(Y, mean(vals))
    end

    return A, Y
end

function fit_z_only(A::Vector{Float64}, Y::Vector{Float64})
    idx = findall(i -> (A[i] > 0 && Y[i] > 0), 1:length(A))
    length(idx) < 2 && return 0.2  # fallback
    x = log.(A[idx])
    y = log.(Y[idx])
    X = hcat(ones(length(x)), x)
    β = X \ y
    return clamp(β[2], 0.01, 0.99)
end

@inline function sar_predict(S0::Float64, A0::Float64, z::Float64, Arem::Float64)
    Arem = clamp(Arem, 1.0, A0)
    return S0 * (Arem / A0)^z
end

# ============================================================
# 11) Scenario design (10 scenarios ordered from "highest divergence" to "lowest")
# ============================================================

struct Scenario
    id::Int
    name::String
    C::Float64
    target_r::Float64
    # occupancy distribution controls:
    occ_mode::Symbol   # :varwide, :varwide_mildsigma, :skew_narrow, :skew_broad, :mid
    sigma_var::Symbol  # :high, :mild, :low
end

# You can edit these to match *exactly* your “past 10 scenarios”.
const SCENARIOS10 = Scenario[
    Scenario(1,  "MaxDiv 1", 0.03, -0.40, :skew_narrow,       :high),
    Scenario(2,  "MaxDiv 2", 0.05, -0.25, :skew_narrow,       :high),
    Scenario(3,  "HighDiv 3",0.07, -0.10, :varwide,           :high),
    Scenario(4,  "HighDiv 4",0.09,  0.05, :varwide,           :mild),
    Scenario(5,  "MidDiv 5", 0.11,  0.20, :varwide_mildsigma, :mild),
    Scenario(6,  "MidDiv 6", 0.13,  0.35, :mid,               :mild),
    Scenario(7,  "LowDiv 7", 0.15,  0.50, :mid,               :low),
    Scenario(8,  "LowDiv 8", 0.18,  0.65, :skew_broad,        :low),
    Scenario(9,  "MinDiv 9", 0.21,  0.78, :skew_broad,        :low),
    Scenario(10, "MinDiv 10",0.24,  0.88, :skew_broad,        :low),
]

# Occupancy target sampler (fraction of grid climatically suitable)
function draw_target_occ(rng::AbstractRNG, scen::Scenario)
    if scen.occ_mode == :varwide
        return 0.03 + rand(rng) * (0.75 - 0.03)
    elseif scen.occ_mode == :varwide_mildsigma
        return 0.05 + rand(rng) * (0.65 - 0.05)
    elseif scen.occ_mode == :skew_narrow
        x = rand(rng, BetaLike(2.0, 8.0))
        return 0.02 + x * (0.55 - 0.02)
    elseif scen.occ_mode == :skew_broad
        x = rand(rng, BetaLike(8.0, 2.0))
        return 0.18 + x * (0.92 - 0.18)
    else # :mid
        x = rand(rng, BetaLike(3.0, 3.0))
        return 0.08 + x * (0.70 - 0.08)
    end
end

function draw_anisotropy(rng::AbstractRNG, sigma_var::Symbol)
    if sigma_var == :high
        return exp(randn(rng) * 0.65)
    elseif sigma_var == :mild
        return exp(randn(rng) * 0.30)
    else
        return exp(randn(rng) * 0.15)
    end
end

# Assign per-species climate masks to match occupancy targets, and (optionally) enforce LCC≥Emin at h=0.
function assign_climate_masks(rng::AbstractRNG, ws::CCWorkspace, env1, env2, mu1, mu2, scen::Scenario)
    clim = Vector{BitVector}(undef, S)
    for sp in 1:S
        got = false
        # attempt by adjusting occupancy target upward if LCC fails
        p0 = draw_target_occ(rng, scen)
        for _ in 1:MAX_OCC_TRIES
            p = p0
            for _ in 1:MAX_SIGMA_RESAMPLE
                an = draw_anisotropy(rng, scen.sigma_var)
                s1, s2 = solve_sigmas_for_occ(env1, env2, mu1[sp], mu2[sp], an, p)
                m = suitability_mask(env1, env2, mu1[sp], mu2[sp], s1, s2, SUIT_THRESH)
                if ENFORCE_BASELINE_VIABILITY
                    _, ok = lcc_size_ge!(ws, m, Emin_patch)
                    if ok
                        clim[sp] = m
                        got = true
                        break
                    end
                else
                    clim[sp] = m
                    got = true
                    break
                end
            end
            got && break
            # if we’re enforcing viability and failed, increase target occupancy a bit and retry
            p0 = clamp(p0 * 1.12, 0.01, 0.98)
        end
        if !got
            # last resort: very broad niche
            clim[sp] = trues(NCELLS)
        end
    end
    return clim
end

# ============================================================
# 12) One replicate simulation: returns curves for
#   - relative-to-baseline (normalized per line)
#   - shared-cohort (raw cohort richness)
# ============================================================

function simulate_replicate(rng::AbstractRNG, ws::CCWorkspace, scen::Scenario, geometry::Symbol; which::Symbol=:consumers)

    nb = round(Int, BASAL_FRAC * S)
    basal_mask = BitVector(falses(S))
    basal_mask[1:nb] .= true

    env1, env2 = make_climate(rng)
    prey = make_metaweb(rng, basal_mask, scen.C)
    mu1, mu2, rgot = assign_centroids_with_target_corr(rng, prey, basal_mask, scen.target_r)

    clim = assign_climate_masks(rng, ws, env1, env2, mu1, mu2, scen)

    order = geometry == :random  ? order_random(rng) :
            geometry == :cluster ? order_clustered(rng) :
            geometry == :front   ? order_front(rng) :
            error("Unknown geometry")

    habitat0 = BitVector(trues(NCELLS))
    A0  = [clim[sp] .& habitat0 for sp in 1:S]
    AB0 = fixed_point_AB(A0, prey, basal_mask)

    # baseline extant sets (connectivity rule)
    extA0  = extant_mask!(ws, A0,  basal_mask; which=which)
    extAB0 = extant_mask!(ws, AB0, basal_mask; which=which)

    S0A  = float(gamma_richness_from_ext(extA0,  basal_mask; which=which))
    S0AB = float(gamma_richness_from_ext(extAB0, basal_mask; which=which))

    # Shared cohort: species extant in both at baseline
    cohort = extA0 .& extAB0
    S0C = float(gamma_richness_from_ext(cohort, basal_mask; which=which, cohort=cohort))  # equals count(cohort within which)

    # SAR fits (area-only) from baseline occupancy
    AptsA,  SptsA  = sar_points_gamma(A0,  basal_mask; which=which)
    zA = fit_z_only(AptsA, SptsA)

    AptsAB, SptsAB = sar_points_gamma(AB0, basal_mask; which=which)
    zAB = fit_z_only(AptsAB, SptsAB)

    # SAR fits for the cohort only (for shared-cohort panels)
    AptsAC,  SptsAC  = sar_points_gamma(A0,  basal_mask; which=which, cohort=cohort)
    zAC = fit_z_only(AptsAC, SptsAC)

    AptsABC, SptsABC = sar_points_gamma(AB0, basal_mask; which=which, cohort=cohort)
    zABC = fit_z_only(AptsABC, SptsABC)

    A0_total = float(NCELLS)
    nH = length(HL)

    # raw curves
    yA  = zeros(Float64, nH)
    yAB = zeros(Float64, nH)
    ySAR_A  = zeros(Float64, nH)
    ySAR_AB = zeros(Float64, nH)

    # shared-cohort curves
    yA_C  = zeros(Float64, nH)
    yAB_C = zeros(Float64, nH)
    ySAR_A_C  = zeros(Float64, nH)
    ySAR_AB_C = zeros(Float64, nH)

    for (k,h) in enumerate(HL)
        kdestroy = round(Int, h * NCELLS)
        habitat = habitat_mask_from_order(order, kdestroy)

        A_pres  = [clim[sp] .& habitat for sp in 1:S]
        AB_pres = fixed_point_AB(A_pres, prey, basal_mask)

        extA  = extant_mask!(ws, A_pres,  basal_mask; which=which)
        extAB = extant_mask!(ws, AB_pres, basal_mask; which=which)

        yA[k]  = gamma_richness_from_ext(extA,  basal_mask; which=which)
        yAB[k] = gamma_richness_from_ext(extAB, basal_mask; which=which)

        # area-only SAR predictions (anchored to each baseline richness)
        Arem = (1.0 - h) * A0_total
        ySAR_A[k]  = sar_predict(S0A,  A0_total, zA,  Arem)
        ySAR_AB[k] = sar_predict(S0AB, A0_total, zAB, Arem)  # “trophic-informed baseline at h=0”

        # Shared cohort: restrict evaluation to baseline-shared species
        extA_coh  = extant_mask!(ws, A_pres,  basal_mask; which=which, cohort=cohort)
        extAB_coh = extant_mask!(ws, AB_pres, basal_mask; which=which, cohort=cohort)

        yA_C[k]  = gamma_richness_from_ext(extA_coh,  basal_mask; which=which, cohort=cohort)
        yAB_C[k] = gamma_richness_from_ext(extAB_coh, basal_mask; which=which, cohort=cohort)

        # cohort-anchored SAR curves
        ySAR_A_C[k]  = sar_predict(S0C, A0_total, zAC,  Arem)
        ySAR_AB_C[k] = sar_predict(S0C, A0_total, zABC, Arem)
    end

    # OPTION A: relative-to-own-baseline for each line
    rel = Dict(
        :A  => (S0A  > 0 ? yA  ./ S0A  : fill(0.0, nH)),
        :AB => (S0AB > 0 ? yAB ./ S0AB : fill(0.0, nH)),
        :SAR_A  => (S0A  > 0 ? ySAR_A  ./ S0A  : fill(0.0, nH)),
        :SAR_AB => (S0AB > 0 ? ySAR_AB ./ S0AB : fill(0.0, nH)),
        :S0A => S0A, :S0AB => S0AB, :rgot => rgot
    )

    # OPTION B: shared cohort (raw)
    shared = Dict(
        :A  => yA_C,
        :AB => yAB_C,
        :SAR_A  => ySAR_A_C,
        :SAR_AB => ySAR_AB_C,
        :S0C => S0C, :cohort_count => count(cohort), :rgot => rgot
    )

    return rel, shared
end

# ============================================================
# 13) Averaging over replicates
# ============================================================
function mean_over_reps(
    dicts::Vector{Dict{Symbol,Vector{Float64}}},
    keys::Vector{Symbol}
)
    out = Dict{Symbol,Vector{Float64}}()

    for k in keys
        mats = reduce(hcat, (d[k] for d in dicts))
        out[k] = vec(mean(mats; dims=2))
    end

    return out
end

# ============================================================
# 14) Plotting helpers
# ============================================================

function add_lines_relative!(ax, x, d)
    lines!(ax, x, d[:A],  label="Mechanistic A")
    lines!(ax, x, d[:AB], label="Mechanistic AB")
    lines!(ax, x, d[:SAR_A],  linestyle=:dash, label="SAR (A baseline)")
    lines!(ax, x, d[:SAR_AB], linestyle=:dash, label="SAR (AB baseline)")
end

function add_lines_shared!(ax, x, d)
    lines!(ax, x, d[:A],  label="Mechanistic A (shared)")
    lines!(ax, x, d[:AB], label="Mechanistic AB (shared)")
    lines!(ax, x, d[:SAR_A],  linestyle=:dash, label="SAR (A, cohort-fit)")
    lines!(ax, x, d[:SAR_AB], linestyle=:dash, label="SAR (AB, cohort-fit)")
end

function scenario_title(scen::Scenario)
    return @sprintf("%s | C=%.2f, r*=%.2f, occ=%s, σvar=%s",
        scen.name, scen.C, scen.target_r, String(scen.occ_mode), String(scen.sigma_var))
end

function make_panel_2x3(
    results_by_geom::Dict{Symbol,Dict{Symbol,Dict{Symbol,Vector{Float64}}}},
    scenA::Scenario,
    scenB::Scenario;
    mode::Symbol = :relative
)
    f = Figure(size=(1700, 900))
    geoms = [:random, :cluster, :front]
    geom_titles = Dict(:random=>"Random loss",
                       :cluster=>"Clustered loss",
                       :front=>"Front loss")

    scens = [scenA, scenB]

    top_label = mode == :relative ?
        "OPTION A — Relative curves (S(h)/S(0) for each line)" :
        "OPTION B — Shared cohort (species extant in both A and AB at h=0)"
    Label(f[0, :], top_label, fontsize=20)

    for (row, scen) in enumerate(scens)
        for (col, g) in enumerate(geoms)
            ax = Axis(
                f[row, col],
                title  = geom_titles[g],
                xlabel = "Habitat loss (proportion destroyed)",
                ylabel = mode == :relative ?
                    "Relative richness (S/S0)" :
                    "Richness (shared cohort)"
            )

            d = results_by_geom[g][Symbol("scen$(scen.id)")]

            if mode == :relative
                add_lines_relative!(ax, HL, d)
                ylims!(ax, 0, 1.05)
            else
                add_lines_shared!(ax, HL, d)
            end

            if row == 1 && col == 3
                axislegend(ax; position=:rb, framevisible=false)
            end
        end

        Label(
            f[row, 0],
            scenario_title(scen),
            fontsize = 13,
            rotation = pi/2,
            padding  = (0, 10, 0, 0)
        )
    end

    return f
end

# ============================================================
# 15) Main runner: compute 10 scenarios × 3 geoms × NREP, then build 10 figures
# ============================================================
function compute_means_for_scenarios(mode::Symbol)
    geoms = [:random, :cluster, :front]
    which = :consumers

    # geom => scen_id => Dict(curve => Vector{Float64})
    out = Dict{Symbol,Dict{Symbol,Dict{Symbol,Vector{Float64}}}}()

    curve_keys = (:A, :AB, :SAR_A, :SAR_AB)

    for g in geoms
        out[g] = Dict{Symbol,Dict{Symbol,Vector{Float64}}}()

        for scen in SCENARIOS10
            reps = Vector{Dict{Symbol,Vector{Float64}}}(undef, NREP)

            Threads.@threads for r in 1:NREP
                tid = Threads.threadid()
                rng = RNGS[tid]
                ws  = WSS[tid]

                seed = BASE_SEED +
                       1_000_000 * scen.id +
                       50_000 * (findfirst(==(g), geoms)) +
                       10_000 * (mode == :relative ? 1 : 2) +
                       r

                Random.seed!(rng, seed)

                rel, shared = simulate_replicate(rng, ws, scen, g; which=which)
                d = (mode == :relative) ? rel : shared

                # keep ONLY vector-valued curve entries
                reps[r] = Dict(k => d[k] for k in curve_keys)
            end

            mean_curves = mean_over_reps(reps, collect(curve_keys))
            out[g][Symbol("scen$(scen.id)")] = mean_curves

            println("Computed $(mode) means: scen=$(scen.id) geom=$(g)")
        end
    end

    return out
end

# Build the 5 pairs (two scenarios per figure, in order)
const PAIRS = [(SCENARIOS10[2i-1], SCENARIOS10[2i]) for i in 1:5]

println("\nComputing OPTION A (relative curves) ...")
means_relative = compute_means_for_scenarios(:relative)

println("\nComputing OPTION B (shared cohort) ...")
means_shared = compute_means_for_scenarios(:shared)

# Save figures: 5 for each option
for (i,(a,b)) in enumerate(PAIRS)
    # OPTION A
    figA = make_panel_2x3(means_relative, a, b; mode=:relative)
    save(joinpath(OUTDIR, @sprintf("fig_rel_pair%02d_scen%02d_%02d.png", i, a.id, b.id)), figA)
    display(figA)

    # OPTION B
    figB = make_panel_2x3(means_shared, a, b; mode=:shared)
    save(joinpath(OUTDIR, @sprintf("fig_shared_pair%02d_scen%02d_%02d.png", i, a.id, b.id)), figB)
    display(figB)
end

println("\nDone. Output directory:\n", OUTDIR)
