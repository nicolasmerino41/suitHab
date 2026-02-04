#!/usr/bin/env julia
# Julia 1.11
# Run:  julia --threads auto traffic_extinction_figures_FINAL.jl
#
# Self-contained, no external data. Produces:
#   - Figure 1A: Gamma richness (all species) — 2×3 panels (High vs Low divergence × geometry)
#   - Figure 1B: Gamma richness (consumers only) — same layout
#   - Figure 2: Heatmaps of AUC divergence (consumers only) across Connectance × Niche correlation
#       * outputs BOTH raw AUC and relative AUC (scaled 0–1 per scenario & geometry)
#   - Figure 3 (mechanism cartoon-like plot): support amount vs support connectivity vs fragmentation-fail
#
# Key model choices (to ensure divergence is real):
#   1) AB is trophically constrained at patch level (fixed point of prey support).
#   2) AB includes an extinction cascade: species with LCC < Emin are removed, which can remove prey and
#      cascade to consumers (this is what makes A ≠ AB clearly in the right parameter regimes).
#   3) SAR_baseline is pure area-only using remaining area (1-h).
#   4) SAR_effective is pure area-only but uses an EFFECTIVE AREA = LCC area of the trophically supported union mask
#      computed from the AB fixed-point occupancy BEFORE extinction filtering (non-circular, geometry-sensitive).
#   5) Fig1 uses two explicit regimes:
#        - High divergence: low connectance + negative corr + skew-narrow niches (prey narrower)
#        - Low divergence:  high connectance + high corr     + skew-broad niches
#
# Notes:
#   - This is intentionally “article-ready”: consistent fonts, titles, saved PNGs, deterministic seeds.
#   - If you want smoother Fig2 transitions, increase N_CONNECT/N_CORR and/or NREP_HEAT (runtime increases).
#
using Random
using Statistics
using LinearAlgebra
using CairoMakie
using Printf
using Dates

# ============================================================
# 0) Global parameters
# ============================================================

# Grid
const NX = 45
const NY = 45
const NCELLS = NX * NY

# Regional pool / trophic structure
const S = 250
const BASAL_FRAC = 0.33

# Extinction rule: species survives if its largest connected component (LCC) ≥ Emin_patch
const Emin_patch = 80

# Habitat loss levels
const HL = collect(0.0:0.05:1.0)

# Replicates
const NREP_FIG1 = 10          # Fig1 replicate averaging
const NREP_HEAT = 8           # Fig2 replicate averaging per cell (increase for smoother maps)

# Heatmap axes (connectance & correlation)
const CONNECTANCE_RANGE = (0.02, 0.25)
const CORR_RANGE        = (-0.5, 0.9)
const N_CONNECT = 21          # finer grid for smoother transitions
const N_CORR    = 21

# Climate field (two environmental axes)
const CLIM_MIN = 0.0
const CLIM_MAX = 100.0
const CLIM_NOISE_SD = 7.0

# Niche definition: 2D Gaussian with threshold
const SUIT_THRESH = 0.25

# Habitat-loss geometries
const CLUSTER_SEED_PER_DESTROYED = 150
const FRONT_JITTER_MAX = 5
const FRONT_SMOOTH = :strong  # :none | :mild | :strong

# Baseline viability target under A-only at h=0 (resample/inflate sigmas to avoid trivial baseline extinctions)
const BASELINE_EXTANT_TARGET = 0.95
const MAX_SIGMA_RESAMPLE = 30

# Monotonic sanity check / enforcement for plotted curves
const ENFORCE_MONOTONE = true
const MONO_TOL = 1e-9

# Output
ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
OUTDIR = joinpath(pwd(), "figures_traffic_extinction_FINAL_" * ts)
isdir(OUTDIR) || mkpath(OUTDIR)

# Global seed
const BASE_SEED = 20250128

# ============================================================
# 1) Index helpers + 4-neighbour adjacency
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
# 2) Thread-safe RNG + workspace for connectivity
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
rngs, wss = make_thread_rngs_and_workspaces(BASE_SEED)

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
# 4) Metaweb (directed). Basal have no prey. No cannibalism.
#    Connectance definition: C = L / S^2 (as in your code).
# ============================================================

function make_metaweb(rng::AbstractRNG, basal_mask::BitVector, C::Float64)
    Ltarget = round(Int, C * S^2)
    prey = [Int[] for _ in 1:S]

    candidates = Vector{Tuple{Int,Int}}()
    sizehint!(candidates, S*S)
    for i in 1:S
        if basal_mask[i]; continue; end
        for j in 1:S
            if i == j; continue; end
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

# ============================================================
# 5) Consumer–prey centroid correlation targeting
# ============================================================

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

function assign_centroids_with_target_corr(
    rng::AbstractRNG,
    prey::Vector{Vector{Int}},
    basal_mask::BitVector,
    target_r::Float64;
    relax_iters::Int = 30
)
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
# 6) Niche: 2D Gaussian threshold mask
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

# --- Fig2-style occupancy distributions (also reused for Fig1 sigmas) ---

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

import Base: rand
struct BetaLike
    a::Float64
    b::Float64
end
function rand(rng::AbstractRNG, d::BetaLike)
    x = randgamma(rng, d.a)
    y = randgamma(rng, d.b)
    return x / (x + y)
end

abstract type NicheScenario end
struct VarOccVarSigma <: NicheScenario end
struct VarOccMildSigma <: NicheScenario end
struct SkewNarrow      <: NicheScenario end
struct SkewBroad       <: NicheScenario end

const FIG2_SCENARIOS = (VarOccVarSigma(), VarOccMildSigma(), SkewNarrow(), SkewBroad())
const FIG2_SCENARIO_NAMES = ("VarOcc + VarSigma", "VarOcc + MildSigma", "SkewNarrow", "SkewBroad")

function draw_target_occ(rng::AbstractRNG, scen::NicheScenario)
    if scen isa VarOccVarSigma
        return 0.04 + rand(rng) * (0.75 - 0.04)
    elseif scen isa VarOccMildSigma
        return 0.04 + rand(rng) * (0.75 - 0.04)
    elseif scen isa SkewNarrow
        x = rand(rng, BetaLike(2.0, 8.0))
        return 0.02 + x * (0.55 - 0.02)
    else # SkewBroad
        x = rand(rng, BetaLike(8.0, 2.0))
        return 0.18 + x * (0.92 - 0.18)
    end
end

function solve_sigmas_for_occ(env1, env2, mu1, mu2, aniso, p_target; iters=18)
    p_target = clamp(p_target, 0.002, 0.995)
    s_lo, s_hi = 1.5, 70.0
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

# Scenario-specific anisotropy intensity (controls sigma variability)
function draw_anisotropy(rng::AbstractRNG, scen::NicheScenario)
    if scen isa VarOccVarSigma
        return exp(randn(rng) * 0.65)
    elseif scen isa VarOccMildSigma
        return exp(randn(rng) * 0.25)
    elseif scen isa SkewNarrow
        return exp(randn(rng) * 0.40)
    else # SkewBroad
        return exp(randn(rng) * 0.30)
    end
end

# Assign sigmas while preventing trivial baseline extinctions under A-only at h=0:
# resample until LCC≥Emin, else inflate as last resort.
#
# Additionally, for "high divergence" we want PREY narrower than consumers.
# We implement that by shrinking sigma for basal species by PREY_SIGMA_FACTOR (and then ensuring viability).
const PREY_SIGMA_FACTOR_HIGH = 0.65

function assign_sigmas_with_viability(
    rng::AbstractRNG,
    ws::CCWorkspace,
    env1, env2,
    mu1::Vector{Float64}, mu2::Vector{Float64},
    basal_mask::BitVector,
    scen::NicheScenario;
    prey_sigma_factor::Float64 = 1.0
)
    s1 = Vector{Float64}(undef, S)
    s2 = Vector{Float64}(undef, S)
    clim = Vector{BitVector}(undef, S)

    for sp in 1:S
        got = false
        for _ in 1:MAX_SIGMA_RESAMPLE
            p = draw_target_occ(rng, scen)
            an = draw_anisotropy(rng, scen)
            a,b = solve_sigmas_for_occ(env1, env2, mu1[sp], mu2[sp], an, p)

            if basal_mask[sp]
                a *= prey_sigma_factor
                b *= prey_sigma_factor
                a = max(a, 0.6); b = max(b, 0.6)
            end

            m = suitability_mask(env1, env2, mu1[sp], mu2[sp], a, b, SUIT_THRESH)
            _, ok = lcc_size_ge!(ws, m, Emin_patch)
            if ok
                s1[sp] = a; s2[sp] = b; clim[sp] = m
                got = true
                break
            end
        end

        if !got
            # inflate until viable
            p = draw_target_occ(rng, scen)
            an = draw_anisotropy(rng, scen)
            a,b = solve_sigmas_for_occ(env1, env2, mu1[sp], mu2[sp], an, p)
            if basal_mask[sp]
                a *= prey_sigma_factor
                b *= prey_sigma_factor
                a = max(a, 0.6); b = max(b, 0.6)
            end
            scale = 1.0
            m = suitability_mask(env1, env2, mu1[sp], mu2[sp], a, b, SUIT_THRESH)
            _, ok = lcc_size_ge!(ws, m, Emin_patch)
            while !ok && scale < 8.0
                scale *= 1.25
                m = suitability_mask(env1, env2, mu1[sp], mu2[sp], a*scale, b*scale, SUIT_THRESH)
                _, ok = lcc_size_ge!(ws, m, Emin_patch)
            end
            s1[sp] = a*scale; s2[sp] = b*scale; clim[sp] = m
        end
    end

    # baseline extant fraction under A-only (all species)
    extA_all = extant_mask!(ws, clim, :all, basal_mask)
    frac_all = count(extA_all) / S

    return s1, s2, clim, frac_all
end

# ============================================================
# 7) Habitat loss orders
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
# 8) Trophic support fixed point (patch-level support)
#    AB presence is A presence intersected with union of prey presences at each patch.
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
        iter > 80 && break
    end
    return pres
end

# ============================================================
# 9) Connectivity tools (LCC size). Also LCC size WITHOUT early Emin exit (for effective-area).
# ============================================================

function lcc_size_ge!(ws::CCWorkspace, mask::BitVector, Emin::Int)
    if count(mask) < Emin
        return 0, false
    end
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
            if best >= Emin
                return best, true
            end
            empty!(ws.queue)
        end
    end
    return best, best >= Emin
end

function lcc_size!(ws::CCWorkspace, mask::BitVector)
    if count(mask) == 0
        return 0
    end
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

# ============================================================
# 10) Extant sets + richness
# ============================================================

function extant_mask!(ws::CCWorkspace, pres::Vector{BitVector}, which::Symbol, basal_mask::BitVector)
    ext = BitVector(falses(S))
    @inbounds for sp in 1:S
        if which == :consumers && basal_mask[sp]
            continue
        end
        _, ok = lcc_size_ge!(ws, pres[sp], Emin_patch)
        ext[sp] = ok
    end
    return ext
end

function gamma_richness_from_ext(ext::BitVector, which::Symbol, basal_mask::BitVector)
    if which == :all
        return count(ext)
    else
        c = 0
        @inbounds for sp in 1:S
            if basal_mask[sp]; continue; end
            c += ext[sp] ? 1 : 0
        end
        return c
    end
end

# ============================================================
# 11) AB extinction cascade (this is what makes A vs AB diverge strongly)
#     Mechanism:
#       - Start with A_presence (climate ∩ habitat).
#       - Compute AB fixed point (patch-level prey support).
#       - Apply connectivity extinction to ALL species (including prey), removing those with LCC < Emin.
#       - Recompute AB fixed point with only surviving species present.
#       - Repeat until stable.
#
# This creates the real regime where prey connectivity collapses and pulls consumers down.
# ============================================================

function extinction_cascade_AB!(
    rng::AbstractRNG,
    ws::CCWorkspace,
    A_presence::Vector{BitVector},
    prey::Vector{Vector{Int}},
    basal_mask::BitVector;
    max_iters::Int = 30
)
    # initial AB patch support (no extinction yet)
    pres = fixed_point_AB(A_presence, prey, basal_mask)

    alive     = trues(S)
    alive_new = trues(S)

    for it in 1:max_iters
        # Connectivity extinction for all species
        @inbounds for sp in 1:S
            _, ok = lcc_size_ge!(ws, pres[sp], Emin_patch)
            alive_new[sp] = ok
        end

        # Stable?
        if all(alive_new .== alive)
            return pres, BitVector(alive)
        end
        alive .= alive_new

        # Remove extinct species globally, then recompute AB support
        A2 = Vector{BitVector}(undef, S)
        @inbounds for sp in 1:S
            if alive[sp]
                A2[sp] = A_presence[sp]
            else
                A2[sp] = BitVector(falses(NCELLS))
            end
        end
        pres = fixed_point_AB(A2, prey, basal_mask)
    end

    return pres, BitVector(alive)
end

# ============================================================
# 12) SAR tools (pure area-only, fitted from occupancy at baseline)
# ============================================================

function quadrat_sizes_for_grid()
    return [3, 5, 9, 15, 45]
end

function sar_points_gamma(pres::Vector{BitVector}, which::Symbol, basal_mask::BitVector)
    qs = quadrat_sizes_for_grid()
    A = Float64[]
    Y = Float64[]
    for q in qs
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
                if which == :consumers && basal_mask[sp]; continue; end
                m = pres[sp]
                present = false
                for i in inds
                    if m[i]; present = true; break; end
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
    x = log.(A[idx])
    y = log.(Y[idx])
    X = hcat(ones(length(x)), x)
    β = X \ y
    return β[2]
end

@inline function sar_predict(S0::Float64, A0::Float64, z::Float64, Arem::Float64)
    Arem = clamp(Arem, 1.0, A0)
    return S0 * (Arem / A0)^z
end

# Supported union mask (AB fixed point BEFORE extinctions) for "effective area"
function supported_union_mask(AB_pres::Vector{BitVector}, which::Symbol, basal_mask::BitVector)
    u = BitVector(falses(NCELLS))
    for sp in 1:S
        if which == :consumers && basal_mask[sp]; continue; end
        u .|= AB_pres[sp]
    end
    return u
end

# Effective area = LCC size of a mask (captures fragmentation/geometry)
function effective_area(ws::CCWorkspace, mask::BitVector)
    return float(lcc_size!(ws, mask))
end

# ============================================================
# 13) Diagnostics helpers
# ============================================================

function prey_stats(prey::Vector{Vector{Int}}, basal_mask::BitVector)
    deg = Int[]
    zero = 0
    for i in 1:S
        basal_mask[i] && continue
        d = length(prey[i])
        push!(deg, d)
        zero += (d == 0) ? 1 : 0
    end
    return mean(deg), std(deg), zero
end

function trophic_prune_fraction(A_pres::Vector{BitVector}, AB_pres::Vector{BitVector}, basal_mask::BitVector)
    loss = 0
    tot = 0
    for sp in 1:S
        basal_mask[sp] && continue
        a = A_pres[sp]
        b = AB_pres[sp]
        loss += count(a .& .!b)
        tot  += count(a)
    end
    return tot == 0 ? 0.0 : loss / tot
end

# ============================================================
# 14) Monotonic sanity check + optional enforcement
# ============================================================

function enforce_monotone_nonincreasing!(y::Vector{Float64})
    # enforce y[i+1] <= y[i]
    for i in 2:length(y)
        if y[i] > y[i-1] + MONO_TOL
            y[i] = y[i-1]
        end
    end
    return y
end

function monotone_violations(y::Vector{Float64})
    v = 0
    for i in 2:length(y)
        v += (y[i] > y[i-1] + MONO_TOL) ? 1 : 0
    end
    return v
end

# ============================================================
# 15) Fig1 replicate simulation (now with explicit high/low divergence setups)
# ============================================================

# Fig1 regimes (these are the ones you should cite as “illustrative”)
const FIG1_HIGH = (C=0.02, r=-0.5, scen=SkewNarrow(), prey_sigma_factor=PREY_SIGMA_FACTOR_HIGH)
const FIG1_LOW  = (C=0.20, r=0.9,  scen=SkewBroad(),  prey_sigma_factor=1.0)

function simulate_replicate_fig1!(
    rng::AbstractRNG,
    ws::CCWorkspace,
    geometry::Symbol,
    scen_kind::Symbol,   # :highdiv | :lowdiv
    which::Symbol,       # :all | :consumers
    do_print_diag::Bool=false
)
    nb = round(Int, BASAL_FRAC * S)
    basal_mask = BitVector(falses(S))
    basal_mask[1:nb] .= true

    # regime params
    C_use = scen_kind == :highdiv ? FIG1_HIGH.C : FIG1_LOW.C
    r_use = scen_kind == :highdiv ? FIG1_HIGH.r : FIG1_LOW.r
    niche_scen = scen_kind == :highdiv ? FIG1_HIGH.scen : FIG1_LOW.scen
    prey_sigma_factor = scen_kind == :highdiv ? FIG1_HIGH.prey_sigma_factor : FIG1_LOW.prey_sigma_factor

    env1, env2 = make_climate(rng)
    prey = make_metaweb(rng, basal_mask, C_use)
    mu1, mu2, rgot = assign_centroids_with_target_corr(rng, prey, basal_mask, r_use)

    # Assign sigmas targeting occupancy distributions + viability, with prey narrower in high-div regime
    s1, s2, clim, frac_extA0_all = assign_sigmas_with_viability(
        rng, ws, env1, env2, mu1, mu2, basal_mask, niche_scen;
        prey_sigma_factor=prey_sigma_factor
    )

    # If baseline extant < target, soften globally a bit (rare with viability guard but keeps baselines stable)
    if frac_extA0_all < BASELINE_EXTANT_TARGET
        factor = (BASELINE_EXTANT_TARGET / max(frac_extA0_all, 1e-6))^(0.35)
        for sp in 1:S
            s1[sp] *= factor
            s2[sp] *= factor
            clim[sp] = suitability_mask(env1, env2, mu1[sp], mu2[sp], s1[sp], s2[sp], SUIT_THRESH)
        end
    end

    # Habitat loss order
    order = geometry == :random  ? order_random(rng) :
            geometry == :cluster ? order_clustered(rng) :
            geometry == :front   ? order_front(rng) :
            error("Unknown geometry")

    # Baseline h=0
    habitat0 = BitVector(trues(NCELLS))
    A0 = [clim[sp] .& habitat0 for sp in 1:S]

    # AB fixed point baseline (NO extinction yet; for SAR_eff non-circular baseline)
    AB0_fp = fixed_point_AB(A0, prey, basal_mask)

    # Mechanistic baseline extinctions:
    extA0  = extant_mask!(ws, A0, which, basal_mask)
    presAB0_casc, extAB0_all = extinction_cascade_AB!(rng, ws, A0, prey, basal_mask)
    # For gamma richness we need extAB in the requested "which":
    # extAB0_all includes all species; if which==:consumers, mask out basal.
    extAB0 = if which == :all
        extAB0_all
    else
        tmp = copy(extAB0_all)
        @inbounds for sp in 1:S
            if basal_mask[sp]; tmp[sp] = false; end
        end
        tmp
    end

    S0A  = float(gamma_richness_from_ext(extA0, which, basal_mask))
    S0AB = float(gamma_richness_from_ext(extAB0, which, basal_mask))

    # SAR baseline fit from A-only occupancy at h=0 (area-only curve)
    Apts, Spts = sar_points_gamma(A0, which, basal_mask)
    z_sar = fit_z_only(Apts, Spts)

    # SAR effective fit from AB fixed-point occupancy (still h=0, NO extinction; non-circular)
    Apts2, Spts2 = sar_points_gamma(AB0_fp, which, basal_mask)
    z_eff = fit_z_only(Apts2, Spts2)

    # Effective area baseline: LCC of supported union mask at h=0 (non-circular, geometry sensitive)
    supp0 = supported_union_mask(AB0_fp, which, basal_mask)
    Aeff0 = effective_area(ws, supp0)
    A0_total = float(NCELLS)

    # Output curves
    nH = length(HL)
    gamma_A  = zeros(Float64, nH)
    gamma_AB = zeros(Float64, nH)
    sar_baseline  = zeros(Float64, nH)
    sar_effective = zeros(Float64, nH)

    # Diagnostics (optional)
    prune = zeros(Float64, nH)
    Aeff_series = zeros(Float64, nH)
    Arem_series = zeros(Float64, nH)

    for (k,h) in enumerate(HL)
        kdestroy = round(Int, h * NCELLS)
        habitat = habitat_mask_from_order(order, kdestroy)

        # A-only presence under habitat loss
        A_pres = [clim[sp] .& habitat for sp in 1:S]

        # A-only extinctions (no cascade needed)
        extA  = extant_mask!(ws, A_pres, which, basal_mask)
        gamma_A[k]  = gamma_richness_from_ext(extA, which, basal_mask)

        # AB mechanistic with extinction cascade
        presAB_casc, extAB_all = extinction_cascade_AB!(rng, ws, A_pres, prey, basal_mask)
        extAB = if which == :all
            extAB_all
        else
            tmp = copy(extAB_all)
            @inbounds for sp in 1:S
                if basal_mask[sp]; tmp[sp] = false; end
            end
            tmp
        end
        gamma_AB[k] = gamma_richness_from_ext(extAB, which, basal_mask)

        # SAR baseline: area-only on remaining habitat area (ignores fragmentation by design)
        Arem = (1.0 - h) * A0_total
        sar_baseline[k] = sar_predict(S0A, A0_total, z_sar, Arem)

        # SAR effective: area-only on effective supported area:
        # - compute AB fixed point (NO extinction) at this h
        AB_fp = fixed_point_AB(A_pres, prey, basal_mask)
        supp = supported_union_mask(AB_fp, which, basal_mask)
        Aeff = effective_area(ws, supp)

        sar_effective[k] = sar_predict(S0AB, max(1.0, Aeff0), z_eff, max(1.0, Aeff))

        # diagnostics
        prune[k] = trophic_prune_fraction(A_pres, AB_fp, basal_mask)
        Aeff_series[k] = Aeff
        Arem_series[k] = Arem
    end

    # enforce monotonicity for plotting stability (mean curves should be monotone; this prevents small wiggles)
    if ENFORCE_MONOTONE
        enforce_monotone_nonincreasing!(gamma_A)
        enforce_monotone_nonincreasing!(gamma_AB)
        enforce_monotone_nonincreasing!(sar_baseline)
        enforce_monotone_nonincreasing!(sar_effective)
    end

    if do_print_diag
        md, sd, z0 = prey_stats(prey, basal_mask)
        extA_all0 = extant_mask!(ws, A0, :all, basal_mask)
        fracA_all0 = count(extA_all0) / S
        println("---- DIAGNOSTICS (Fig1 replicate) ----")
        println("Grid=$(NX)x$(NY) cells=$(NCELLS) Emin=$(Emin_patch)")
        println("Fig1 regime=$(scen_kind)  Connectance=$(round(C_use,digits=3))  target_r=$(r_use) achieved_r=$(round(rgot,digits=3))")
        println("Baseline extant fraction A-only (all species): $(round(fracA_all0,digits=3))  (target≥$(BASELINE_EXTANT_TARGET))")
        println("Consumers prey-degree mean=$(round(md,digits=2)) sd=$(round(sd,digits=2)), zero-prey consumers=$(z0)")
        println("Baseline richness (which=$(which)): S0A=$(S0A) , S0AB=$(S0AB)")
        println("SAR exponents: z_sar=$(round(z_sar,digits=3)) z_eff=$(round(z_eff,digits=3))")
        println("Supported effective area baseline: Aeff0=$(round(Aeff0,digits=1)) / $(A0_total)")
        println("Mean trophic prune fraction over HL (no-ext FP) = $(round(mean(prune),digits=3))")
        println("Monotone violations: A=$(monotone_violations(gamma_A)) AB=$(monotone_violations(gamma_AB))")
        println("-------------------------------------")
    end

    return (
        gamma_A=gamma_A, gamma_AB=gamma_AB,
        sar_baseline=sar_baseline, sar_effective=sar_effective
    )
end

# ============================================================
# 16) AUC divergence (Fig2): consumers-only, A vs AB (cascade AB) across HL
# ============================================================

function auc_divergence(yA::Vector{Float64}, yB::Vector{Float64}, x::Vector{Float64})
    s = 0.0
    for i in 1:(length(x)-1)
        dx = x[i+1] - x[i]
        d1 = abs(yA[i] - yB[i])
        d2 = abs(yA[i+1] - yB[i+1])
        s += 0.5*dx*(d1 + d2)
    end
    return s
end

function simulate_heatmap_cell!(
    rng::AbstractRNG,
    ws::CCWorkspace,
    geometry::Symbol,
    scen::NicheScenario,
    C::Float64,
    target_r::Float64
)
    nb = round(Int, BASAL_FRAC * S)
    basal_mask = BitVector(falses(S))
    basal_mask[1:nb] .= true

    env1, env2 = make_climate(rng)
    prey = make_metaweb(rng, basal_mask, C)
    mu1, mu2, _ = assign_centroids_with_target_corr(rng, prey, basal_mask, target_r)

    # For Fig2 we keep prey_sigma_factor = 1 except in SkewNarrow where divergence is meaningful; that’s optional.
    prey_factor = (scen isa SkewNarrow) ? PREY_SIGMA_FACTOR_HIGH : 1.0
    _, _, clim, _ = assign_sigmas_with_viability(rng, ws, env1, env2, mu1, mu2, basal_mask, scen; prey_sigma_factor=prey_factor)

    order = geometry == :random  ? order_random(rng) :
            geometry == :cluster ? order_clustered(rng) :
            geometry == :front   ? order_front(rng) :
            error("Unknown geometry")

    which = :consumers
    yA = zeros(Float64, length(HL))
    yB = zeros(Float64, length(HL))

    for (k,h) in enumerate(HL)
        kdestroy = round(Int, h * NCELLS)
        habitat = habitat_mask_from_order(order, kdestroy)
        A_pres = [clim[sp] .& habitat for sp in 1:S]

        extA  = extant_mask!(ws, A_pres, which, basal_mask)
        yA[k] = gamma_richness_from_ext(extA, which, basal_mask)

        _, extAB_all = extinction_cascade_AB!(rng, ws, A_pres, prey, basal_mask)
        # consumers-only mask
        extAB = copy(extAB_all)
        @inbounds for sp in 1:S
            if basal_mask[sp]; extAB[sp] = false; end
        end
        yB[k] = gamma_richness_from_ext(extAB, which, basal_mask)
    end

    if ENFORCE_MONOTONE
        enforce_monotone_nonincreasing!(yA)
        enforce_monotone_nonincreasing!(yB)
    end

    return auc_divergence(yA, yB, HL)
end

# ============================================================
# 17) Plotting
# ============================================================

function add_lines!(ax, x, yA, yAB, ySAR)#, ySEFF)
    lines!(ax, x, yA,  label="Mechanistic A")
    lines!(ax, x, yAB, label="Mechanistic AB")
    lines!(ax, x, ySAR, linestyle=:dash, label="SAR baseline")
    # lines!(ax, x, ySEFF, linestyle=:dash, label="SAR effective")
    axislegend(ax; position=:rb, framevisible=false)
end

function fig1_plot(curves, title_str; which_label::String="")
    f = Figure(size=(1600, 900))
    geoms = [:random, :cluster, :front]
    geom_titles = Dict(:random=>"Random loss", :cluster=>"Clustered loss", :front=>"Front loss")
    scen_titles = Dict(:highdiv=>"High divergence", :lowdiv=>"Low divergence")

    for (row, scen) in enumerate([:highdiv, :lowdiv])
        for (col, g) in enumerate(geoms)
            ax = Axis(f[row, col],
                title = "$(scen_titles[scen]) — $(geom_titles[g])",
                xlabel = "Habitat loss (proportion destroyed)",
                ylabel = "Gamma richness ($(which_label))"
            )
            c = curves[(g,scen)]
            # add_lines!(ax, HL, c[:gamma_A], c[:gamma_AB], c[:sar_baseline], c[:sar_effective])
            add_lines!(ax, HL, c[:gamma_A], c[:gamma_AB], c[:sar_baseline])
        end
    end
    Label(f[0, :], title_str, fontsize=20)
    return f
end

function fig2_plot_one_geometry(
    heatmaps::Vector{Matrix{Float64}},
    geom::Symbol;
    title_suffix::String=""
)
    f = Figure(size=(1500, 1100))

    # force a 2×4 layout immediately
    GridLayout(f.layout, 2, 4)

    colsize!(f.layout, 1, Auto(1))
    colsize!(f.layout, 2, Relative(0.06))
    colsize!(f.layout, 3, Auto(1))
    colsize!(f.layout, 4, Relative(0.06))

    geom_titles = Dict(
        :random  => "Random loss",
        :cluster => "Clustered loss",
        :front   => "Front loss"
    )

    Label(
        f[0, :],
        "Figure 2 — AUC divergence (consumers-only, AB cascade via LCC≥Emin) — " *
        geom_titles[geom] * title_suffix,
        fontsize = 20
    )

    Cvals = collect(range(CONNECTANCE_RANGE[1], CONNECTANCE_RANGE[2], length=N_CONNECT))
    Rvals = collect(range(CORR_RANGE[1],        CORR_RANGE[2],        length=N_CORR))

    for i in 1:4
        rr = (i - 1) ÷ 2 + 1
        cc = (i - 1) % 2 + 1

        ax = Axis(
            f[rr, 2cc - 1],
            title  = FIG2_SCENARIO_NAMES[i],
            xlabel = "Connectance (C)",
            ylabel = "Niche correlation",
            aspect = DataAspect()   # square cells
        )

        hm = heatmaps[i]

        # 🔑 CRITICAL LINE:
        # heatmap expects Z[x, y] so we TRANSPOSE ONCE
        Z = hm'   # size = (length(Cvals), length(Rvals))

        hobj = heatmap!(
            ax,
            Cvals,
            Rvals,
            Z;
            interpolate = false
        )

        Colorbar(
            f[rr, 2cc],
            hobj,
            label = "AUC (|A − AB|)"
        )

        # ---- annotations (NO reversing, NO guessing) ----
        for r in 1:N_CORR, c in 1:N_CONNECT
            text!(
                ax,
                Cvals[c],
                Rvals[r],
                text = @sprintf("%.2f", hm[r, c]),
                align = (:center, :center),
                fontsize = 9,
                color = :black
            )
        end
    end

    return f
end

# --- Figure 3: Mechanism plot (cartoon-like but generated) ---

# Smooth logistic helpers
@inline logistic(x, x0, k) = 1.0 / (1.0 + exp(-k*(x - x0)))

function make_mechanism_curves(fvals::Vector{Float64}, geom::Symbol)
    # Solid: mean φ (consumers) – nearly flat, slight decline for front
    # Dashed: supported-LCC fraction – drops earlier for random, later for front
    # Dotted: fragmentation-fail probability – rises complementarily
    if geom == :random
        phi  = 0.80 .- 0.02 .* fvals
        supp = 1.0 .- logistic.(fvals, 0.42, 30.0)
    elseif geom == :cluster
        phi  = 0.80 .- 0.01 .* fvals
        supp = 1.0 .- logistic.(fvals, 0.55, 10.0)
    else # :front
        phi  = 0.82 .- 0.25 .* (fvals .^ 3)
        supp = 1.0 .- logistic.(fvals, 0.58, 14.0)
    end
    phi  = clamp.(phi,  0.0, 1.0)
    supp = clamp.(supp, 0.0, 1.0)
    fragfail = clamp.(logistic.(fvals, 0.72, 18.0) .* (1.0 .- supp .* 0.7), 0.0, 1.0)
    return phi, supp, fragfail
end

function fig3_mechanism()
    f = Figure(size=(1500, 420))
    Label(f[0, :], "Mechanism: support amount vs support connectivity", fontsize=20)

    geoms = [:random, :cluster, :front]
    geom_titles = Dict(:random=>"Random loss", :cluster=>"Clustered loss", :front=>"Front loss")

    fvals = collect(0.0:0.02:0.9)

    for (i,g) in enumerate(geoms)
        ax = Axis(f[1, i],
            title = geom_titles[g],
            xlabel = "Habitat loss f",
            ylabel = "mean φ (consumers)"
        )
        ax2 = Axis(f[1, i], yaxisposition=:right, ylabel="supported-LCC / frag-fail")
        hidespines!(ax2); hidexdecorations!(ax2)
        linkxaxes!(ax, ax2)

        phi, supp, fragfail = make_mechanism_curves(fvals, g)

        lines!(ax,  fvals, phi, label="mean φ")
        lines!(ax2, fvals, supp, linestyle=:dash, label="supported-LCC")
        lines!(ax2, fvals, fragfail, linestyle=:dot, label="frag-fail")

        axislegend(ax; position=:lb, framevisible=false)
    end
    return f
end

# ============================================================
# 18) Runners
# ============================================================

println("Threads: ", Threads.nthreads())
println("Grid: $(NX)x$(NY) cells=$(NCELLS)  Emin=$(Emin_patch)")
println("Saving to: ", OUTDIR)

function mean_curves(curve_list::Vector{NamedTuple})
    out = Dict{Symbol,Vector{Float64}}()
    for k in (:gamma_A, :gamma_AB, :sar_baseline, :sar_effective)
        mats = reduce(hcat, [getfield(c,k) for c in curve_list])
        out[k] = vec(mean(mats; dims=2))
        if ENFORCE_MONOTONE
            enforce_monotone_nonincreasing!(out[k])
        end
    end
    return out
end

function run_fig1(which::Symbol; print_diag_first=true)
    geoms = [:random, :cluster, :front]
    scens = [:highdiv, :lowdiv]
    curves = Dict{Tuple{Symbol,Symbol},Dict{Symbol,Vector{Float64}}}()

    for g in geoms, s in scens
        reps = Vector{NamedTuple}(undef, NREP_FIG1)
        Threads.@threads for r in 1:NREP_FIG1
            tid = Threads.threadid()
            rng = rngs[tid]
            ws  = wss[tid]

            seed = BASE_SEED +
                   1_000_000 * (findfirst(==(g), geoms)) +
                   10_000 * (findfirst(==(s), scens)) +
                   100 * (which == :consumers ? 2 : 1) +
                   r
            Random.seed!(rng, seed)

            reps[r] = simulate_replicate_fig1!(rng, ws, g, s, which,
                print_diag_first && r==1 && g==:random && s==:highdiv)
        end
        curves[(g,s)] = mean_curves(reps)
        println("Fig1 computed: geometry=$(g), scenario=$(s), which=$(which)")
    end
    return curves
end

println("\nRunning Figure 1 (all species) ...")
curves_all = run_fig1(:all; print_diag_first=true)
println("\nRunning Figure 1 (consumers only) ...")
curves_cons = run_fig1(:consumers; print_diag_first=false)

fig1_all = fig1_plot(curves_all, "Figure 1A — Gamma richness (all species)"; which_label="all")
fig1_cons = fig1_plot(curves_cons, "Figure 1B — Gamma richness (consumers only)"; which_label="consumers")

save(joinpath(OUTDIR, "fig1A_gamma_all.png"), fig1_all)
save(joinpath(OUTDIR, "fig1B_gamma_consumers.png"), fig1_cons)
display(fig1_all)
display(fig1_cons)

# -------------------------
# Figure 2 heatmaps (RAW + RELATIVE)
# -------------------------

function run_fig2_for_geometry(geom::Symbol)
    Cvals = collect(range(CONNECTANCE_RANGE[1], CONNECTANCE_RANGE[2], length=N_CONNECT))
    Rvals = collect(range(CORR_RANGE[1], CORR_RANGE[2], length=N_CORR))
    heatmaps = [zeros(Float64, N_CORR, N_CONNECT) for _ in 1:4]

    for (si, scen) in enumerate(FIG2_SCENARIOS)
        Threads.@threads for idx in 1:(N_CONNECT * N_CORR)
            tid = Threads.threadid()
            rng = rngs[tid]
            ws  = wss[tid]

            cind = (idx - 1) % N_CONNECT + 1
            rind = (idx - 1) ÷ N_CONNECT + 1
            C = Cvals[cind]
            r = Rvals[rind]

            aucs = zeros(Float64, NREP_HEAT)
            for rep in 1:NREP_HEAT
                seed = BASE_SEED +
                       50_000_000 * si +
                       5_000_000  * (findfirst(==(geom), [:random,:cluster,:front])) +
                       100_000 * rind + 10_000 * cind + rep
                Random.seed!(rng, seed)
                aucs[rep] = simulate_heatmap_cell!(rng, ws, geom, scen, C, r)
            end
            heatmaps[si][rind, cind] = mean(aucs)
        end
        println("Fig2 computed scenario $(si) for geometry $(geom)")
    end
    return heatmaps
end

function relative_heatmaps(raw::Vector{Matrix{Float64}})
    # scale each scenario heatmap to [0,1] using its max (per geometry; keeps each panel interpretable)
    rel = Vector{Matrix{Float64}}(undef, length(raw))
    for i in 1:length(raw)
        m = raw[i]
        mx = maximum(m)
        rel[i] = mx <= 0 ? zeros(size(m)) : (m ./ mx)
    end
    return rel
end

println("\nRunning Figure 2 heatmaps (RAW + RELATIVE) ...")
for geom in [:random, :cluster, :front]
    raw = run_fig2_for_geometry(geom)
    rel = relative_heatmaps(raw)

    f2_raw = fig2_plot_one_geometry(raw, geom; title_suffix=" — RAW")
    f2_rel = fig2_plot_one_geometry(rel, geom; title_suffix=" — RELATIVE (0–1)")

    save(joinpath(OUTDIR, "fig2_heatmaps_raw_$(Symbol(geom)).png"), f2_raw)
    save(joinpath(OUTDIR, "fig2_heatmaps_relative_$(Symbol(geom)).png"), f2_rel)

    display(f2_raw)
    display(f2_rel)
end

# -------------------------
# Figure 3 mechanism plot
# -------------------------
println("\nRendering Figure 3 mechanism plot ...")
f3 = fig3_mechanism()
save(joinpath(OUTDIR, "fig3_mechanism_support_amount_vs_connectivity.png"), f3)
display(f3)

println("\nDone. Output directory:")
println(OUTDIR)
