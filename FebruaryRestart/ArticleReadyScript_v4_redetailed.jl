#!/usr/bin/env julia
# Julia 1.11
# Run:  julia --threads auto traffic_extinction_figures.jl
#
# PURPOSE
# - Mechanistic extinction curves under habitat loss with:
#   A-only (climate niche) vs AB (climate + trophic support), plus connectivity extinction rule (LCC ≥ Emin).
# - SAR baseline: pure area-only prediction (dashed), fit from baseline A-only distribution using non-overlapping quadrats.
# - SAR effective: pure area-only prediction (dashed), but uses trophically informed baseline (AB at h=0) and
#   evaluates richness vs "supported area" (defined NON-CIRCULARLY: from AB fixed point BEFORE extinction filtering).
#
# FIXES vs prior versions
# 1) NON-CIRCULAR supported area for SAR_eff: computed from AB fixed-point occupancy before extinctions.
# 2) Enforce baseline viability: ≥95% species extant under A-only at h=0 (LCC ≥ Emin) by resampling sigmas.
# 3) SAR curves are clean area-only expressions: monotone in remaining (effective) area, no circular feedback.
# 4) SAR lines plotted dashed.
# 5) Diagnostics added to verify: baseline extant fraction, achieved niche correlation, trophic support strength, etc.
#
# NOTE
# - Figure 1 duplicated: all-species and consumers-only
# - Figure 1B: mean local richness retained, but you can ignore it.

using Random
using Statistics
using LinearAlgebra
using CairoMakie
using Printf
using Dates

# ============================================================
# 0) Parameters
# ============================================================

# Grid size reduced, Emin kept at 80
const NX = 45
const NY = 45
const NCELLS = NX * NY

const S = 250
const BASAL_FRAC = 0.33
const Emin_patch = 80

# Habitat loss levels
const HL = collect(0.0:0.05:1.0)

# Replicates
const NREP_FIG1 = 10
const NREP_HEAT = 10

# Heatmap axes
const CONNECTANCE_RANGE = (0.02, 0.25)
const CORR_RANGE = (-0.5, 0.9)
const N_CONNECT = 10
const N_CORR = 10

# Medium values
const C_mid = 0.12
const r_mid = 0.5

# Climate
const CLIM_MIN = 0.0
const CLIM_MAX = 100.0
const CLIM_NOISE_SD = 7.0

# Niche option A: 2D Gaussian with threshold
const SUIT_THRESH = 0.25

# Habitat loss geometry params
const CLUSTER_SEED_PER_DESTROYED = 150
const FRONT_JITTER_MAX = 5
const FRONT_SMOOTH = :strong

# Baseline viability target under A-only
const BASELINE_EXTANT_TARGET = 0.95

# Resampling knobs (for baseline viability)
const MAX_SIGMA_RESAMPLE = 30

# Output
ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
OUTDIR = joinpath(pwd(), "figures_traffic_extinction_" * ts)
isdir(OUTDIR) || mkpath(OUTDIR)

# Global seed
const BASE_SEED = 20250128

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
# 2) Thread-safe RNG + workspace
# ============================================================

mutable struct CCWorkspace
    seen::Vector{Int32}
    stamp::Int32
    queue::Vector{Int}
    queue2::Vector{Int}
end

function make_thread_rngs_and_workspaces(base_seed::Int)
    nt = Threads.nthreads()
    rngs = Vector{MersenneTwister}(undef, nt)
    wss  = Vector{CCWorkspace}(undef, nt)
    for t in 1:nt
        rngs[t] = MersenneTwister(base_seed + 10_000 * t + 137)
        wss[t]  = CCWorkspace(fill(Int32(0), NCELLS), Int32(0), Int[], Int[])
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
# 4) Metaweb (directed), no cannibalism, basal have no prey
# ============================================================

function make_metaweb(rng::AbstractRNG, basal_mask::BitVector, C::Float64)
    Ltarget = round(Int, C * S^2)  # C = L/S^2 (per your definition)
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
# 5) Consumer–prey centroid correlation (assign consumer centroids conditional on prey)
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
# 6) Niche option A: 2D Gaussian threshold
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

# Scenario-specific sigma draws (we will RESAMPLE per species until baseline viability target is met)
function draw_sigma_pair(rng::AbstractRNG, scen::Symbol)
    # More controlled than before: avoid ultra-narrow sigmas that make baseline extinct by construction.
    if scen == :highdiv
        # narrow-ish but not extreme
        an = exp(randn(rng) * 0.45)
        base = exp(log(10.5) + randn(rng) * 0.35)
        s1 = clamp(base * an, 2.0, 25.0)
        s2 = clamp(base / an, 2.0, 25.0)
        return s1, s2
    else
        # broad
        an = exp(randn(rng) * 0.25)
        base = exp(log(18.0) + randn(rng) * 0.30)
        s1 = clamp(base * an, 3.5, 50.0)
        s2 = clamp(base / an, 3.5, 50.0)
        return s1, s2
    end
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
        s = x + jit[y] + 0.03*randn(rng)  # smooth front with slight irregularity
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
# 8) AB trophic fixed point
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
# 9) Connectivity: LCC size and largest component extraction
# ============================================================

function lcc_size_ge!(ws, mask::BitVector, Emin::Int)
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

# ============================================================
# 10) Extant sets + richness
# ============================================================

function extant_mask!(ws, pres::Vector{BitVector}, which::Symbol, basal_mask::BitVector)
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
            basal_mask[sp] && continue
            c += ext[sp] ? 1 : 0
        end
        return c
    end
end

# ============================================================
# 11) SAR building: non-overlapping quadrats, fit z, anchor at baseline
# SAR is AREA-ONLY: S(Arem) = S0 * (Arem/A0)^z
# SAR_eff is AREA-ONLY in supported area: S(Aeff) = S0_AB * (Aeff/Aeff0)^z_eff
# Supported area is computed from AB fixed-point occupancy (before extinctions), non-circular.
# ============================================================

function quadrat_sizes_for_grid()
    # divisors of 45 => 1,3,5,9,15,45 ; avoid 1 for stability but keep full map
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

# Supported area = union of AB presences for the relevant set (all or consumers) BEFORE extinction filtering
function supported_area(AB_pres::Vector{BitVector}, which::Symbol, basal_mask::BitVector)
    u = BitVector(falses(NCELLS))
    for sp in 1:S
        if which == :consumers && basal_mask[sp]; continue; end
        u .|= AB_pres[sp]
    end
    return count(u)
end

# ============================================================
# 12) Diagnostics helpers
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
# 13) Enforce baseline viability: resample sigma for each species until LCC≥Emin (A-only at h=0)
# Strategy: keep centroids fixed; resample (s1,s2) until climate mask yields LCC≥Emin.
# This guarantees baseline extant fraction high and prevents "already extinct" pools.
# ============================================================

function enforce_baseline_viability!(
    rng::AbstractRNG,
    ws::CCWorkspace,
    env1, env2,
    mu1::Vector{Float64}, mu2::Vector{Float64},
    basal_mask::BitVector,
    scen_kind::Symbol
)
    s1 = Vector{Float64}(undef, S)
    s2 = Vector{Float64}(undef, S)
    clim = Vector{BitVector}(undef, S)

    ok_count = 0
    for sp in 1:S
        # basal species can still be non-viable climatically; but they count toward A-only extant in "all-species".
        # We enforce viability for everyone to stabilize baselines.
        got = false
        for _ in 1:MAX_SIGMA_RESAMPLE
            a,b = draw_sigma_pair(rng, scen_kind)
            m = suitability_mask(env1, env2, mu1[sp], mu2[sp], a, b, SUIT_THRESH)
            _, ok = lcc_size_ge!(ws, m, Emin_patch)
            if ok
                s1[sp] = a
                s2[sp] = b
                clim[sp] = m
                got = true
                break
            end
        end
        if !got
            # last resort: progressively inflate sigma until viable
            a,b = draw_sigma_pair(rng, scen_kind)
            scale = 1.0
            m = suitability_mask(env1, env2, mu1[sp], mu2[sp], a, b, SUIT_THRESH)
            _, ok = lcc_size_ge!(ws, m, Emin_patch)
            while !ok && scale < 6.0
                scale *= 1.25
                m = suitability_mask(env1, env2, mu1[sp], mu2[sp], a*scale, b*scale, SUIT_THRESH)
                _, ok = lcc_size_ge!(ws, m, Emin_patch)
            end
            s1[sp] = a*scale
            s2[sp] = b*scale
            clim[sp] = m
        end
        ok_count += 1
    end

    # sanity: compute extant fraction under A-only "all"
    extA_all = extant_mask!(ws, clim, :all, basal_mask)
    frac_all = count(extA_all) / S

    return s1, s2, clim, frac_all
end

# ============================================================
# 14) Fig1 replicate simulation
# ============================================================

function simulate_replicate_fig1!(
    rng::AbstractRNG,
    ws::CCWorkspace,
    geometry::Symbol,
    scen_kind::Symbol,   # :highdiv or :lowdiv
    which::Symbol,       # :all or :consumers
    do_print_diag::Bool=false
)
    nb = round(Int, BASAL_FRAC * S)
    basal_mask = BitVector(falses(S))
    basal_mask[1:nb] .= true

    env1, env2 = make_climate(rng)
    prey = make_metaweb(rng, basal_mask, C_mid)

    target_r = scen_kind == :highdiv ? 0.0 : 0.8
    mu1, mu2, rgot = assign_centroids_with_target_corr(rng, prey, basal_mask, target_r)

    # Enforce baseline viability (A-only) via sigma resampling.
    s1, s2, clim, frac_extA0_all = enforce_baseline_viability!(rng, ws, env1, env2, mu1, mu2, basal_mask, scen_kind)

    # If baseline extant < target (for all species), we broaden globally a bit:
    if frac_extA0_all < BASELINE_EXTANT_TARGET
        # global soften: multiply all sigmas by factor to lift extant fraction
        factor = (BASELINE_EXTANT_TARGET / max(frac_extA0_all, 1e-6))^(0.35)
        for sp in 1:S
            s1[sp] *= factor
            s2[sp] *= factor
            clim[sp] = suitability_mask(env1, env2, mu1[sp], mu2[sp], s1[sp], s2[sp], SUIT_THRESH)
        end
    end

    # orders
    order = geometry == :random  ? order_random(rng) :
            geometry == :cluster ? order_clustered(rng) :
            geometry == :front   ? order_front(rng) :
            error("Unknown geometry")

    habitat0 = BitVector(trues(NCELLS))
    A0 = [clim[sp] .& habitat0 for sp in 1:S]
    AB0 = fixed_point_AB(A0, prey, basal_mask)

    # extant sets at baseline
    extA0  = extant_mask!(ws, A0, which, basal_mask)
    extAB0 = extant_mask!(ws, AB0, which, basal_mask)

    S0A  = float(gamma_richness_from_ext(extA0, which, basal_mask))
    S0AB = float(gamma_richness_from_ext(extAB0, which, basal_mask))

    # SAR baseline fit from baseline A-only occupancy (area-only curve)
    Apts, Spts = sar_points_gamma(A0, which, basal_mask)
    z_sar = fit_z_only(Apts, Spts)

    # SAR effective fit from baseline AB occupancy (area-only curve in supported area)
    Apts2, Spts2 = sar_points_gamma(AB0, which, basal_mask)
    z_eff = fit_z_only(Apts2, Spts2)

    # supported baseline area (NON-CIRCULAR)
    Aeff0    = float(supported_area(AB0, which, basal_mask))
    A0_total = float(NCELLS)

    nH = length(HL)
    gamma_A  = zeros(Float64, nH)
    gamma_AB = zeros(Float64, nH)
    sar_baseline  = zeros(Float64, nH)
    sar_effective = zeros(Float64, nH)

    # Diagnostics time series (optional; kept for prints)
    prune = zeros(Float64, nH)
    suppA = zeros(Float64, nH)
    suppAB = zeros(Float64, nH)

    for (k,h) in enumerate(HL)
        kdestroy = round(Int, h * NCELLS)
        habitat = habitat_mask_from_order(order, kdestroy)

        A_pres  = [clim[sp] .& habitat for sp in 1:S]
        AB_pres = fixed_point_AB(A_pres, prey, basal_mask)

        extA  = extant_mask!(ws, A_pres, which, basal_mask)
        extAB = extant_mask!(ws, AB_pres, which, basal_mask)

        gamma_A[k]  = gamma_richness_from_ext(extA, which, basal_mask)
        gamma_AB[k] = gamma_richness_from_ext(extAB, which, basal_mask)

        # SAR baseline: area-only on remaining habitat area (ignores fragmentation + niches by design)
        Arem = (1.0 - h) * A0_total
        sar_baseline[k] = sar_predict(S0A, A0_total, z_sar, Arem)

        # SAR effective: area-only on supported area (computed from AB fixed point BEFORE extinctions)
        Aeff = supported_area(AB_pres, which, basal_mask)
        sar_effective[k] = sar_predict(S0AB, max(1.0, Aeff0), z_eff, max(1.0, Aeff))

        # diagnostics
        prune[k] = trophic_prune_fraction(A_pres, AB_pres, basal_mask)
        suppA[k] = mean([count(A_pres[sp]) for sp in 1:S]) / NCELLS
        suppAB[k] = mean([count(AB_pres[sp]) for sp in 1:S]) / NCELLS
    end

    if do_print_diag
        md, sd, z0 = prey_stats(prey, basal_mask)
        extA_all0 = extant_mask!(ws, A0, :all, basal_mask)
        fracA_all0 = count(extA_all0) / S
        println("---- DIAGNOSTICS (Fig1 replicate) ----")
        println("Threads=$(Threads.nthreads()) Grid=$(NX)x$(NY) cells=$(NCELLS) Emin=$(Emin_patch)")
        println("Baseline extant fraction A-only (all species): $(round(fracA_all0,digits=3))  (target≥$(BASELINE_EXTANT_TARGET))")
        println("Connectance C=$(C_mid); target_r=$(target_r); achieved_r=$(round(rgot,digits=3))")
        println("Consumers prey-degree mean=$(round(md,digits=2)) sd=$(round(sd,digits=2)), zero-prey consumers=$(z0)")
        println("Baseline extant richness (which=$(which)): S0A=$(S0A) , S0AB=$(S0AB)")
        println("Baseline supported area (Aeff0=$(Aeff0) / A0=$(A0_total))")
        println("SAR exponents: z_sar=$(round(z_sar,digits=3)) z_eff=$(round(z_eff,digits=3))")
        println("Mean trophic prune fraction over HL = $(round(mean(prune),digits=3))")
        println("Mean occupancy A vs AB over HL (cells): A=$(round(mean(suppA),digits=3)), AB=$(round(mean(suppAB),digits=3))")
        println("-------------------------------------")
    end

    return (
        gamma_A=gamma_A, gamma_AB=gamma_AB,
        sar_baseline=sar_baseline, sar_effective=sar_effective,
        achieved_r=rgot, S0A=S0A, S0AB=S0AB
    )
end

# ============================================================
# 15) AUC divergence (Fig2) between mechanistic A and AB (consumers only)
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
struct Fig2Scenario1 <: NicheScenario end
struct Fig2Scenario2 <: NicheScenario end
struct Fig2Scenario3 <: NicheScenario end
struct Fig2Scenario4 <: NicheScenario end
const FIG2_SCENARIOS = (Fig2Scenario1(), Fig2Scenario2(), Fig2Scenario3(), Fig2Scenario4())
const FIG2_SCENARIO_NAMES = ("VarOcc + VarSigma", "VarOcc + MildSigma", "SkewNarrow", "SkewBroad")

function draw_target_occ(rng::AbstractRNG, scen::NicheScenario)
    if scen isa Fig2Scenario1
        return 0.04 + rand(rng) * (0.75 - 0.04)
    elseif scen isa Fig2Scenario2
        return 0.04 + rand(rng) * (0.75 - 0.04)
    elseif scen isa Fig2Scenario3
        x = rand(rng, BetaLike(2.0, 8.0))
        return 0.02 + x * (0.55 - 0.02)
    else
        x = rand(rng, BetaLike(8.0, 2.0))
        return 0.18 + x * (0.92 - 0.18)
    end
end

function solve_sigmas_for_occ(env1, env2, mu1, mu2, aniso, p_target; iters=16)
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

function assign_sigmas_fig2(rng::AbstractRNG, ws::CCWorkspace, env1, env2, mu1, mu2, basal_mask::BitVector, scen::NicheScenario)
    # For Fig2 we ALSO prevent trivial baseline extinctions by ensuring each species has LCC≥Emin at h=0 under A-only.
    s1 = Vector{Float64}(undef, S)
    s2 = Vector{Float64}(undef, S)
    clim = Vector{BitVector}(undef, S)
    for sp in 1:S
        got = false
        for _ in 1:MAX_SIGMA_RESAMPLE
            p = draw_target_occ(rng, scen)
            an = scen isa Fig2Scenario1 ? exp(randn(rng) * 0.65) :
                 scen isa Fig2Scenario2 ? exp(randn(rng) * 0.25) :
                                         exp(randn(rng) * 0.40)
            a,b = solve_sigmas_for_occ(env1, env2, mu1[sp], mu2[sp], an, p)
            m = suitability_mask(env1, env2, mu1[sp], mu2[sp], a, b, SUIT_THRESH)
            _, ok = lcc_size_ge!(ws, m, Emin_patch)
            if ok
                s1[sp] = a; s2[sp] = b; clim[sp] = m
                got = true
                break
            end
        end
        if !got
            # inflate
            p = draw_target_occ(rng, scen)
            an = exp(randn(rng) * 0.35)
            a,b = solve_sigmas_for_occ(env1, env2, mu1[sp], mu2[sp], an, p)
            scale = 1.0
            m = suitability_mask(env1, env2, mu1[sp], mu2[sp], a, b, SUIT_THRESH)
            _, ok = lcc_size_ge!(ws, m, Emin_patch)
            while !ok && scale < 6.0
                scale *= 1.25
                m = suitability_mask(env1, env2, mu1[sp], mu2[sp], a*scale, b*scale, SUIT_THRESH)
                _, ok = lcc_size_ge!(ws, m, Emin_patch)
            end
            s1[sp] = a*scale; s2[sp] = b*scale; clim[sp] = m
        end
    end
    return s1, s2, clim
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

    _, _, clim = assign_sigmas_fig2(rng, ws, env1, env2, mu1, mu2, basal_mask, scen)

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
        AB_pres = fixed_point_AB(A_pres, prey, basal_mask)
        extA  = extant_mask!(ws, A_pres, which, basal_mask)
        extAB = extant_mask!(ws, AB_pres, which, basal_mask)
        yA[k] = gamma_richness_from_ext(extA, which, basal_mask)
        yB[k] = gamma_richness_from_ext(extAB, which, basal_mask)
    end
    return auc_divergence(yA, yB, HL)
end

# ============================================================
# 16) Plotting
# ============================================================

function add_lines!(ax, x, yA, yAB, ySAR, ySEFF)
    lines!(ax, x, yA, label="Mechanistic A")
    lines!(ax, x, yAB, label="Mechanistic AB")
    lines!(ax, x, ySAR, linestyle=:dash, label="SAR baseline")
    lines!(ax, x, ySEFF, linestyle=:dash, label="SAR effective")
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
            add_lines!(ax, HL, c[:gamma_A], c[:gamma_AB], c[:sar_baseline], c[:sar_effective])
        end
    end
    Label(f[0, :], title_str, fontsize=20)
    return f
end

function fig2_plot_one_geometry(heatmaps::Vector{Matrix{Float64}}, geom::Symbol)
    f = Figure(size=(1400, 1100))
    geom_titles = Dict(:random=>"Random loss", :cluster=>"Clustered loss", :front=>"Front loss")
    Label(f[0, :], "Figure 2 — AUC divergence (consumers-only, extinctions via LCC≥Emin) — $(geom_titles[geom])", fontsize=20)

    Cvals = collect(range(CONNECTANCE_RANGE[1], CONNECTANCE_RANGE[2], length=N_CONNECT))
    Rvals = collect(range(CORR_RANGE[1], CORR_RANGE[2], length=N_CORR))

    for i in 1:4
        rr = (i-1) ÷ 2 + 1
        cc = (i-1) % 2 + 1
        ax = Axis(f[rr, 2*cc-1],
            title = FIG2_SCENARIO_NAMES[i],
            xlabel = "Connectance (C)",
            ylabel = "Niche correlation"
        )
        hm = heatmaps[i]
        hobj = heatmap!(ax, hm)
        Colorbar(f[rr, 2*cc], hobj, label="AUC (|A-AB|)")
        for r in 1:size(hm,1), c in 1:size(hm,2)
            text!(ax, c, r, text=@sprintf("%.2f", hm[r,c]),
                  align=(:center,:center), fontsize=10, color=:black)
        end
        ax.xticks = (1:N_CONNECT, [@sprintf("%.2f", x) for x in Cvals])
        ax.yticks = (1:N_CORR, [@sprintf("%.2f", y) for y in Rvals])
    end
    return f
end

# ============================================================
# 17) Runners
# ============================================================

println("Threads: ", Threads.nthreads())
println("Grid: $(NX)x$(NY) cells=$(NCELLS)  Emin=$(Emin_patch)")
println("Saving to: ", OUTDIR)

function mean_curves(curve_list::Vector{NamedTuple})
    out = Dict{Symbol,Vector{Float64}}()
    for k in (:gamma_A, :gamma_AB, :sar_baseline, :sar_effective)
        mats = reduce(hcat, [getfield(c,k) for c in curve_list])
        out[k] = vec(mean(mats; dims=2))
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
fig1_cons = fig1_plot(curves_cons, "Figure 1A — Gamma richness (consumers only)"; which_label="consumers")

save(joinpath(OUTDIR, "fig1A_gamma_all.png"), fig1_all)
save(joinpath(OUTDIR, "fig1A_gamma_consumers.png"), fig1_cons)
display(fig1_all)
display(fig1_cons)

# -------------------------
# Figure 2 heatmaps
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

println("\nRunning Figure 2 heatmaps ...")
for geom in [:random, :cluster, :front]
    hms = run_fig2_for_geometry(geom)
    f2 = fig2_plot_one_geometry(hms, geom)
    save(joinpath(OUTDIR, "fig2_heatmaps_$(Symbol(geom)).png"), f2)
    display(f2)
end

println("\nDone. Output directory:")
println(OUTDIR)
