#!/usr/bin/env julia
# Julia 1.11
# Run: julia --threads auto biotic_divergence_sensitivity_plots.jl
#
# PURPOSE
# ------------------------------------------------------------
# Sensitivity plots (line charts) for A vs AB divergence metrics, holding S=250 fixed.
#
# For each sensitivity variable (x-axis), we evaluate a small set of representative
# (C, r) "cells" (4–8, far apart in parameter space) and plot the mean metric across replicates.
#
# Sensitivities included:
#  1) Grid size (NX=NY varied)            [Emin fixed absolute at baseline]
#  2) Emin absolute (LCC threshold)      [grid fixed at baseline]
#  3) Emin as fraction of landscape      [grid fixed at baseline; Emin = round(frac*NCELLS)]
#  4) SUIT_THRESH (niche threshold)      [grid & Emin fixed at baseline]
#
# Model choices match your pipeline:
# - Environment: random or autocorrelated field, rescaled to [0,100]
# - Niches: 1D Gaussian threshold mask
# - Networks: random / modular / heavy-tail / cascade
# - Mechanistic niche correlation: consumer μ correlated with mean prey μ via relaxation
# - One-prey rule: AB fixed point
# - Connectivity filter applied post-hoc, keeping ALL components >= Emin
#
# Outputs: PNG figures in OUTDIR
using Random
using Statistics
using LinearAlgebra
using CairoMakie
using Printf
using Dates

# ============================================================
# 0) Global parameters (tune here)
# ============================================================

# --- Species pool (FIXED as requested)
const S = 250
const BASAL_FRAC = 0.30

# --- Environmental field domain
const E_MIN = 0.0
const E_MAX = 100.0

# --- Environmental autocorrelation
const AUTOCORR_ITERS = 18
const AUTOCORR_ALPHA = 0.55  # 0..1 (higher = smoother)

# --- Niche suitability threshold (baseline; will be varied in sensitivity 4)
const SUIT_THRESH_BASE = 0.25

# --- Connectivity filter (baseline; will be varied in sensitivity 2/3)
const USE_CONNECTIVITY_FILTER = true
const Emin_patch_base = 60

# --- Network-family knobs (match your script)
const N_MODULES = 6
const MODULAR_IN_BIAS = 6.0
const HEAVYTAIL_GAMMA = 2.2
const HEAVYTAIL_KMAX_FRAC = 0.35
const CASCADE_LAMBDA = 2.5

# --- Mechanistic niche-correlation builder knobs
const RELAX_ITERS = 30
const MU_NOISE_SD = 1.8

# --- Replicates per point (keep modest; increase for final)
const NREP = 8

# --- Seeds
const BASE_SEED = 20260202

# --- Output directory
# ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
OUTDIR = joinpath(pwd(), "output_biotic_divergence_sensitivity_")
isdir(OUTDIR) || mkpath(OUTDIR)

# --- Which environments / networks / regimes to run (EDIT THESE)
const ENVKINDS = [:random, :autocorr]
const NETFAMS_TO_RUN = [:random, :modular, :heavytail, :cascade]
const REGIMES_TO_RUN = 1:4

# --- Representative (C, r) "cells" to plot as lines (EDIT IF YOU WANT)
# Choose 8 points spanning the space.
const CELLS_CR = [
    (0.02, 0.10),
    (0.02, 0.80),
    (0.055, 0.10),
    (0.055, 0.45),
    (0.055, 0.80),
    (0.09, 0.10),
    (0.09, 0.45),
    (0.09, 0.80),
]

# --- Sensitivity x-axes (EDIT IF YOU WANT)
const GRID_SIZES = [20, 40, 60, 80, 100]                        # NX=NY
const EMIN_ABS_VALUES = [0, 20, 40, 60, 90, 120]       # absolute Emin
const EMIN_FRACS = [0.0, 0.005, 0.01, 0.015, 0.02, 0.03]  # fraction of NCELLS
const SUIT_THRESH_VALUES = [0.15, 0.20, 0.25, 0.30, 0.35]

# ============================================================
# 1) Grid + neighbors (4-neighbour)
# ============================================================
struct Grid
    NX::Int
    NY::Int
    NCELLS::Int
    neigh4::Vector{Vector{Int}}
end

@inline linidx(x::Int, y::Int, NX::Int) = (y - 1) * NX + x
@inline x_of(i::Int, NX::Int) = ((i - 1) % NX) + 1
@inline y_of(i::Int, NX::Int) = ((i - 1) ÷ NX) + 1

function build_grid(NX::Int, NY::Int)
    nc = NX * NY
    neigh = Vector{Vector{Int}}(undef, nc)
    for i in 1:nc
        x = x_of(i, NX); y = y_of(i, NX)
        nb = Int[]
        if x > 1;  push!(nb, linidx(x-1,y,NX)); end
        if x < NX; push!(nb, linidx(x+1,y,NX)); end
        if y > 1;  push!(nb, linidx(x,y-1,NX)); end
        if y < NY; push!(nb, linidx(x,y+1,NX)); end
        neigh[i] = nb
    end
    return Grid(NX, NY, nc, neigh)
end

# ============================================================
# 2) Connectivity utilities (keep ALL components >= Emin)
# ============================================================
mutable struct CCWorkspace
    seen::Vector{Int32}
    stamp::Int32
    queue::Vector{Int}
end

function make_workspace(ncells::Int)
    return CCWorkspace(fill(Int32(0), ncells), Int32(0), Int[])
end

function keep_components_ge_Emin!(ws::CCWorkspace, grid::Grid, mask::BitVector, Emin::Int, out::BitVector)
    fill!(out, false)
    count(mask) == 0 && return out

    ws.stamp += 1
    stamp = ws.stamp
    seen = ws.seen
    empty!(ws.queue)

    comp_nodes = Int[]

    @inbounds for i in 1:grid.NCELLS
        if mask[i] && seen[i] != stamp
            empty!(comp_nodes)

            # BFS collect component
            seen[i] = stamp
            push!(ws.queue, i)
            qpos = 1

            while qpos <= length(ws.queue)
                v = ws.queue[qpos]; qpos += 1
                push!(comp_nodes, v)
                for nb in grid.neigh4[v]
                    if mask[nb] && seen[nb] != stamp
                        seen[nb] = stamp
                        push!(ws.queue, nb)
                    end
                end
            end

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

function apply_connectivity_filter!(ws::CCWorkspace, grid::Grid, mask::BitVector, Emin::Int, tmpout::BitVector)
    if !USE_CONNECTIVITY_FILTER || Emin <= 0
        return mask
    end
    if count(mask) < Emin
        fill!(tmpout, false)
        return tmpout
    end
    keep_components_ge_Emin!(ws, grid, mask, Emin, tmpout)
    return tmpout
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

function smooth_field_once!(grid::Grid, E::Vector{Float64}, tmp::Vector{Float64}, α::Float64)
    @inbounds for i in 1:grid.NCELLS
        m = E[i]
        nb = grid.neigh4[i]
        for j in nb
            m += E[j]
        end
        m /= (1 + length(nb))
        tmp[i] = (1-α)*E[i] + α*m
    end
    copyto!(E, tmp)
    return E
end

function make_environment(rng::AbstractRNG, grid::Grid, kind::Symbol)
    E = randn(rng, grid.NCELLS)
    if kind == :autocorr
        tmp = similar(E)
        for _ in 1:AUTOCORR_ITERS
            smooth_field_once!(grid, E, tmp, AUTOCORR_ALPHA)
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
    m = BitVector(undef, length(E))
    @inbounds for i in eachindex(E)
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

function draw_sigmas(rng::AbstractRNG, regime::BreadthRegime, Sp::Int)
    σ = Vector{Float64}(undef, Sp)
    for i in 1:Sp
        val = exp(log(regime.meanσ) + regime.logsd * randn(rng))
        σ[i] = clamp(val, 1.5, 45.0)
    end
    return σ
end

# ============================================================
# 5) Network builders (4 families)
# ============================================================
function consumers_and_basal(Sp::Int)
    nb = round(Int, BASAL_FRAC * Sp)
    basal_mask = BitVector(falses(Sp))
    basal_mask[1:nb] .= true
    consumers = collect((nb+1):Sp)
    return nb, basal_mask, consumers
end

function ensure_min1_prey!(rng::AbstractRNG, prey::Vector{Vector{Int}}, basal_mask::BitVector, Sp::Int)
    for i in 1:Sp
        basal_mask[i] && continue
        if isempty(prey[i])
            candidates = findall(basal_mask)
            if isempty(candidates)
                cand = [j for j in 1:Sp if j != i]
                push!(prey[i], cand[rand(rng, 1:length(cand))])
            else
                push!(prey[i], candidates[rand(rng, 1:length(candidates))])
            end
        end
    end
end

function build_metaweb_random(rng::AbstractRNG, C::Float64, basal_mask::BitVector, Sp::Int)
    prey = [Int[] for _ in 1:Sp]
    consumers = findall(!, basal_mask)
    Ltarget = round(Int, C * Sp^2)

    for i in consumers
        cand = [j for j in 1:Sp if j != i]
        push!(prey[i], cand[rand(rng, 1:length(cand))])
    end

    L = length(consumers)
    while L < Ltarget
        i = consumers[rand(rng, 1:length(consumers))]
        j = rand(rng, 1:Sp)
        (j == i) && continue
        if j ∉ prey[i]
            push!(prey[i], j)
            L += 1
        end
    end
    return prey
end

function build_metaweb_modular(rng::AbstractRNG, C::Float64, basal_mask::BitVector, Sp::Int)
    prey = [Int[] for _ in 1:Sp]
    consumers = findall(!, basal_mask)

    MODULE = Vector{Int}(undef, Sp)
    for i in 1:Sp
        MODULE[i] = 1 + (i - 1) % N_MODULES
    end

    Ltarget = round(Int, C * Sp^2)

    function sample_prey(i::Int)
        mi = MODULE[i]
        inmod = Int[]
        outmod = Int[]
        for j in 1:Sp
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
        if j ∉ prey[i]
            push!(prey[i], j)
            L += 1
        end
    end
    return prey
end

function build_metaweb_heavytail(rng::AbstractRNG, C::Float64, basal_mask::BitVector, Sp::Int)
    prey = [Int[] for _ in 1:Sp]
    consumers = findall(!, basal_mask)
    nC = length(consumers)
    Ltarget = round(Int, C * Sp^2)

    w = zeros(Float64, nC)
    for k in 1:nC
        u = rand(rng)
        w[k] = u^(-1/(HEAVYTAIL_GAMMA-1))
    end
    w ./= sum(w)

    deg = ones(Int, nC)
    remaining = max(0, Ltarget - nC)
    for _ in 1:remaining
        rrr = rand(rng)
        acc = 0.0
        idx = 1
        for k in 1:nC
            acc += w[k]
            if rrr <= acc
                idx = k
                break
            end
        end
        deg[idx] += 1
    end

    kmax = max(2, round(Int, HEAVYTAIL_KMAX_FRAC * (Sp-1)))
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
        cand = [j for j in 1:Sp if j != i]
        shuffle!(rng, cand)
        d = min(deg[kk], length(cand))
        prey[i] = cand[1:d]
    end

    ensure_min1_prey!(rng, prey, basal_mask, Sp)
    return prey
end

function build_metaweb_cascade(rng::AbstractRNG, C::Float64, basal_mask::BitVector, Sp::Int)
    prey = [Int[] for _ in 1:Sp]
    consumers = findall(!, basal_mask)
    nC = length(consumers)
    Ltarget = round(Int, C * Sp^2)

    ranks = zeros(Float64, Sp)
    for attempt in 1:200
        for i in 1:Sp
            ranks[i] = rand(rng)
        end
        ok = true
        for i in consumers
            lower = findall(j -> ranks[j] < ranks[i] && j != i, 1:Sp)
            isempty(lower) && (ok = false; break)
        end
        ok && break
        attempt == 200 && error("Failed to sample cascade ranks")
    end

    function sample_lower_prey(i::Int)
        ri = ranks[i]
        lower = Int[]
        w = Float64[]
        for j in 1:Sp
            (j == i) && continue
            if ranks[j] < ri
                push!(lower, j)
                push!(w, exp(-CASCADE_LAMBDA * (ri - ranks[j])))
            end
        end
        sw = sum(w)
        rrr = rand(rng) * sw
        acc = 0.0
        for k in 1:length(lower)
            acc += w[k]
            if rrr <= acc
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
        if j ∉ prey[i]
            push!(prey[i], j)
            L += 1
        end
    end
    return prey
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

function prey_means(mu::Vector{Float64}, prey::Vector{Vector{Int}}, basal_mask::BitVector, Sp::Int)
    pm = fill(NaN, Sp)
    for i in 1:Sp
        if basal_mask[i] || isempty(prey[i])
            continue
        end
        pm[i] = mean(mu[prey[i]])
    end
    return pm
end

function mechanistic_corr(mu::Vector{Float64}, prey::Vector{Vector{Int}}, basal_mask::BitVector, Sp::Int)
    pm = prey_means(mu, prey, basal_mask, Sp)
    idx = findall(i -> !basal_mask[i] && !isnan(pm[i]), 1:Sp)
    length(idx) < 6 && return 0.0
    return pearson_r(mu[idx], pm[idx])
end

function assign_mus_with_target_corr!(
    rng::AbstractRNG,
    mu::Vector{Float64},
    prey::Vector{Vector{Int}},
    basal_mask::BitVector,
    target_r::Float64,
    Sp::Int;
    relax_iters::Int = RELAX_ITERS
)
    consumers = findall(!, basal_mask)

    function relax(alpha::Float64)
        m = copy(mu)
        for _ in 1:relax_iters
            pm = prey_means(m, prey, basal_mask, Sp)
            @inbounds for i in consumers
                isnan(pm[i]) && continue
                m[i] = clamp((1-alpha)*m[i] + alpha*pm[i] + MU_NOISE_SD*randn(rng), E_MIN, E_MAX)
            end
        end
        return mechanistic_corr(m, prey, basal_mask, Sp), m
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
function fixed_point_AB(A_pres::Vector{BitVector}, prey::Vector{Vector{Int}}, basal_mask::BitVector, ncells::Int)
    Sp = length(A_pres)  # <- infer S from inputs

    pres = [copy(A_pres[i]) for i in 1:Sp]
    newp = [BitVector(falses(ncells)) for _ in 1:Sp]

    changed = true
    iter = 0
    while changed
        iter += 1
        changed = false

        for i in 1:Sp
            if basal_mask[i] || isempty(prey[i])
                newp[i] = pres[i]
            else
                u = BitVector(falses(ncells))
                @inbounds for j in prey[i]
                    u .|= pres[j]
                end
                newp[i] = A_pres[i] .& u
            end
        end

        for i in 1:Sp
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
# 8) Metrics (consumers-only)
# ============================================================
function gamma_richness_cons(pres::Vector{BitVector}, basal_mask::BitVector, Sp::Int)
    c = 0
    for i in 1:Sp
        basal_mask[i] && continue
        c += (count(pres[i]) > 0) ? 1 : 0
    end
    return c
end

function mean_jaccard_mismatch(A::Vector{BitVector}, AB::Vector{BitVector}, basal_mask::BitVector, Sp::Int)
    vals = Float64[]
    for i in 1:Sp
        basal_mask[i] && continue
        Ai = A[i]
        count(Ai) == 0 && continue
        ABi = AB[i]
        inter = count(Ai .& ABi)
        uni   = count(Ai .| ABi)
        J = uni == 0 ? 1.0 : (inter / uni)
        push!(vals, 1 - J)
    end
    return isempty(vals) ? NaN : mean(vals)
end

function frac_affected(A::Vector{BitVector}, AB::Vector{BitVector}, basal_mask::BitVector, Sp::Int)
    num = 0
    den = 0
    for i in 1:Sp
        basal_mask[i] && continue
        Ai = A[i]
        count(Ai) == 0 && continue
        den += 1
        num += (Ai != AB[i]) ? 1 : 0
    end
    return den == 0 ? NaN : num / den
end

# ============================================================
# 9) One replicate at (grid, envkind, netfamily, regime, C, target_r, Emin, suit_thresh)
# ============================================================
function simulate_one!(rng::AbstractRNG, grid::Grid, ws::CCWorkspace,
                       envkind::Symbol, netfamily::Symbol, regime::BreadthRegime,
                       C::Float64, target_r::Float64, Emin::Int, suit_thresh::Float64,
                       Sp::Int)

    nb, basal_mask, _ = consumers_and_basal(Sp)

    E = make_environment(rng, grid, envkind)

    prey = if netfamily == :random
        build_metaweb_random(rng, C, basal_mask, Sp)
    elseif netfamily == :modular
        build_metaweb_modular(rng, C, basal_mask, Sp)
    elseif netfamily == :heavytail
        build_metaweb_heavytail(rng, C, basal_mask, Sp)
    elseif netfamily == :cascade
        build_metaweb_cascade(rng, C, basal_mask, Sp)
    else
        error("Unknown netfamily: $netfamily")
    end

    σ = draw_sigmas(rng, regime, Sp)
    μ = Vector{Float64}(undef, Sp)
    for i in 1:Sp
        μ[i] = rand(rng) * (E_MAX - E_MIN) + E_MIN
    end

    achieved_r = assign_mus_with_target_corr!(rng, μ, prey, basal_mask, target_r, Sp)

    A_raw = Vector{BitVector}(undef, Sp)
    for i in 1:Sp
        A_raw[i] = suitability_mask_1d(E, μ[i], σ[i], suit_thresh)
    end

    AB_raw = fixed_point_AB(A_raw, prey, basal_mask, grid.NCELLS)

    A  = Vector{BitVector}(undef, Sp)
    AB = Vector{BitVector}(undef, Sp)
    tmp = BitVector(falses(grid.NCELLS))

    for i in 1:Sp
        Af = apply_connectivity_filter!(ws, grid, A_raw[i], Emin, tmp)
        A[i] = (Af === tmp) ? copy(tmp) : A_raw[i]

        ABf = apply_connectivity_filter!(ws, grid, AB_raw[i], Emin, tmp)
        AB[i] = (ABf === tmp) ? copy(tmp) : AB_raw[i]
    end

    SA  = gamma_richness_cons(A, basal_mask, Sp)
    SAB = gamma_richness_cons(AB, basal_mask, Sp)
    dSrel = (SA == 0) ? NaN : (1.0 - (SAB / SA))

    mjm = mean_jaccard_mismatch(A, AB, basal_mask, Sp)
    fa  = frac_affected(A, AB, basal_mask, Sp)

    return (dSrel=dSrel, mean_jaccard_mismatch=mjm, frac_affected=fa, achieved_r=achieved_r)
end

# ============================================================
# 10) Sensitivity runner
# ============================================================
const NETNAMES = Dict(
    :random=>"Random",
    :modular=>"Modular",
    :heavytail=>"Heavy-tail",
    :cascade=>"Cascade"
)

function cell_label(C::Float64, r::Float64)
    return @sprintf("C=%.3f, r=%.2f", C, r)
end

# Deterministic seed per condition (so lines are stable across runs)
function seed_for(rep::Int; env::Symbol, net::Symbol, reg_i::Int, C::Float64, r::Float64,
                  NX::Int, Emin::Int, suit::Float64, tag::Int)
    # tag distinguishes which sensitivity experiment (so they don't reuse identical streams)
    # fold reals into ints
    cI = round(Int, 10_000 * C)
    rI = round(Int, 10_000 * r)
    sI = round(Int, 10_000 * suit)
    return BASE_SEED +
           1_000_000 * tag +
           100_000   * (findfirst(==(env), ENVKINDS)) +
           10_000    * (findfirst(==(net), NETFAMS_TO_RUN)) +
           1_000     * reg_i +
           10        * NX +
           3         * Emin +
           7         * cI +
           11        * rI +
           13        * sI +
           rep
end

# Run one sensitivity experiment:
# xvals: vector of x values
# cfg(x) returns (grid::Grid, Emin::Int, suit_thresh::Float64)
function run_sensitivity(xvals, cfg; tag::Int, label_x::String)
    # Output structure:
    # results[(env, net, reg_i, metric, cell_index)] = Vector{Float64}(length(xvals))
    metrics = [:dSrel, :mean_jaccard_mismatch, :frac_affected, :achieved_r]
    results = Dict{Tuple{Symbol,Symbol,Int,Symbol,Int}, Vector{Float64}}()

    for env in ENVKINDS, net in NETFAMS_TO_RUN, reg_i in REGIMES_TO_RUN, m in metrics, ci in 1:length(CELLS_CR)
        results[(env, net, reg_i, m, ci)] = fill(NaN, length(xvals))
    end

    for (xi, x) in enumerate(xvals)
        grid, Emin, suit_thresh = cfg(x)
        ws = make_workspace(grid.NCELLS)
        rng = MersenneTwister(BASE_SEED + 17)

        for env in ENVKINDS, net in NETFAMS_TO_RUN, reg_i in REGIMES_TO_RUN
            reg = regimes[reg_i]

            for (ci, (C, r)) in enumerate(CELLS_CR)
                dS = Float64[]; jm = Float64[]; fa = Float64[]; ar = Float64[]

                for rep in 1:NREP
                    seed = seed_for(rep; env=env, net=net, reg_i=reg_i, C=C, r=r,
                                    NX=grid.NX, Emin=Emin, suit=suit_thresh, tag=tag)
                    Random.seed!(rng, seed)
                    out = simulate_one!(rng, grid, ws, env, net, reg, C, r, Emin, suit_thresh)
                    push!(dS, out.dSrel)
                    push!(jm, out.mean_jaccard_mismatch)
                    push!(fa, out.frac_affected)
                    push!(ar, out.achieved_r)
                end

                results[(env, net, reg_i, :dSrel, ci)][xi] = mean(dS)
                results[(env, net, reg_i, :mean_jaccard_mismatch, ci)][xi] = mean(jm)
                results[(env, net, reg_i, :frac_affected, ci)][xi] = mean(fa)
                results[(env, net, reg_i, :achieved_r, ci)][xi] = mean(ar)
            end
        end

        println(@sprintf("Done x=%s (%d/%d)  Emin=%d  suit=%.3f  grid=%dx%d",
                         string(x), xi, length(xvals), Emin, suit_thresh, grid.NX, grid.NY))
    end

    return results
end

function run_sensitivity_threaded(xvals, cfg; tag::Int, label_x::String)
    metrics = (:dSrel, :mean_jaccard_mismatch, :frac_affected, :achieved_r)

    # Preallocate results exactly like before
    results = Dict{Tuple{Symbol,Symbol,Int,Symbol,Int}, Vector{Float64}}()
    for env in ENVKINDS, net in NETFAMS_TO_RUN, reg_i in REGIMES_TO_RUN, m in metrics, ci in 1:length(CELLS_CR)
        results[(env, net, reg_i, m, ci)] = fill(NaN, length(xvals))
    end

    # Build a static list of all parameter combinations except x (we’ll thread over the full grid of jobs)
    combos = collect(Iterators.product(ENVKINDS, NETFAMS_TO_RUN, REGIMES_TO_RUN, 1:length(CELLS_CR)))

    # Main loop over x (sequential) so your progress printing stays nice and cfg() is called once per x
    for (xi, x) in enumerate(xvals)
        grid, Emin, suit_thresh, Sp = cfg(x)

        # Thread-local accumulation arrays for this xi (so threads do not touch `results` while working)
        # dims: env × net × reg × ci  (we’ll store four metrics in separate arrays)
        nenv = length(ENVKINDS)
        nnet = length(NETFAMS_TO_RUN)
        nreg = length(REGIMES_TO_RUN)
        nci  = length(CELLS_CR)

        dS_mean = Array{Float64}(undef, nenv, nnet, nreg, nci)
        jm_mean = Array{Float64}(undef, nenv, nnet, nreg, nci)
        fa_mean = Array{Float64}(undef, nenv, nnet, nreg, nci)
        ar_mean = Array{Float64}(undef, nenv, nnet, nreg, nci)

        # Create index maps so we can store into dense arrays
        env_to_i = Dict(env => i for (i, env) in enumerate(ENVKINDS))
        net_to_i = Dict(net => i for (i, net) in enumerate(NETFAMS_TO_RUN))
        reg_to_i = Dict(reg => i for (i, reg) in enumerate(REGIMES_TO_RUN))

        # Threaded over (env, net, reg_i, ci) combos for this fixed x/grid
        Threads.@threads for k in eachindex(combos)
            env, net, reg_i, ci = combos[k]
            reg = regimes[reg_i]
            C, r = CELLS_CR[ci]

            # Each job gets its own workspace and RNG -> thread-safe
            ws  = make_workspace(grid.NCELLS)
            rng = MersenneTwister()  # independent per job

            # Online mean (avoid allocating dS/jm/fa/ar vectors)
            mdS = 0.0; mjm = 0.0; mfa = 0.0; mar = 0.0

            for rep in 1:NREP
                seed = seed_for(rep; env=env, net=net, reg_i=reg_i, C=C, r=r,
                                NX=grid.NX, Emin=Emin, suit=suit_thresh, tag=tag)
                Random.seed!(rng, seed)

                out = simulate_one!(rng, grid, ws, env, net, reg, C, r, Emin, suit_thresh, Sp)

                # numerically stable online mean update
                t = rep
                mdS += (out.dSrel - mdS) / t
                mjm += (out.mean_jaccard_mismatch - mjm) / t
                mfa += (out.frac_affected - mfa) / t
                mar += (out.achieved_r - mar) / t
            end

            ie = env_to_i[env]
            in = net_to_i[net]
            ir = reg_to_i[reg_i]

            dS_mean[ie,in,ir,ci] = mdS
            jm_mean[ie,in,ir,ci] = mjm
            fa_mean[ie,in,ir,ci] = mfa
            ar_mean[ie,in,ir,ci] = mar
        end

        # Commit results for this xi sequentially (no thread contention on Dict/arrays)
        for env in ENVKINDS, net in NETFAMS_TO_RUN, reg_i in REGIMES_TO_RUN, ci in 1:length(CELLS_CR)
            ie = env_to_i[env]; in = net_to_i[net]; ir = reg_to_i[reg_i]
            results[(env, net, reg_i, :dSrel, ci)][xi] = dS_mean[ie,in,ir,ci]
            results[(env, net, reg_i, :mean_jaccard_mismatch, ci)][xi] = jm_mean[ie,in,ir,ci]
            results[(env, net, reg_i, :frac_affected, ci)][xi] = fa_mean[ie,in,ir,ci]
            results[(env, net, reg_i, :achieved_r, ci)][xi] = ar_mean[ie,in,ir,ci]
        end

        println(@sprintf("Done x=%s (%d/%d)  Emin=%d  suit=%.3f  grid=%dx%d",
                         string(x), xi, length(xvals), Emin, suit_thresh, grid.NX, grid.NY))
    end

    return results
end

# ============================================================
# 11) Plotting: line plots (3 metrics + achieved_r diagnostic)
# ============================================================
function axislegend_outside_right!(gl::GridLayout, ax;
                                   row::Int,
                                   col_axis::Int=1,
                                   col_legend::Int=2,
                                   nbanks::Int=2,
                                   framevisible::Bool=false,
                                   title::String="")

    # Make sure column 2 exists and has room (safe to call repeatedly)
    try
        colsize!(gl, col_axis, Relative(0.72))
        colsize!(gl, col_legend, Relative(0.28))
    catch
        # ignore if colsize! not available in your Makie version
    end

    # Create legend from axis' labeled plots
    leg = Legend(gl[row, col_legend], ax; nbanks=nbanks, framevisible=framevisible, title=title)

    # Optional nicer alignment
    leg.tellheight = false
    leg.tellwidth  = true

    return leg
end


function plot_sensitivity_lines(xvals, results; expname::String, xlabel::String)
    metrics_main = [
        # (:dSrel, "Relative richness loss (1 - S_AB/S_A)"),
        (:mean_jaccard_mismatch, "Mean Jaccard mismatch (mean(1 - J))"),
        # (:frac_affected, "Fraction affected (A_i != AB_i)"),
    ]
    # metric_diag = (:achieved_r, "Achieved mechanistic niche corr")

    for env in ENVKINDS, net in NETFAMS_TO_RUN, reg_i in REGIMES_TO_RUN
        regname = regimes[reg_i].name
        netname = NETNAMES[net]
        envname = (env == :random) ? "Random env" : "Autocorr env"

        f = Figure(size=(1200, 620))
        gl = f[1,1] = GridLayout()
        rowgap!(gl, 12); colgap!(gl, 10)

        Label(gl[1,1], "$(expname) — $(envname) — $(netname) — $(regname)",
              fontsize=12, halign=:left)

        # Main 3 metrics
        for (pi, (m, ttl)) in enumerate(metrics_main)
            ax = Axis(gl[pi+1, 1], xlabel=xlabel, ylabel=string(m), title=ttl)

            plots = Any[]   # collect the line objects for the legend

            for ci in 1:length(CELLS_CR)
                y = results[(env, net, reg_i, m, ci)]
                p = lines!(ax, xvals, y, label=cell_label(CELLS_CR[ci]...))
                push!(plots, p)
                ylims!(ax, 0, 1)
            end

            axislegend_outside_right!(gl, ax; row=pi+1, nbanks=2, framevisible=false)
        end

        # # Achieved r diagnostic
        # axd = Axis(gl[5, 1], xlabel=xlabel, ylabel="achieved_r", title=metric_diag[2])
        # for ci in 1:length(CELLS_CR)
        #     y = results[(env, net, reg_i, metric_diag[1], ci)]
        #     lines!(axd, xvals, y, label=cell_label(CELLS_CR[ci]...))
        # end
        # axislegend(axd; position=:rb, nbanks=2, framevisible=false)

        fname = replace(lowercase(expname), " " => "_")
        out = joinpath(OUTDIR, "sens_$(fname)_$(env)_$(net)_reg$(reg_i).png")
        save(out, f)
        display(f)
    end
end

# ============================================================
# 12) Run all requested sensitivity experiments
# ============================================================
println("OUTDIR: ", OUTDIR)
println("S base: ", S, "  basal_frac=", BASAL_FRAC, "  reps=", NREP)
println("Cells (C,r): ", CELLS_CR)
println("Running envs: ", ENVKINDS, "  nets: ", NETFAMS_TO_RUN, "  regimes: ", collect(REGIMES_TO_RUN))
println()

# Baseline grid for non-grid sensitivities
const GRID_BASE = 60

# Sensitivity 1: Grid size (Emin fixed absolute at baseline; suit_thresh baseline)
# res_grid = run_sensitivity_threaded(GRID_SIZES, x -> begin
#         NX = Int(x); NY = Int(x)
#         grid = build_grid(NX, NY)
#         Emin = Emin_patch_base
#         suit = SUIT_THRESH_BASE
#         Sp   = S_BASE
#         return (grid, Emin, suit, Sp)
#     end; tag=1, label_x="grid")

using Serialization
# serialize(joinpath(OUTDIR, "sens_grid_size.jls"), res_grid)

res_grid = deserialize(joinpath(OUTDIR, "sens_grid_size.jls"))
plot_sensitivity_lines(
    GRID_SIZES, res_grid;
    expname="Grid size sensitivity (Emin fixed)",
    xlabel="Grid size (NX=NY)"
)

# Sensitivity 2: Emin absolute (grid fixed; suit_thresh baseline)
# res_emin_abs = run_sensitivity_threaded(EMIN_ABS_VALUES, x -> begin
#         grid = build_grid(GRID_BASE, GRID_BASE)
#         Emin = Int(x)
#         suit = SUIT_THRESH_BASE,
#         Sp   = S_BASE
#         return (grid, Emin, suit)
#     end; tag=2, label_x="emin_abs")
# serialize(joinpath(OUTDIR, "sens_emin_abs.jls"), res_emin_abs)

res_emin_abs = deserialize(joinpath(OUTDIR, "sens_emin_abs.jls"))
plot_sensitivity_lines(EMIN_ABS_VALUES, res_emin_abs;
    expname="Emin absolute sensitivity (grid fixed)",
    xlabel="Emin (LCC threshold, cells)")

# Sensitivity 3: Emin proportional (grid fixed; suit_thresh baseline)
# res_emin_frac = run_sensitivity_threaded(EMIN_FRACS, x -> begin
#         grid = build_grid(GRID_BASE, GRID_BASE)
#         Emin = round(Int, Float64(x) * grid.NCELLS)
#         suit = SUIT_THRESH_BASE,
#         Sp   = S_BASE
#         return (grid, Emin, suit)
#     end; tag=3, label_x="emin_frac")
# serialize(joinpath(OUTDIR, "sens_emin_frac.jls"), res_emin_frac)

res_emin_frac = deserialize(joinpath(OUTDIR, "sens_emin_frac.jls"))
plot_sensitivity_lines(EMIN_FRACS, res_emin_frac;
    expname="Emin fraction sensitivity (grid fixed)",
    xlabel="Emin fraction of landscape (Emin = frac * NCELLS)")

# Sensitivity 4: SUIT_THRESH (grid & Emin fixed)
# res_suit = run_sensitivity_threaded(SUIT_THRESH_VALUES, x -> begin
#         grid = build_grid(GRID_BASE, GRID_BASE)
#         Emin = Emin_patch_base
#         suit = Float64(x),
#         Sp   = S_BASE
#         return (grid, Emin, suit)
#     end; tag=4, label_x="suit_thresh")
# serialize(joinpath(OUTDIR, "sens_suit_thresh.jls"), res_suit)

res_suit = deserialize(joinpath(OUTDIR, "sens_suit_thresh.jls"))
plot_sensitivity_lines(SUIT_THRESH_VALUES, res_suit;
    expname="Suitability threshold sensitivity (grid & Emin fixed)",
    xlabel="SUIT_THRESH")

# Sensitivity 5: Species per cell (grid & Emin fixed)
const S_VALUES = 150:50:400

# res_S = run_sensitivity_threaded(S_VALUES, x -> begin
#         grid = build_grid(GRID_BASE, GRID_BASE)
#         Emin = Emin_patch_base
#         suit = SUIT_THRESH_BASE
#         Sp   = Int(x)
#         return (grid, Emin, suit, Sp)
#     end; tag=5, label_x="S")
# serialize(joinpath(OUTDIR, "sens_S.jls"), res_S)

res_S = deserialize(joinpath(OUTDIR, "sens_S.jls"))
plot_sensitivity_lines(S_VALUES, res_S;
    expname="Species richness sensitivity (grid & Emin fixed)",
    xlabel="S")

println("\nDone. Figures saved to:")
println(OUTDIR)

######## EXTRAAA
# ============================================================
# QUALITATIVE INVARIANCE TEST for Sensitivity 4 (SUIT_THRESH)
# Checks whether regime ordering changes across suit_thresh values
# using mean_jaccard_mismatch aggregated over CELLS_CR.
# ============================================================

# res_suit might be a Dict already, or a struct holding results in .results
results_suit = hasproperty(res_suit, :results) ? getproperty(res_suit, :results) : res_suit

const METRIC_Q = :mean_jaccard_mismatch

# Find baseline index (closest value in SUIT_THRESH_VALUES)
ib = argmin(abs.(SUIT_THRESH_VALUES .- SUIT_THRESH_BASE))
@info "Baseline suit_thresh" SUIT_THRESH_BASE "closest grid value" SUIT_THRESH_VALUES[ib] "index" ib

# mean across CELLS_CR at a given threshold index k
function regime_score_at_k(results, env, net, reg_i, metric, k)
    vals = Float64[]
    for ci in 1:length(CELLS_CR)
        y = results[(env, net, reg_i, metric, ci)]
        v = y[k]
        isfinite(v) && push!(vals, v)
    end
    return isempty(vals) ? NaN : mean(vals)
end

# rank regimes by score at threshold index k (higher mismatch = "worse" = rank 1)
function regime_order_at_k(results, env, net, regset, metric, k)
    scores = [(reg_i, regime_score_at_k(results, env, net, reg_i, metric, k)) for reg_i in regset]
    # drop NaNs if any
    scores = filter(x -> isfinite(x[2]), scores)
    sort!(scores, by = x -> x[2], rev=true)
    return [x[1] for x in scores], Dict(x[1] => x[2] for x in scores)
end

# Kendall-style agreement score in [-1,1] for two orderings over the same set
function kendall_tau(orderA::Vector{Int}, orderB::Vector{Int})
    # assumes both contain same items (possibly subset if NaNs dropped)
    posA = Dict(orderA[i] => i for i in eachindex(orderA))
    posB = Dict(orderB[i] => i for i in eachindex(orderB))
    items = intersect(keys(posA), keys(posB)) |> collect

    n = length(items)
    n < 2 && return NaN

    conc = 0
    disc = 0
    for i in 1:n-1, j in i+1:n
        a = items[i]; b = items[j]
        sa = posA[a] < posA[b]
        sb = posB[a] < posB[b]
        if sa == sb
            conc += 1
        else
            disc += 1
        end
    end
    denom = conc + disc
    return denom == 0 ? NaN : (conc - disc) / denom
end

# Collect outputs for inspection / CSV
using DataFrames
rows = DataFrame(env=String[], net=String[], suit=Float64[],
                 reg=String[], score=Float64[], rank=Int[],
                 baseline_same=Bool[], tau_vs_base=Float64[])

for env in ENVKINDS, net in NETFAMS_TO_RUN

    base_order, base_scores = regime_order_at_k(results_suit, env, net, REGIMES_TO_RUN, METRIC_Q, ib)

    @info "Baseline regime order (worst→best)" env net [(regimes[r].name, base_scores[r]) for r in base_order]

    n_changed = 0
    changed_suits = Float64[]

    for k in eachindex(SUIT_THRESH_VALUES)
        suit = SUIT_THRESH_VALUES[k]

        ord, sc = regime_order_at_k(results_suit, env, net, REGIMES_TO_RUN, METRIC_Q, k)

        # Compare only on common regimes (in case some NaNs drop regimes at some k)
        common = intersect(base_order, ord)
        base_common = [r for r in base_order if r in common]
        ord_common  = [r for r in ord if r in common]

        same = (base_common == ord_common)
        τ = kendall_tau(base_common, ord_common)

        if !same
            n_changed += 1
            push!(changed_suits, suit)
        end

        # store per-regime ranks/scores for CSV
        # build rank dict for current ordering
        rankpos = Dict(ord[i] => i for i in eachindex(ord))
        for r in ord
            push!(rows, (String(env), String(net), Float64(suit),
                         regimes[r].name, sc[r], rankpos[r],
                         same, τ))
        end
    end

    @info "Regime-order changes across suit_thresh grid" env net n_changed
    if n_changed > 0
        @info "Threshold values with order changes" env net changed_suits
    end
end

CSV.write(joinpath(OUTDIR, "qualitative_invariance_suit_thresh.csv"), rows)
@info "Saved invariance table" joinpath(OUTDIR, "qualitative_invariance_suit_thresh.csv")

# ============================================================
# DIAGNOSTICS: what regime swaps happen vs suit_thresh?
# ============================================================
results_suit = hasproperty(res_suit, :results) ? getproperty(res_suit, :results) : res_suit
const METRIC_Q = :mean_jaccard_mismatch

# baseline index closest to SUIT_THRESH_BASE
ib = argmin(abs.(SUIT_THRESH_VALUES .- SUIT_THRESH_BASE))
base_suit = SUIT_THRESH_VALUES[ib]
@info "Baseline suit_thresh (closest in grid)" base_suit

function regime_score_at_k(results, env, net, reg_i, metric, k)
    vals = Float64[]
    for ci in 1:length(CELLS_CR)
        y = results[(env, net, reg_i, metric, ci)]
        v = y[k]
        isfinite(v) && push!(vals, v)
    end
    return isempty(vals) ? NaN : mean(vals)
end

function regime_order_scores_at_k(results, env, net, regset, metric, k)
    scores = [(reg_i, regime_score_at_k(results, env, net, reg_i, metric, k)) for reg_i in regset]
    scores = filter(x -> isfinite(x[2]), scores)
    sort!(scores, by=x->x[2], rev=true)  # worst mismatch first
    order = [x[1] for x in scores]
    sc = Dict(x[1] => x[2] for x in scores)
    return order, sc
end

# Small helper: print order with scores and gaps
function print_order(env, net, suit, order, sc)
    pairs = [(regimes[r].name, sc[r]) for r in order]
    println("\n(env=$(env), net=$(net), suit=$(round(suit,digits=4)))")
    for (i,(nm,val)) in enumerate(pairs)
        gap = (i < length(pairs)) ? (val - pairs[i+1][2]) : NaN
        println(rpad("  $(i)) $(nm)", 30), "  score=$(round(val,digits=6))",
                (isfinite(gap) ? "  gap_to_next=$(round(gap,digits=6))" : ""))
    end
end

# Identify changed suit values and print what swapped
for env in ENVKINDS, net in NETFAMS_TO_RUN

    base_order, base_sc = regime_order_scores_at_k(results_suit, env, net, REGIMES_TO_RUN, METRIC_Q, ib)
    print_order(env, net, base_suit, base_order, base_sc)

    # find indices where ordering differs from baseline (on common regimes)
    changed_idx = Int[]
    for k in eachindex(SUIT_THRESH_VALUES)
        ord, _ = regime_order_scores_at_k(results_suit, env, net, REGIMES_TO_RUN, METRIC_Q, k)
        common = intersect(base_order, ord)
        base_common = [r for r in base_order if r in common]
        ord_common  = [r for r in ord if r in common]
        (base_common == ord_common) || push!(changed_idx, k)
    end

    if isempty(changed_idx)
        println("\n  → No order changes across suit_thresh grid.")
        continue
    end

    changed_suits = SUIT_THRESH_VALUES[changed_idx]
    println("\n  → Order changes at suit_thresh values: ", changed_suits)

    # print each changed case in detail
    for k in changed_idx
        ord, sc = regime_order_scores_at_k(results_suit, env, net, REGIMES_TO_RUN, METRIC_Q, k)
        print_order(env, net, SUIT_THRESH_VALUES[k], ord, sc)
    end

    # ---- Stable window around baseline ----
    # contiguous k around ib where order matches baseline
    function order_matches_base(k)
        ord, _ = regime_order_scores_at_k(results_suit, env, net, REGIMES_TO_RUN, METRIC_Q, k)
        common = intersect(base_order, ord)
        base_common = [r for r in base_order if r in common]
        ord_common  = [r for r in ord if r in common]
        return base_common == ord_common
    end

    kL = ib
    while kL > firstindex(SUIT_THRESH_VALUES) && order_matches_base(kL-1)
        kL -= 1
    end
    kR = ib
    while kR < lastindex(SUIT_THRESH_VALUES) && order_matches_base(kR+1)
        kR += 1
    end

    println("\n  Stable window around baseline: suit in [",
            SUIT_THRESH_VALUES[kL], ", ", SUIT_THRESH_VALUES[kR], "]",
            " (indices $kL:$kR)")
end
