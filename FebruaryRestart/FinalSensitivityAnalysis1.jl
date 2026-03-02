#!/usr/bin/env julia
# Julia 1.11
# Run: julia --threads auto sensitivity_SI_summary.jl
#
# PURPOSE
# ------------------------------------------------------------
# One script that:
# 1) Runs *all* sensitivity experiments for the A vs AB baseline-divergence pipeline
#    (no habitat loss; divergence = Jaccard mismatch between A and AB distributions).
# 2) Produces an SI-ready *summary table* (CSV + pretty printed TSV) with:
#       - Δr (lowC: lowr − highr)
#       - ΔC (lowr: lowC − highC)
#       - Δr/ΔC ratio
#       - slopes vs r and vs C
#       - % monotonic-in-r checks passed
#       - Kendall τ of regime ordering vs baseline (optional but included)
#    for BOTH mean mismatch and q90 mismatch.
# 3) Produces a single “tornado plot” summarizing robustness:
#       - For each nuisance knob: range of *relative change* in median(Δr/ΔC) vs baseline.
#
# CHANGES REQUESTED
# ------------------------------------------------------------
# (1) Removed Emin ABS analysis — only Emin FRAC remains.
# (2) Fixed "S sweep with fixed mean degree" to avoid NaN / no-C-range:
#     - Compute baseline mean degree per baseline connectance value: k̂(C0)=C0*S_BASE
#     - For each richness-modified case Sp, map each baseline C0 to:
#           Cnew(C0,Sp)=k̂(C0)/Sp
#       (clamped to [1e-6, 0.5])
#     - This preserves the *range* of connectance values across CELLS_ALL while keeping mean degree fixed.

using Random
using Statistics
using LinearAlgebra
using CairoMakie
using Printf
using Dates
using DataFrames
using CSV
using Pkg

# ============================================================
# 0) BASELINE + SENSITIVITY GRID
# ============================================================

# ----- Baseline (paper) settings
const S_BASE = 250
const BASAL_FRAC = 0.30
const GRID_BASE = 60            # NX=NY
const EMIN_BASE = 60            # absolute cells (baseline used in baseline + other knobs except emin_frac)
const SUIT_BASE = 0.25
const ENVKINDS = (:random, :autocorr)
const NETFAMS  = (:random, :modular, :heavytail, :cascade)

# ----- Environmental domain
const E_MIN = 0.0
const E_MAX = 100.0
const AUTOCORR_ITERS = 18
const AUTOCORR_ALPHA = 0.55

# ----- Niche breadth regimes (match pipeline)
struct BreadthRegime
    name::String
    meanσ::Float64
    logsd::Float64
end
const REGIMES = [
    BreadthRegime("Narrow + LowVar",  7.5, 0.20),
    BreadthRegime("Narrow + HighVar", 7.5, 0.55),
    BreadthRegime("Broad + LowVar",   16.0, 0.20),
    BreadthRegime("Broad + HighVar",  16.0, 0.55),
]

# ----- Network knobs (match pipeline)
const N_MODULES = 6
const MODULAR_IN_BIAS = 6.0
const HEAVYTAIL_GAMMA = 2.2
const HEAVYTAIL_KMAX_FRAC = 0.35
const CASCADE_LAMBDA = 2.5

# ----- Mechanistic corr relaxation knobs
const RELAX_ITERS = 30
const MU_NOISE_SD = 1.8

# ----- Representative (C,r) “cells” for summaries (keep these fixed)
const CELLS_EXTREMES = (
    (0.02, 0.10),  # lowC lowr
    (0.02, 0.80),  # lowC highr
    (0.09, 0.10),  # highC lowr
    (0.09, 0.80),  # highC highr
)
const CELLS_EXTRA = (
    (0.055, 0.10),
    (0.055, 0.45),
    (0.055, 0.80),
    (0.02,  0.45),
    (0.09,  0.45),
)
const CELLS_ALL = vcat(collect(CELLS_EXTREMES), collect(CELLS_EXTRA))

# ----- Sensitivity x-axes (tune density here)
const GRID_SIZES         = [20, 40, 60, 80, 100]                    # NX=NY
const EMIN_FRACS         = [0.0, 0.005, 0.01, 0.015, 0.02, 0.03]    # Emin = frac*NCELLS
const SUIT_THRESH_VALUES = [0.15, 0.20, 0.25, 0.30, 0.35]
const S_VALUES_KHAT      = collect(150:50:400)                      # S varied, mean degree fixed (per baseline C values)

# ----- Replicates
const NREP = 8

# ----- Seeds / output
const BASE_SEED = 20260202
OUTDIR = joinpath(pwd(), "output_sensitivity_SI_")
isdir(OUTDIR) || mkpath(OUTDIR)

# ============================================================
# 1) GRID + CONNECTIVITY (keep ALL components >= Emin)
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

mutable struct CCWorkspace
    seen::Vector{Int32}
    stamp::Int32
    queue::Vector{Int}
end

make_workspace(ncells::Int) = CCWorkspace(fill(Int32(0), ncells), Int32(0), Int[])

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
    if Emin <= 0
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
# 2) ENVIRONMENT
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
# 3) NICHES + MISMATCH (mean + q90)
# ============================================================

@inline function suitability_mask_1d(E::Vector{Float64}, μ::Float64, σ::Float64, thresh::Float64)
    lim = sqrt(-2.0 * log(thresh))
    invσ = 1.0 / max(σ, 1e-6)
    m = BitVector(undef, length(E))
    @inbounds for i in eachindex(E)
        z = (E[i] - μ) * invσ
        m[i] = abs(z) <= lim
    end
    return m
end

function draw_sigmas(rng::AbstractRNG, reg::BreadthRegime, Sp::Int)
    σ = Vector{Float64}(undef, Sp)
    for i in 1:Sp
        val = exp(log(reg.meanσ) + reg.logsd * randn(rng))
        σ[i] = clamp(val, 1.5, 45.0)
    end
    return σ
end

function consumers_and_basal(Sp::Int)
    nb = round(Int, BASAL_FRAC * Sp)
    basal_mask = BitVector(falses(Sp))
    basal_mask[1:nb] .= true
    consumers = collect((nb+1):Sp)
    return nb, basal_mask, consumers
end

function per_species_mismatch(A::Vector{BitVector}, AB::Vector{BitVector}, basal_mask::BitVector)
    Sp = length(A)
    m = Float64[]
    sizehint!(m, Sp)
    for i in 1:Sp
        basal_mask[i] && continue
        Ai = A[i]
        count(Ai) == 0 && continue
        inter = count(Ai .& AB[i])
        uni   = count(Ai .| AB[i])
        J = (uni == 0) ? 1.0 : (inter / uni)
        push!(m, 1 - J)
    end
    return m
end

mismatch_stats(m::Vector{Float64}) = isempty(m) ? (mean=NaN, q90=NaN) : (mean=mean(m), q90=quantile(m, 0.90))

# ============================================================
# 4) NETWORKS (4 families)
# ============================================================

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
# 5) MECHANISTIC NICHE CORRELATION BUILDER
# ============================================================

function pearson_r(a::Vector{Float64}, b::Vector{Float64})
    ma = mean(a); mb = mean(b)
    sa = std(a); sb = std(b)
    (sa == 0 || sb == 0) && return 0.0
    return mean((a .- ma) .* (b .- mb)) / (sa * sb)
end

function prey_means(mu::Vector{Float64}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    Sp = length(mu)
    pm = fill(NaN, Sp)
    for i in 1:Sp
        if basal_mask[i] || isempty(prey[i]); continue; end
        pm[i] = mean(mu[prey[i]])
    end
    return pm
end

function mechanistic_corr(mu::Vector{Float64}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    pm = prey_means(mu, prey, basal_mask)
    idx = findall(i -> !basal_mask[i] && !isnan(pm[i]), 1:length(mu))
    length(idx) < 6 && return 0.0
    return pearson_r(mu[idx], pm[idx])
end

function assign_mus_with_target_corr!(rng::AbstractRNG, mu::Vector{Float64},
                                     prey::Vector{Vector{Int}}, basal_mask::BitVector,
                                     target_r::Float64; relax_iters::Int=RELAX_ITERS)
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
        if rmid < target_r; lo = mid; else; hi = mid; end
    end

    copyto!(mu, best_mu)
    return best_r
end

# ============================================================
# 6) AB FIXED POINT (one-prey rule)
# ============================================================

function fixed_point_AB(A_pres::Vector{BitVector}, prey::Vector{Vector{Int}}, basal_mask::BitVector, ncells::Int)
    Sp = length(A_pres)
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
# 7) ONE REPLICATE: return mean + q90 mismatch (+ achieved_r)
# ============================================================

function simulate_one!(rng::AbstractRNG, grid::Grid, ws::CCWorkspace,
                       envkind::Symbol, netfamily::Symbol, reg::BreadthRegime,
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

    σ = draw_sigmas(rng, reg, Sp)
    μ = Vector{Float64}(undef, Sp)
    for i in 1:Sp
        μ[i] = rand(rng) * (E_MAX - E_MIN) + E_MIN
    end

    achieved_r = assign_mus_with_target_corr!(rng, μ, prey, basal_mask, target_r)

    A_raw = Vector{BitVector}(undef, Sp)
    for i in 1:Sp
        A_raw[i] = suitability_mask_1d(E, μ[i], σ[i], suit_thresh)
    end

    AB_raw = fixed_point_AB(A_raw, prey, basal_mask, grid.NCELLS)

    tmp = BitVector(falses(grid.NCELLS))
    A  = Vector{BitVector}(undef, Sp)
    AB = Vector{BitVector}(undef, Sp)

    for i in 1:Sp
        Af = apply_connectivity_filter!(ws, grid, A_raw[i], Emin, tmp)
        A[i] = (Af === tmp) ? copy(tmp) : A_raw[i]

        ABf = apply_connectivity_filter!(ws, grid, AB_raw[i], Emin, tmp)
        AB[i] = (ABf === tmp) ? copy(tmp) : AB_raw[i]
    end

    m = per_species_mismatch(A, AB, basal_mask)
    st = mismatch_stats(m)

    return (mean_mismatch=st.mean, q90_mismatch=st.q90, achieved_r=achieved_r)
end

# ============================================================
# 8) SUMMARY METRICS FOR A GIVEN SETTING
# ============================================================

# Deterministic seed per condition
function seed_for(rep::Int; env::Symbol, net::Symbol, reg_i::Int, C::Float64, r::Float64,
                  NX::Int, Emin::Int, suit::Float64, Sp::Int, tag::Int)
    cI = round(Int, 10_000 * C)
    rI = round(Int, 10_000 * r)
    sI = round(Int, 10_000 * suit)
    return BASE_SEED +
           1_000_000 * tag +
           100_000   * (findfirst(==(env), ENVKINDS)) +
           10_000    * (findfirst(==(net), NETFAMS)) +
           1_000     * reg_i +
           10        * NX +
           3         * Emin +
           7         * cI +
           11        * rI +
           13        * sI +
           17        * Sp +
           rep
end

# Fit slopes over CELLS_ALL for a given (env,net,reg)
# Returns slope_r (vs r) and slope_C (vs C) from least squares.
function slopes_over_cells(CRs::Vector{Tuple{Float64,Float64}}, vals::Vector{Float64})
    # model: vals = a + bC*C + br*r
    C = Float64[x[1] for x in CRs]
    r = Float64[x[2] for x in CRs]
    y = Float64.(vals)
    good = isfinite.(y)
    if sum(good) < 5
        return (slopeC=NaN, slopeR=NaN)
    end
    X = hcat(ones(sum(good)), C[good], r[good])
    β = X \ y[good]
    return (slopeC=β[2], slopeR=β[3])
end

# Simple monotonic-in-r check:
# For each unique C among CELLS_ALL, check mean mismatch at lowr < highr (i.e., mismatch decreases with r).
function monotonic_in_r(CRs::Vector{Tuple{Float64,Float64}}, vals::Vector{Float64})
    Cs = unique([x[1] for x in CRs])
    passed = 0
    total = 0
    for C in Cs
        idx = findall(i -> CRs[i][1] == C && isfinite(vals[i]), 1:length(CRs))
        length(idx) < 2 && continue
        rr = [CRs[i][2] for i in idx]
        yy = [vals[i] for i in idx]
        imin = argmin(rr); imax = argmax(rr)
        total += 1
        if yy[imax] <= yy[imin]  # higher r => lower mismatch
            passed += 1
        end
    end
    return total == 0 ? NaN : passed / total
end

# Kendall tau for regime ordering vs baseline (per env/net)
function kendall_tau(orderA::Vector{Int}, orderB::Vector{Int})
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
        if sa == sb; conc += 1 else disc += 1 end
    end
    denom = conc + disc
    return denom == 0 ? NaN : (conc - disc) / denom
end

# Compute ordering of regimes by score (higher mismatch = "worse" first)
function regime_order(scores::Vector{Float64})
    idx = collect(1:length(scores))
    good = isfinite.(scores)
    idx = idx[good]
    sc  = scores[good]
    ord = idx[sortperm(sc; rev=true)]
    return ord
end

# Core: compute summary stats for one knob setting
# NOTE: C_mapper allows per-cell mapping of baseline C0 -> C used in simulation.
#       If C_mapper is nothing, uses baseline C0.
function summarize_setting(; tag::Int, NX::Int, Emin::Int, suit::Float64, Sp::Int,
                           C_mapper::Union{Nothing,Function}=nothing)
    grid = build_grid(NX, NX)

    nreg = length(REGIMES)
    ncells = length(CELLS_ALL)

    rows = NamedTuple[]

    regime_scores_mean = Dict{Tuple{Symbol,Symbol}, Vector{Float64}}()
    regime_scores_q90  = Dict{Tuple{Symbol,Symbol}, Vector{Float64}}()
    for env in ENVKINDS, net in NETFAMS
        regime_scores_mean[(env,net)] = fill(NaN, nreg)
        regime_scores_q90[(env,net)]  = fill(NaN, nreg)
    end

    for env in ENVKINDS, net in NETFAMS, (reg_i, reg) in enumerate(REGIMES)

        mean_cell = Vector{Float64}(undef, ncells)
        q90_cell  = Vector{Float64}(undef, ncells)

        Threads.@threads for ci in 1:ncells
            C0, r0 = CELLS_ALL[ci]
            C = isnothing(C_mapper) ? C0 : C_mapper(C0)

            ws  = make_workspace(grid.NCELLS)
            rng = MersenneTwister()

            mvals = Float64[]
            qvals = Float64[]
            sizehint!(mvals, NREP)
            sizehint!(qvals, NREP)

            for rep in 1:NREP
                seed = seed_for(rep; env=env, net=net, reg_i=reg_i,
                                C=C, r=r0, NX=NX, Emin=Emin, suit=suit, Sp=Sp, tag=tag)
                Random.seed!(rng, seed)
                out = simulate_one!(rng, grid, ws, env, net, reg, C, r0, Emin, suit, Sp)
                push!(mvals, out.mean_mismatch)
                push!(qvals, out.q90_mismatch)
            end

            mean_cell[ci] = mean(mvals)
            q90_cell[ci]  = mean(qvals)
        end

        function find_cell_idx(Ct, rt)
            for (i,(C,r)) in enumerate(CELLS_ALL)
                if isapprox(C, Ct; atol=1e-12) && isapprox(r, rt; atol=1e-12)
                    return i
                end
            end
            return nothing
        end

        i_LL = find_cell_idx(CELLS_EXTREMES[1][1], CELLS_EXTREMES[1][2])
        i_LH = find_cell_idx(CELLS_EXTREMES[2][1], CELLS_EXTREMES[2][2])
        i_HL = find_cell_idx(CELLS_EXTREMES[3][1], CELLS_EXTREMES[3][2])

        Δr_mean = (i_LL===nothing || i_LH===nothing) ? NaN : (mean_cell[i_LL] - mean_cell[i_LH])
        Δr_q90  = (i_LL===nothing || i_LH===nothing) ? NaN : (q90_cell[i_LL]  - q90_cell[i_LH])

        ΔC_mean = (i_LL===nothing || i_HL===nothing) ? NaN : (mean_cell[i_LL] - mean_cell[i_HL])
        ΔC_q90  = (i_LL===nothing || i_HL===nothing) ? NaN : (q90_cell[i_LL]  - q90_cell[i_HL])

        ratio_mean = (isfinite(Δr_mean) && isfinite(ΔC_mean) && abs(ΔC_mean)>1e-12) ? (Δr_mean/ΔC_mean) : NaN
        ratio_q90  = (isfinite(Δr_q90)  && isfinite(ΔC_q90)  && abs(ΔC_q90)>1e-12)  ? (Δr_q90/ΔC_q90)   : NaN

        sl_mean = slopes_over_cells(CELLS_ALL, mean_cell)
        sl_q90  = slopes_over_cells(CELLS_ALL, q90_cell)

        mono_mean = monotonic_in_r(CELLS_ALL, mean_cell)
        mono_q90  = monotonic_in_r(CELLS_ALL, q90_cell)

        key = (env, net)
        regime_scores_mean[key][reg_i] = mean(mean_cell[isfinite.(mean_cell)])
        regime_scores_q90[key][reg_i]  = mean(q90_cell[isfinite.(q90_cell)])

        push!(rows, (
            env=String(env),
            net=String(net),
            regime=REGIMES[reg_i].name,
            NX=NX,
            Emin=Emin,
            suit=suit,
            S=Sp,
            Δr_mean=Δr_mean,
            ΔC_mean=ΔC_mean,
            ratio_mean=ratio_mean,
            slopeC_mean=sl_mean.slopeC,
            slopeR_mean=sl_mean.slopeR,
            mono_mean=mono_mean,
            Δr_q90=Δr_q90,
            ΔC_q90=ΔC_q90,
            ratio_q90=ratio_q90,
            slopeC_q90=sl_q90.slopeC,
            slopeR_q90=sl_q90.slopeR,
            mono_q90=mono_q90
        ))
    end

    return DataFrame(rows), regime_scores_mean, regime_scores_q90
end

# ============================================================
# 9) RUN ALL SENSITIVITIES + BUILD SI TABLE
# ============================================================
println("Running BASELINE...")
base_df, base_regscores_mean, base_regscores_q90 =
    summarize_setting(tag=0, NX=GRID_BASE, Emin=EMIN_BASE, suit=SUIT_BASE, Sp=S_BASE)

base_order_mean = Dict{Tuple{Symbol,Symbol}, Vector{Int}}()
base_order_q90  = Dict{Tuple{Symbol,Symbol}, Vector{Int}}()
for env in ENVKINDS, net in NETFAMS
    base_order_mean[(env,net)] = regime_order(base_regscores_mean[(env,net)])
    base_order_q90[(env,net)]  = regime_order(base_regscores_q90[(env,net)])
end

function aggregate_SI(df::DataFrame; knob::String, knob_value)
    function med_iqr(v)
        v = v[isfinite.(v)]
        isempty(v) && return (NaN, NaN, NaN)
        q1 = quantile(v, 0.25)
        q2 = quantile(v, 0.50)
        q3 = quantile(v, 0.75)
        return (q2, q1, q3)
    end

    m_med, m_q1, m_q3 = med_iqr(df.ratio_mean)
    q_med, q_q1, q_q3 = med_iqr(df.ratio_q90)

    dr_med, _, _ = med_iqr(df.Δr_mean)
    dC_med, _, _ = med_iqr(df.ΔC_mean)

    mono_m = mean(df.mono_mean[isfinite.(df.mono_mean)])
    mono_q = mean(df.mono_q90[isfinite.(df.mono_q90)])

    return (
        knob=knob,
        value=Float64(knob_value),
        ratio_mean_median=m_med,
        ratio_mean_q25=m_q1,
        ratio_mean_q75=m_q3,
        ratio_q90_median=q_med,
        ratio_q90_q25=q_q1,
        ratio_q90_q75=q_q3,
        Δr_mean_median=dr_med,
        ΔC_mean_median=dC_med,
        mono_mean_rate=mono_m,
        mono_q90_rate=mono_q
    )
end

function tau_summary(regscores_mean, regscores_q90)
    τm = Float64[]
    τq = Float64[]
    for env in ENVKINDS, net in NETFAMS
        ord_m = regime_order(regscores_mean[(env,net)])
        ord_q = regime_order(regscores_q90[(env,net)])
        push!(τm, kendall_tau(base_order_mean[(env,net)], ord_m))
        push!(τq, kendall_tau(base_order_q90[(env,net)], ord_q))
    end
    τm = τm[isfinite.(τm)]
    τq = τq[isfinite.(τq)]
    return (tau_mean=isempty(τm) ? NaN : median(τm),
            tau_q90 =isempty(τq) ? NaN : median(τq))
end

si_rows = NamedTuple[]
append!(si_rows, [(; aggregate_SI(base_df; knob="baseline", knob_value=0.0)..., tau_mean=1.0, tau_q90=1.0)])

# ----- 1) Grid size sensitivity (Emin fixed absolute; suit fixed)
println("\nRunning GRID SIZE sensitivity...")
for NX in GRID_SIZES
    df, rsm, rsq = summarize_setting(tag=1, NX=NX, Emin=EMIN_BASE, suit=SUIT_BASE, Sp=S_BASE)
    τ = tau_summary(rsm, rsq)
    push!(si_rows, (; aggregate_SI(df; knob="grid_size", knob_value=NX)..., tau_mean=τ.tau_mean, tau_q90=τ.tau_q90))
end

# ----- 2) Emin fraction (grid fixed; suit fixed)
println("\nRunning EMIN FRAC sensitivity...")
for frac in EMIN_FRACS
    grid = build_grid(GRID_BASE, GRID_BASE)
    Emin = round(Int, frac * grid.NCELLS)
    df, rsm, rsq = summarize_setting(tag=2, NX=GRID_BASE, Emin=Emin, suit=SUIT_BASE, Sp=S_BASE)
    τ = tau_summary(rsm, rsq)
    push!(si_rows, (; aggregate_SI(df; knob="emin_frac", knob_value=frac)..., tau_mean=τ.tau_mean, tau_q90=τ.tau_q90))
end

# ----- 3) Suitability threshold (grid fixed; Emin fixed)
println("\nRunning SUIT_THRESH sensitivity...")
for suit in SUIT_THRESH_VALUES
    df, rsm, rsq = summarize_setting(tag=3, NX=GRID_BASE, Emin=EMIN_BASE, suit=suit, Sp=S_BASE)
    τ = tau_summary(rsm, rsq)
    push!(si_rows, (; aggregate_SI(df; knob="suit_thresh", knob_value=suit)..., tau_mean=τ.tau_mean, tau_q90=τ.tau_q90))
end

# ----- 4) S sweep with fixed mean degree per baseline connectance value
# For each baseline cell connectance C0, define baseline mean degree k̂(C0)=C0*S_BASE.
# For richness-modified case Sp, map C0 -> Cnew = k̂(C0)/Sp to preserve *range* of C across cells.
println("\nRunning S (mean degree fixed per baseline C) sensitivity...")
for Sp in S_VALUES_KHAT
    # define mapping C0 -> Cnew for this Sp
    C_mapper = (C0::Float64) -> begin
        khat = C0 * S_BASE
        clamp(khat / Sp, 1e-6, 0.5)
    end
    df, rsm, rsq = summarize_setting(tag=4, NX=GRID_BASE, Emin=EMIN_BASE, suit=SUIT_BASE, Sp=Sp, C_mapper=C_mapper)
    τ = tau_summary(rsm, rsq)
    push!(si_rows, (; aggregate_SI(df; knob="S_khat_fixed", knob_value=Sp)..., tau_mean=τ.tau_mean, tau_q90=τ.tau_q90))
end

# Build SI table dataframe
si = DataFrame(si_rows)

# Save outputs
csv_path = joinpath(OUTDIR, "SI_sensitivity_summary.csv")
CSV.write(csv_path, si)

tsv_path = joinpath(OUTDIR, "SI_sensitivity_summary.tsv")
open(tsv_path, "w") do io
    println(io, join(names(si), '\t'))
    for r in eachrow(si)
        println(io, join([string(r[c]) for c in names(si)], '\t'))
    end
end

println("\nSaved SI summary table:")
println("  ", csv_path)
println("  ", tsv_path)

# ============================================================
# 10) TORNADO PLOT (relative change in median ratio_mean vs baseline)
# ============================================================
base_ratio = si[(si.knob .== "baseline"), :ratio_mean_median][1]

knobs = unique(si.knob)
knobs = filter(k -> k != "baseline", knobs)

tornado_rows = NamedTuple[]
for k in knobs
    vals = si[(si.knob .== k), :ratio_mean_median]
    vals = vals[isfinite.(vals)]
    if isempty(vals) || !isfinite(base_ratio) || abs(base_ratio) < 1e-12
        push!(tornado_rows, (knob=k, lo=NaN, hi=NaN))
        continue
    end
    rel = vals ./ base_ratio .- 1.0
    push!(tornado_rows, (knob=k, lo=minimum(rel), hi=maximum(rel)))
end
tor = DataFrame(tornado_rows)

tor.width = tor.hi .- tor.lo
sort!(tor, :width, rev=true)

begin
    fig = Figure(size=(1000, 520))
    ax = Axis(fig[1,1],
        title="Sensitivity analysis: relative change in median(Δr/ΔC) vs baseline",
        xlabel="Relative change (median ratio / baseline − 1)",
        ylabel="",
        yticklabelsize=13
    )

    ypos = 1:nrow(tor)
    for (i, row) in enumerate(eachrow(tor))
        if isfinite(row.lo) && isfinite(row.hi)
            lines!(ax, [row.lo, row.hi], [i, i], linewidth=8)
            scatter!(ax, [0.0], [i], markersize=10)
        end
    end
    pretty = Dict(
        "S_khat_fixed" => "Richness (mean degree fixed)",
        "emin_frac"    => "Emin Fraction",
        "suit_thresh"  => "Suitability Threshold",
        "grid_size"    => "Grid Size",
    )

    labels = [get(pretty, String(k), String(k)) for k in tor.knob]
    ax.yticks = (collect(ypos), labels)
    vlines!(ax, [0.0], linestyle=:dash, linewidth=2)

    save(joinpath(OUTDIR, "tornado_sensitivity_ratio_mean.png"), fig)
    display(fig)
end

println("\nSaved tornado plot:")
println("  ", joinpath(OUTDIR, "tornado_sensitivity_ratio_mean.png"))

println("\nDONE. Output directory:")
println("  ", OUTDIR)