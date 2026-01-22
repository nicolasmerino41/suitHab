# suitHab_H1_H2_H3.jl
#
# FULL synthetic pipeline (static) organised by hypotheses:
#   H1: Trophic support creates systematic SAR error (AB vs A vs SAR), even static.
#   H2: Geometry matters via fragmentation of trophically supported habitat (patch criterion reveals it).
#   H3: Strength of geometry×biotic amplification is controlled by redundancy/synchrony & trophic depth,
#       WITHOUT changing the prey rule.
#
# Produces (displays, does not save):
#   Figure 1 (H1): A-only vs AB vs SAR, total-area & patch criteria, per geometry.
#   Figure 2 (H2): Mechanisms (phi, supported-LCC, frag-fail) + SAR-error vs supported-LCC scatter.
#   Figure 3a (H3): Sweep k_prey × select_sigma -> heatmaps (AUC SAR-error patch; peak frag-fail).
#   Figure 3b (H3): Sweep Lmax × basal_frac -> heatmaps (AUC SAR-error patch; peak frag-fail).
#
# Threaded + thread-safe:
#   - No shared RNG; each replicate uses deterministic per-rep seed (independent of thread schedule).
#   - No shared mutable buffers; each thread has its own scratch buffers and accumulators.
#
# Run with, e.g.:
#   JULIA_NUM_THREADS=12 julia suitHab_static_full_pipeline_H1_H2_H3_threads.jl

using Random
using Statistics
using Printf
using CairoMakie
using Base.Threads

# ----------------------------
# Deterministic per-replicate RNG (thread-schedule independent)
# ----------------------------
@inline function splitmix64(x::UInt64)
    z = x + 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

@inline function rep_rng(master_seed::Int, rep::Int)
    s = splitmix64(UInt64(master_seed) ⊻ (UInt64(rep) * 0x9e3779b97f4a7c15))
    return MersenneTwister(Int(mod(s, UInt64(typemax(Int)))))
end

# ----------------------------
# Grid indexing
# ----------------------------
@inline idx(i::Int, j::Int, n2::Int) = (i - 1) * n2 + j

# ----------------------------
# Simple smoothing (periodic) + z-score for env fields
# ----------------------------
function smooth_field!(Z::Matrix{Float64}; iters::Int=25)
    n1, n2 = size(Z)
    tmp = similar(Z)
    for _ in 1:iters
        @inbounds for i in 1:n1, j in 1:n2
            ip = (i == n1) ? 1 : i + 1
            im = (i == 1)  ? n1 : i - 1
            jp = (j == n2) ? 1 : j + 1
            jm = (j == 1)  ? n2 : j - 1
            tmp[i,j] = (Z[i,j] + Z[im,j] + Z[ip,j] + Z[i,jm] + Z[i,jp]) / 5.0
        end
        Z .= tmp
    end
    return Z
end

function zscore!(Z::Matrix{Float64})
    μ = mean(Z)
    σ = std(Z)
    Z .-= μ
    Z ./= (σ > 0 ? σ : 1.0)
    return Z
end

# ----------------------------
# Neighborhoods (4-neighborhood)
# ----------------------------
function build_neighbors(n1::Int, n2::Int)
    N = n1*n2
    up    = zeros(Int, N)
    down  = zeros(Int, N)
    left  = zeros(Int, N)
    right = zeros(Int, N)
    @inbounds for r in 1:n1, c in 1:n2
        k = idx(r,c,n2)
        up[k]    = (r > 1)  ? idx(r-1,c,n2) : 0
        down[k]  = (r < n1) ? idx(r+1,c,n2) : 0
        left[k]  = (c > 1)  ? idx(r,c-1,n2) : 0
        right[k] = (c < n2) ? idx(r,c+1,n2) : 0
    end
    return up, down, left, right
end

# Stamp-based BFS LCC on a subset
function lcc_max_size(present_flag::BitVector, present_idx::Vector{Int},
                      up::Vector{Int}, down::Vector{Int}, left::Vector{Int}, right::Vector{Int},
                      seen::Vector{Int}, queue::Vector{Int}, stamp::Base.RefValue{Int})

    isempty(present_idx) && return 0

    stamp[] += 1
    s = stamp[]
    if s == typemax(Int)
        fill!(seen, 0)
        stamp[] = 1
        s = 1
    end

    maxsz = 0
    @inbounds for start in present_idx
        if seen[start] == s
            continue
        end
        head = 1
        tail = 1
        queue[1] = start
        seen[start] = s
        sz = 0
        while head <= tail
            v = queue[head]; head += 1
            sz += 1

            nb = up[v]
            if nb != 0 && present_flag[nb] && seen[nb] != s
                tail += 1; queue[tail] = nb; seen[nb] = s
            end
            nb = down[v]
            if nb != 0 && present_flag[nb] && seen[nb] != s
                tail += 1; queue[tail] = nb; seen[nb] = s
            end
            nb = left[v]
            if nb != 0 && present_flag[nb] && seen[nb] != s
                tail += 1; queue[tail] = nb; seen[nb] = s
            end
            nb = right[v]
            if nb != 0 && present_flag[nb] && seen[nb] != s
                tail += 1; queue[tail] = nb; seen[nb] = s
            end
        end
        maxsz = max(maxsz, sz)
    end
    return maxsz
end

# ----------------------------
# Weighted sampling helpers
# ----------------------------
@inline function weighted_pick_index(rng::AbstractRNG, w::Vector{Float64})
    s = sum(w)
    if s <= 0
        return rand(rng, 1:length(w))
    end
    u = rand(rng) * s
    acc = 0.0
    @inbounds for k in 1:length(w)
        acc += w[k]
        if acc >= u
            return k
        end
    end
    return length(w)
end

function sample_weighted_no_replace(rng::AbstractRNG, cands::Vector{Int}, w::Vector{Float64}, k::Int)
    k = min(k, length(cands))
    chosen = Int[]
    k == 0 && return chosen
    remaining = collect(1:length(cands))
    for _ in 1:k
        ww = w[remaining]
        pick_local = weighted_pick_index(rng, ww)
        pick_idx = remaining[pick_local]
        push!(chosen, cands[pick_idx])
        deleteat!(remaining, pick_local)
        isempty(remaining) && break
    end
    return chosen
end

# ----------------------------
# AUC (trapz) on fgrid
# ----------------------------
function auc_trapz(x::Vector{Float64}, y::Vector{Float64})
    s = 0.0
    @inbounds for i in 1:length(x)-1
        dx = x[i+1] - x[i]
        s += dx * (y[i] + y[i+1]) / 2.0
    end
    return s
end

# ----------------------------
# Environment + abiotic maps
# ----------------------------
function make_environment(rng::AbstractRNG; n1::Int, n2::Int, smooth_iters::Int=30)
    env1 = randn(rng, n1, n2)
    env2 = randn(rng, n1, n2)
    smooth_field!(env1; iters=smooth_iters)
    smooth_field!(env2; iters=smooth_iters)
    zscore!(env1); zscore!(env2)
    return env1, env2
end

function make_abiotic_maps(rng::AbstractRNG, env1::Matrix{Float64}, env2::Matrix{Float64};
                           S::Int, niche_sigma::Float64, niche_cut::Float64)
    n1, n2 = size(env1)
    N = n1*n2
    mu1 = randn(rng, S)
    mu2 = randn(rng, S)
    A = BitMatrix(undef, S, N)
    @inbounds for i in 1:S
        for r in 1:n1, c in 1:n2
            d1 = env1[r,c] - mu1[i]
            d2 = env2[r,c] - mu2[i]
            suit = exp(-(d1*d1 + d2*d2) / (2.0*niche_sigma^2))
            A[i, idx(r,c,n2)] = (suit >= niche_cut)
        end
    end
    return A, mu1, mu2
end

# ----------------------------
# Metaweb (TL + diets)
#   - Consumers pick prey among lower TL species.
#   - match_sigma: consumer↔prey environmental similarity bias
#   - select_sigma: prey guild cohesion around a seed prey
# ----------------------------
function build_metaweb(rng::AbstractRNG;
                       S::Int, basal_frac::Float64, Lmax::Int,
                       k_prey::Int, match_sigma::Float64, select_sigma::Float64,
                       mu1::Vector{Float64}, mu2::Vector{Float64})

    nb = max(1, round(Int, basal_frac*S))
    perm = randperm(rng, S)
    basals = perm[1:nb]
    is_basal = falses(S)
    is_basal[basals] .= true

    TL = fill(1, S)
    @inbounds for i in 1:S
        if !is_basal[i]
            TL[i] = rand(rng, 2:Lmax)
        end
    end

    preylist = [Int[] for _ in 1:S]

    @inbounds for i in 1:S
        TL[i] == 1 && continue

        cands = [j for j in 1:S if TL[j] < TL[i]]
        isempty(cands) && (cands = collect(basals))

        wmatch = Vector{Float64}(undef, length(cands))
        for (t, j) in enumerate(cands)
            d1 = mu1[i]-mu1[j]; d2 = mu2[i]-mu2[j]
            wmatch[t] = exp(-(d1*d1 + d2*d2) / (2.0*match_sigma^2))
        end

        g = cands[weighted_pick_index(rng, wmatch)]
        if k_prey == 1
            preylist[i] = [g]
            continue
        end

        cands2 = [j for j in cands if j != g]
        isempty(cands2) && (preylist[i] = [g]; continue)

        w = Vector{Float64}(undef, length(cands2))
        for (t, j) in enumerate(cands2)
            d1g = mu1[g]-mu1[j]; d2g = mu2[g]-mu2[j]
            wguild = exp(-(d1g*d1g + d2g*d2g) / (2.0*select_sigma^2))
            d1i = mu1[i]-mu1[j]; d2i = mu2[i]-mu2[j]
            wfeas = exp(-(d1i*d1i + d2i*d2i) / (2.0*match_sigma^2))
            w[t] = wguild*wfeas
        end

        others = sample_weighted_no_replace(rng, cands2, w, k_prey-1)
        preylist[i] = vcat([g], others)
    end

    return TL, preylist
end

# ----------------------------
# Habitat loss geometries (nested orders)
# ----------------------------
function make_loss_order(rng::AbstractRNG, geometry::Symbol, env1::Matrix{Float64})
    n1, n2 = size(env1)
    N = n1*n2
    if geometry == :random
        return randperm(rng, N)
    end
    score = Vector{Float64}(undef, N)
    if geometry == :cluster
        Z = randn(rng, n1, n2)
        smooth_field!(Z; iters=35)
        zscore!(Z)
        @inbounds for r in 1:n1, c in 1:n2
            score[idx(r,c,n2)] = Z[r,c]
        end
    elseif geometry == :front
        @inbounds for r in 1:n1, c in 1:n2
            score[idx(r,c,n2)] = env1[r,c]
        end
    else
        error("Unknown geometry: $geometry")
    end
    return sortperm(score; rev=true)
end

function keep_from_order(ord::Vector{Int}, f::Float64, N::Int)
    keep = trues(N)
    nrem = clamp(round(Int, f*N), 0, N)
    @inbounds for t in 1:nrem
        keep[ord[t]] = false
    end
    return keep
end

# ----------------------------
# Core: AB computation under fixed prey rule: ">= 1 prey present in cell"
# Also builds consumer footprint (cells where any consumer is supported)
# ----------------------------
function compute_P_and_counts!(P::BitMatrix,
                              Acount::Vector{Int}, ABcount::Vector{Int},
                              A::BitMatrix, keep::BitVector,
                              TL::Vector{Int}, preylist::Vector{Vector{Int}},
                              cons_present::BitVector, cons_idx::Vector{Int},
                              N::Int)

    S = size(A,1)
    fill!(Acount, 0)
    fill!(ABcount, 0)
    P .= false
    fill!(cons_present, false)
    empty!(cons_idx)

    @inbounds for i in 1:S
        c = 0
        for k in 1:N
            c += (keep[k] & A[i,k]) ? 1 : 0
        end
        Acount[i] = c
    end

    order = sortperm(TL) # bottom-up
    @inbounds for i in order
        if TL[i] == 1
            c = 0
            for k in 1:N
                v = keep[k] & A[i,k]
                P[i,k] = v
                c += v ? 1 : 0
            end
            ABcount[i] = c
        else
            prey = preylist[i]
            c = 0
            for k in 1:N
                sup = false
                for pj in prey
                    sup |= P[pj,k]
                    sup && break
                end
                v = keep[k] & A[i,k] & sup
                P[i,k] = v
                if v
                    c += 1
                    if !cons_present[k]
                        cons_present[k] = true
                        push!(cons_idx, k)
                    end
                end
            end
            ABcount[i] = c
        end
    end

    return nothing
end

# ----------------------------
# Patch sizes per species: max connected component on A-only and AB
# Uses "shortcut": if total count < Epatch, maxpatch=count (no BFS).
# Scratch "present_flag" is cleared only on indices touched.
# ----------------------------
function fill_present_from_A!(present_flag::BitVector, present_idx::Vector{Int},
                             A::BitMatrix, i::Int, keep::BitVector, N::Int)
    empty!(present_idx)
    @inbounds for k in 1:N
        if keep[k] & A[i,k]
            present_flag[k] = true
            push!(present_idx, k)
        end
    end
    return length(present_idx)
end

function fill_present_from_P!(present_flag::BitVector, present_idx::Vector{Int},
                             P::BitMatrix, i::Int, N::Int)
    empty!(present_idx)
    @inbounds for k in 1:N
        if P[i,k]
            present_flag[k] = true
            push!(present_idx, k)
        end
    end
    return length(present_idx)
end

function clear_present!(present_flag::BitVector, present_idx::Vector{Int})
    @inbounds for k in present_idx
        present_flag[k] = false
    end
    return nothing
end

function compute_maxpatches!(Amax::Vector{Int}, ABmax::Vector{Int},
                            A::BitMatrix, P::BitMatrix, keep::BitVector,
                            Acount::Vector{Int}, ABcount::Vector{Int},
                            Epatch::Int,
                            present_flag::BitVector, present_idx::Vector{Int},
                            up::Vector{Int}, down::Vector{Int}, left::Vector{Int}, right::Vector{Int},
                            seen::Vector{Int}, queue::Vector{Int}, stamp::Base.RefValue{Int},
                            N::Int)
    S = length(Acount)
    @inbounds for i in 1:S
        if Acount[i] < Epatch
            Amax[i] = Acount[i]
        else
            fill_present_from_A!(present_flag, present_idx, A, i, keep, N)
            Amax[i] = lcc_max_size(present_flag, present_idx, up, down, left, right, seen, queue, stamp)
            clear_present!(present_flag, present_idx)
        end

        if ABcount[i] < Epatch
            ABmax[i] = ABcount[i]
        else
            fill_present_from_P!(present_flag, present_idx, P, i, N)
            ABmax[i] = lcc_max_size(present_flag, present_idx, up, down, left, right, seen, queue, stamp)
            clear_present!(present_flag, present_idx)
        end
    end
    return nothing
end

# ----------------------------
# SAR fitting at f=0 on A-only, under two presence definitions
#   - total-area presence: count >= Emin_total
#   - patch presence:      maxpatch >= Epatch
# Squares are contiguous; sampled at various area fractions.
# ----------------------------
function sar_fit_total(rng::AbstractRNG, A::BitMatrix;
                       n1::Int, n2::Int, Emin_total::Int,
                       area_fracs::Vector{Float64}, samples::Int)
    S, N = size(A)
    keep_sq = falses(N)
    a_eff = Float64[]
    Sbar  = Float64[]

    for a in area_fracs
        side = max(3, round(Int, sqrt(a) * min(n1, n2)))
        side = min(side, min(n1, n2))
        frac_eff = (side*side) / N

        rich = zeros(Float64, samples)
        for s in 1:samples
            fill!(keep_sq, false)
            r0 = rand(rng, 1:(n1 - side + 1))
            c0 = rand(rng, 1:(n2 - side + 1))
            @inbounds for r in r0:(r0+side-1), c in c0:(c0+side-1)
                keep_sq[idx(r,c,n2)] = true
            end

            cntS = 0
            @inbounds for i in 1:S
                cells = 0
                for k in 1:N
                    cells += (keep_sq[k] & A[i,k]) ? 1 : 0
                end
                cntS += (cells >= Emin_total)
            end
            rich[s] = cntS
        end

        push!(a_eff, frac_eff)
        push!(Sbar, mean(rich))
    end

    xs = Float64[]; ys = Float64[]
    @inbounds for i in eachindex(a_eff)
        if Sbar[i] > 0 && a_eff[i] > 0
            push!(xs, log(a_eff[i]))
            push!(ys, log(Sbar[i]))
        end
    end
    length(xs) < 2 && return (NaN, NaN, NaN)

    x̄ = mean(xs); ȳ = mean(ys)
    z = sum((xs .- x̄) .* (ys .- ȳ)) / sum((xs .- x̄).^2)
    b = ȳ - z*x̄
    c = exp(b)
    yhat = b .+ z .* xs
    ssr = sum((ys .- yhat).^2)
    sst = sum((ys .- ȳ).^2)
    r2 = (sst > 0) ? (1 - ssr/sst) : NaN
    return (c, z, r2)
end

function sar_fit_patch(rng::AbstractRNG, A::BitMatrix;
                       n1::Int, n2::Int, Epatch::Int,
                       area_fracs::Vector{Float64}, samples::Int,
                       up::Vector{Int}, down::Vector{Int}, left::Vector{Int}, right::Vector{Int},
                       present_flag::BitVector, present_idx::Vector{Int},
                       seen::Vector{Int}, queue::Vector{Int}, stamp::Base.RefValue{Int})

    S, N = size(A)
    keep_sq = falses(N)
    a_eff = Float64[]
    Sbar  = Float64[]

    for a in area_fracs
        side = max(3, round(Int, sqrt(a) * min(n1, n2)))
        side = min(side, min(n1, n2))
        frac_eff = (side*side) / N

        rich = zeros(Float64, samples)
        for s in 1:samples
            fill!(keep_sq, false)
            r0 = rand(rng, 1:(n1 - side + 1))
            c0 = rand(rng, 1:(n2 - side + 1))
            @inbounds for r in r0:(r0+side-1), c in c0:(c0+side-1)
                keep_sq[idx(r,c,n2)] = true
            end

            cntS = 0
            @inbounds for i in 1:S
                empty!(present_idx)
                for k in 1:N
                    if keep_sq[k] & A[i,k]
                        present_flag[k] = true
                        push!(present_idx, k)
                    end
                end
                if length(present_idx) >= Epatch
                    maxsz = lcc_max_size(present_flag, present_idx, up, down, left, right, seen, queue, stamp)
                    cntS += (maxsz >= Epatch)
                end
                clear_present!(present_flag, present_idx)
            end
            rich[s] = cntS
        end

        push!(a_eff, frac_eff)
        push!(Sbar, mean(rich))
    end

    xs = Float64[]; ys = Float64[]
    @inbounds for i in eachindex(a_eff)
        if Sbar[i] > 0 && a_eff[i] > 0
            push!(xs, log(a_eff[i]))
            push!(ys, log(Sbar[i]))
        end
    end
    length(xs) < 2 && return (NaN, NaN, NaN)

    x̄ = mean(xs); ȳ = mean(ys)
    z = sum((xs .- x̄) .* (ys .- ȳ)) / sum((xs .- x̄).^2)
    b = ȳ - z*x̄
    c = exp(b)
    yhat = b .+ z .* xs
    ssr = sum((ys .- yhat).^2)
    sst = sum((ys .- ȳ).^2)
    r2 = (sst > 0) ? (1 - ssr/sst) : NaN
    return (c, z, r2)
end

# ----------------------------
# Configuration
# ----------------------------
struct RunConfig
    n1::Int
    n2::Int
    S::Int
    basal_frac::Float64
    Lmax::Int
    niche_sigma::Float64
    niche_cut::Float64
    k_prey::Int
    select_sigma::Float64
    match_sigma::Float64
    Emin_total::Int
    Epatch::Int
    fgrid::Vector{Float64}
    geoms::Vector{Symbol}
    sar_area_fracs::Vector{Float64}
    sar_samples::Int
end

# ----------------------------
# Thread-local scratch buffers (thread-safe)
# ----------------------------
mutable struct ThreadBuffers
    # per-species counts
    Acount::Vector{Int}
    ABcount::Vector{Int}
    Amax::Vector{Int}
    ABmax::Vector{Int}
    # AB presence matrix
    P::BitMatrix
    # BFS scratch
    present_flag::BitVector
    present_idx::Vector{Int}
    cons_present::BitVector
    cons_idx::Vector{Int}
    seen::Vector{Int}
    queue::Vector{Int}
    stamp::Base.RefValue{Int}
end

function make_thread_buffers(cfg::RunConfig)
    N = cfg.n1 * cfg.n2
    return ThreadBuffers(
        zeros(Int, cfg.S),
        zeros(Int, cfg.S),
        zeros(Int, cfg.S),
        zeros(Int, cfg.S),
        falses(cfg.S, N),
        falses(N),
        Int[],
        falses(N),
        Int[],
        zeros(Int, N),
        Vector{Int}(undef, N),
        Ref(0)
    )
end

# ----------------------------
# One replicate: compute curves for each geometry + SAR predictions (rep-specific)
# ----------------------------
function run_one_replicate!(cfg::RunConfig, rng::AbstractRNG,
                            up::Vector{Int}, down::Vector{Int}, left::Vector{Int}, right::Vector{Int},
                            buf::ThreadBuffers)

    n1, n2, S = cfg.n1, cfg.n2, cfg.S
    N = n1*n2
    T = length(cfg.fgrid)
    G = length(cfg.geoms)

    # Build env and abiotic maps + metaweb
    env1, env2 = make_environment(rng; n1=n1, n2=n2, smooth_iters=30)
    A, mu1, mu2 = make_abiotic_maps(rng, env1, env2; S=S, niche_sigma=cfg.niche_sigma, niche_cut=cfg.niche_cut)
    TL, preylist = build_metaweb(rng; S=S, basal_frac=cfg.basal_frac, Lmax=cfg.Lmax,
                                 k_prey=cfg.k_prey, match_sigma=cfg.match_sigma, select_sigma=cfg.select_sigma,
                                 mu1=mu1, mu2=mu2)

    # SAR fits at f=0 on A-only (rep-specific)
    cT, zT, r2T = sar_fit_total(rng, A; n1=n1, n2=n2, Emin_total=cfg.Emin_total,
                               area_fracs=cfg.sar_area_fracs, samples=cfg.sar_samples)

    cP, zP, r2P = sar_fit_patch(rng, A; n1=n1, n2=n2, Epatch=cfg.Epatch,
                               area_fracs=cfg.sar_area_fracs, samples=cfg.sar_samples,
                               up=up, down=down, left=left, right=right,
                               present_flag=buf.present_flag, present_idx=buf.present_idx,
                               seen=buf.seen, queue=buf.queue, stamp=buf.stamp)

    # Precompute nested loss orders per geometry
    ords = [make_loss_order(rng, g, env1) for g in cfg.geoms]

    # Baseline at f=0 (same keep0)
    keep0 = trues(N)
    compute_P_and_counts!(buf.P, buf.Acount, buf.ABcount, A, keep0, TL, preylist, buf.cons_present, buf.cons_idx, N)
    compute_maxpatches!(buf.Amax, buf.ABmax, A, buf.P, keep0, buf.Acount, buf.ABcount, cfg.Epatch,
                        buf.present_flag, buf.present_idx, up, down, left, right, buf.seen, buf.queue, buf.stamp, N)

    S0A_tot  = count(>=(cfg.Emin_total), buf.Acount)
    S0AB_tot = count(>=(cfg.Emin_total), buf.ABcount)
    S0A_pat  = count(>=(cfg.Epatch),     buf.Amax)
    S0AB_pat = count(>=(cfg.Epatch),     buf.ABmax)

    # SAR prediction curves for this replicate (area-only; same for all geometries)
    SAR_total = zeros(Float64, T)
    SAR_patch = zeros(Float64, T)
    @inbounds for (t, f) in enumerate(cfg.fgrid)
        h = max(1e-9, 1.0 - f)
        if isfinite(cT) && isfinite(zT)
            SpT = clamp(cT * h^zT, 0.0, float(S0A_tot))
            SAR_total[t] = S0A_tot - SpT
        else
            SAR_total[t] = NaN
        end
        if isfinite(cP) && isfinite(zP)
            SpP = clamp(cP * h^zP, 0.0, float(S0A_pat))
            SAR_patch[t] = S0A_pat - SpP
        else
            SAR_patch[t] = NaN
        end
    end

    # Outputs per geometry
    EA_tot  = zeros(Float64, G, T)
    EAB_tot = zeros(Float64, G, T)
    EA_pat  = zeros(Float64, G, T)
    EAB_pat = zeros(Float64, G, T)
    phi     = zeros(Float64, G, T)
    lccF    = zeros(Float64, G, T)
    fragfail= zeros(Float64, G, T)

    # Sweep f for each geometry
    for (gj, g) in enumerate(cfg.geoms)
        ord = ords[gj]
        for (t, f) in enumerate(cfg.fgrid)
            keep = (f == 0.0) ? keep0 : keep_from_order(ord, f, N)

            compute_P_and_counts!(buf.P, buf.Acount, buf.ABcount, A, keep, TL, preylist, buf.cons_present, buf.cons_idx, N)
            compute_maxpatches!(buf.Amax, buf.ABmax, A, buf.P, keep, buf.Acount, buf.ABcount, cfg.Epatch,
                                buf.present_flag, buf.present_idx, up, down, left, right, buf.seen, buf.queue, buf.stamp, N)

            SA_tot  = count(>=(cfg.Emin_total), buf.Acount)
            SAB_tot = count(>=(cfg.Emin_total), buf.ABcount)
            SA_pat  = count(>=(cfg.Epatch),     buf.Amax)
            SAB_pat = count(>=(cfg.Epatch),     buf.ABmax)

            EA_tot[gj,t]  = S0A_tot  - SA_tot
            EAB_tot[gj,t] = S0AB_tot - SAB_tot
            EA_pat[gj,t]  = S0A_pat  - SA_pat
            EAB_pat[gj,t] = S0AB_pat - SAB_pat

            # mean φ among consumers
            sφ = 0.0; nφ = 0
            @inbounds for i in 1:S
                if TL[i] > 1 && buf.Acount[i] > 0
                    sφ += buf.ABcount[i] / buf.Acount[i]
                    nφ += 1
                end
            end
            phi[gj,t] = (nφ > 0 ? sφ/nφ : NaN)

            # LCC fraction on consumer-supported overlap footprint
            tot_cons_cells = length(buf.cons_idx)
            if tot_cons_cells == 0
                lccF[gj,t] = 0.0
            else
                maxsz = lcc_max_size(buf.cons_present, buf.cons_idx, up, down, left, right, buf.seen, buf.queue, buf.stamp)
                lccF[gj,t] = maxsz / tot_cons_cells
            end

            # frag-fail rate among consumers: survive total-area but fail patch
            ff = 0; nn = 0
            @inbounds for i in 1:S
                if TL[i] > 1
                    nn += 1
                    if (buf.ABcount[i] >= cfg.Emin_total) && (buf.ABmax[i] < cfg.Epatch)
                        ff += 1
                    end
                end
            end
            fragfail[gj,t] = (nn > 0 ? ff/nn : NaN)
        end
    end

    return (EA_tot=EA_tot, EAB_tot=EAB_tot, EA_pat=EA_pat, EAB_pat=EAB_pat,
            phi=phi, lcc=lccF, fragfail=fragfail,
            SAR_total=SAR_total, SAR_patch=SAR_patch)
end

# ----------------------------
# Multi-replicate (threaded) aggregation: returns mean curves
# ----------------------------
function run_many(cfg::RunConfig; seed::Int=1234, reps::Int=40)
    println("\nRunning reps=$reps with Threads.nthreads()=$(nthreads()) ...")
    n1,n2,S = cfg.n1,cfg.n2,cfg.S
    N = n1*n2
    up, down, left, right = build_neighbors(n1,n2)

    G = length(cfg.geoms)
    T = length(cfg.fgrid)

    # Thread-local buffers
    tbufs = [make_thread_buffers(cfg) for _ in 1:nthreads()]

    # Thread-local accumulators (avoid locks)
    EA_tot_sum  = [zeros(Float64, G, T) for _ in 1:nthreads()]
    EAB_tot_sum = [zeros(Float64, G, T) for _ in 1:nthreads()]
    EA_pat_sum  = [zeros(Float64, G, T) for _ in 1:nthreads()]
    EAB_pat_sum = [zeros(Float64, G, T) for _ in 1:nthreads()]
    phi_sum     = [zeros(Float64, G, T) for _ in 1:nthreads()]
    lcc_sum     = [zeros(Float64, G, T) for _ in 1:nthreads()]
    ff_sum      = [zeros(Float64, G, T) for _ in 1:nthreads()]
    SAR_T_sum   = [zeros(Float64, T) for _ in 1:nthreads()]
    SAR_P_sum   = [zeros(Float64, T) for _ in 1:nthreads()]
    nSAR_T      = zeros(Int, nthreads())
    nSAR_P      = zeros(Int, nthreads())

    @threads for rep in 1:reps
        tid = threadid()
        rng = rep_rng(seed, rep)
        out = run_one_replicate!(cfg, rng, up, down, left, right, tbufs[tid])

        EA_tot_sum[tid]  .+= out.EA_tot
        EAB_tot_sum[tid] .+= out.EAB_tot
        EA_pat_sum[tid]  .+= out.EA_pat
        EAB_pat_sum[tid] .+= out.EAB_pat
        phi_sum[tid]     .+= out.phi
        lcc_sum[tid]     .+= out.lcc
        ff_sum[tid]      .+= out.fragfail

        if all(isfinite, out.SAR_total)
            SAR_T_sum[tid] .+= out.SAR_total
            nSAR_T[tid] += 1
        end
        if all(isfinite, out.SAR_patch)
            SAR_P_sum[tid] .+= out.SAR_patch
            nSAR_P[tid] += 1
        end
    end

    # Reduce thread-local sums
    EA_tot  = zeros(Float64, G, T)
    EAB_tot = zeros(Float64, G, T)
    EA_pat  = zeros(Float64, G, T)
    EAB_pat = zeros(Float64, G, T)
    phi     = zeros(Float64, G, T)
    lcc     = zeros(Float64, G, T)
    ff      = zeros(Float64, G, T)
    SAR_T   = zeros(Float64, T)
    SAR_P   = zeros(Float64, T)

    for tid in 1:nthreads()
        EA_tot  .+= EA_tot_sum[tid]
        EAB_tot .+= EAB_tot_sum[tid]
        EA_pat  .+= EA_pat_sum[tid]
        EAB_pat .+= EAB_pat_sum[tid]
        phi     .+= phi_sum[tid]
        lcc     .+= lcc_sum[tid]
        ff      .+= ff_sum[tid]
        SAR_T   .+= SAR_T_sum[tid]
        SAR_P   .+= SAR_P_sum[tid]
    end

    EA_tot  ./= reps
    EAB_tot ./= reps
    EA_pat  ./= reps
    EAB_pat ./= reps
    phi     ./= reps
    lcc     ./= reps
    ff      ./= reps

    nT = sum(nSAR_T); nP = sum(nSAR_P)
    if nT > 0
        SAR_T ./= nT
    else
        SAR_T .= NaN
    end
    if nP > 0
        SAR_P ./= nP
    else
        SAR_P .= NaN
    end

    return (EA_tot=EA_tot, EAB_tot=EAB_tot, EA_pat=EA_pat, EAB_pat=EAB_pat,
            phi=phi, lcc=lcc, fragfail=ff, SAR_total=SAR_T, SAR_patch=SAR_P)
end

# ----------------------------
# Figure 1 (H1): decomposition A-only vs AB vs SAR (total & patch)
# ----------------------------
function fig_H1(cfg::RunConfig, agg)
    G = length(cfg.geoms)
    T = length(cfg.fgrid)
    fig = Figure(size=(1750, 950))

    for (j,g) in enumerate(cfg.geoms)
        ax = Axis(fig[1,j],
            xlabel="habitat loss f",
            ylabel=(j==1 ? "incremental extinctions" : ""),
            title="TOTAL-AREA criterion (Emin=$(cfg.Emin_total))   geometry=$(String(g))"
        )
        lines!(ax, cfg.fgrid, vec(agg.EA_tot[j,:]);  label="A-only obs", linewidth=3)
        lines!(ax, cfg.fgrid, vec(agg.EAB_tot[j,:]); label="AB obs (biotic)", linewidth=3)
        lines!(ax, cfg.fgrid, agg.SAR_total;         label="SAR baseline (area-only)", linestyle=:dash, linewidth=3)
        axislegend(ax; position=:lt)
    end

    for (j,g) in enumerate(cfg.geoms)
        ax = Axis(fig[2,j],
            xlabel="habitat loss f",
            ylabel=(j==1 ? "incremental extinctions" : ""),
            title="PATCH criterion (Epatch=$(cfg.Epatch))   geometry=$(String(g))"
        )
        lines!(ax, cfg.fgrid, vec(agg.EA_pat[j,:]);  label="A-only obs", linewidth=3)
        lines!(ax, cfg.fgrid, vec(agg.EAB_pat[j,:]); label="AB obs (biotic)", linewidth=3)
        lines!(ax, cfg.fgrid, agg.SAR_patch;         label="SAR baseline (area-only)", linestyle=:dash, linewidth=3)
        axislegend(ax; position=:lt)
    end

    display(fig)
    # return fig
end

# ----------------------------
# Figure 2 (H2): mechanisms + SAR-error vs supported-LCC
# ----------------------------
function fig_H2(cfg::RunConfig, agg)
    G = length(cfg.geoms)
    fig = Figure(size=(1750, 900))

    # Mechanism panels
    for (j,g) in enumerate(cfg.geoms)
        axL = Axis(fig[1,j],
            xlabel="habitat loss f",
            ylabel=(j==1 ? "mean φ (consumers)" : ""),
            title="Mechanisms (geometry=$(String(g)))"
        )
        axR = Axis(fig[1,j],
            yaxisposition=:right,
            ylabel=(j==3 ? "LCC / frag-fail rate" : ""),
            xgridvisible=false, ygridvisible=false
        )
        linkxaxes!(axL, axR)
        hidespines!(axR)
        hidexdecorations!(axR)

        ylims!(axL, (0.0, 1.0))
        ylims!(axR, (0.0, 1.0))

        lines!(axL, cfg.fgrid, vec(agg.phi[j,:]);       linewidth=3, label="mean φ")
        lines!(axR, cfg.fgrid, vec(agg.lcc[j,:]);       linewidth=3, linestyle=:dash, label="supported-LCC fraction")
        lines!(axR, cfg.fgrid, vec(agg.fragfail[j,:]);  linewidth=3, linestyle=:dot,  label="frag-fail rate")

        axislegend(axL; position=:rb)
    end

    # Scatter: SAR error vs LCC for both criteria
    axS = Axis(fig[2,1:3],
        xlabel="supported-LCC fraction (consumer footprint)",
        ylabel="SAR error vs biotic reality (AB_inc − SAR_pred)",
        title="Mechanistic link: fragmentation of trophically supported habitat organizes SAR error"
    )

    for (j,g) in enumerate(cfg.geoms)
        x = vec(agg.lcc[j,:])
        errT = vec(agg.EAB_tot[j,:]) .- agg.SAR_total
        errP = vec(agg.EAB_pat[j,:]) .- agg.SAR_patch

        scatter!(axS, x, errT; markersize=10, label="$(String(g)) total")
        scatter!(axS, x, errP; markersize=10, marker=:utriangle, label="$(String(g)) patch")
        lines!(axS, x, errT; linewidth=2)
        lines!(axS, x, errP; linewidth=2, linestyle=:dash)
    end
    axislegend(axS; position=:rb, nbanks=2)

    display(fig)
    # return fig
end

# ----------------------------
# H3 summaries (threaded) for parameter sweeps:
#   returns per-geometry:
#     - AUC_SARerr_patch  = ∫ (EAB_patch - SAR_patch) df
#     - peak_fragfail     = max_t fragfail(t)
# ----------------------------
function summarize_H3(cfg::RunConfig; seed::Int=999, reps::Int=20)
    agg = run_many(cfg; seed=seed, reps=reps)

    G = length(cfg.geoms)
    auc_err_patch = zeros(Float64, G)
    peak_ff = zeros(Float64, G)

    for j in 1:G
        errP = vec(agg.EAB_pat[j,:]) .- agg.SAR_patch
        auc_err_patch[j] = auc_trapz(cfg.fgrid, errP)
        peak_ff[j] = maximum(vec(agg.fragfail[j,:]))
    end
    return (auc_err_patch=auc_err_patch, peak_ff=peak_ff)
end

# ----------------------------
# Figure 3a (H3): k_prey × select_sigma heatmaps
# ----------------------------
function fig_H3_sweep1(base::RunConfig; seed::Int=2025, reps::Int=18,
                       k_list = [1,2,3,4,6,8],
                       sel_list = [0.15,0.25,0.4,0.6,0.9,1.5])

    G = length(base.geoms)
    nk = length(k_list)
    ns = length(sel_list)

    M_auc = zeros(Float64, G, nk, ns)
    M_ff  = zeros(Float64, G, nk, ns)

    println("\nH3 Sweep1: k_prey × select_sigma (reps=$reps) ...")

    settings = [(ik,is) for ik in 1:nk for is in 1:ns]
    @threads for sidx in 1:length(settings)
        ik, is = settings[sidx]
        cfg = RunConfig(base.n1, base.n2, base.S,
                        base.basal_frac, base.Lmax,
                        base.niche_sigma, base.niche_cut,
                        k_list[ik], sel_list[is], base.match_sigma,
                        base.Emin_total, base.Epatch,
                        base.fgrid, base.geoms,
                        base.sar_area_fracs, base.sar_samples)

        out = summarize_H3(cfg; seed=seed + 1000*ik + is, reps=reps)
        for gj in 1:G
            M_auc[gj,ik,is] = out.auc_err_patch[gj]
            M_ff[gj,ik,is]  = out.peak_ff[gj]
        end
    end

    fig = Figure(size=(1750, 900))
    for (gj,g) in enumerate(base.geoms)
        ax1 = Axis(fig[1,gj], title="AUC(SAR error) PATCH — geometry=$(String(g))",
                   xlabel="select_sigma", ylabel="k_prey")
        hm1 = heatmap!(ax1, 1:ns, 1:nk, M_auc[gj,:,:])
        ax1.xticks = (1:ns, string.(sel_list))
        ax1.yticks = (1:nk, string.(k_list))
        Colorbar(fig[1,gj+3], hm1)

        ax2 = Axis(fig[2,gj], title="peak frag-fail rate — geometry=$(String(g))",
                   xlabel="select_sigma", ylabel="k_prey")
        hm2 = heatmap!(ax2, 1:ns, 1:nk, M_ff[gj,:,:])
        ax2.xticks = (1:ns, string.(sel_list))
        ax2.yticks = (1:nk, string.(k_list))
        Colorbar(fig[2,gj+3], hm2)
    end

    Label(fig[0,1:6],
          "H3 Sweep1: redundancy/synchrony control (fixed rules). Larger AUC means SAR underestimates AB more.",
          fontsize=16)

    display(fig)
    # return fig
end

# ----------------------------
# Figure 3b (H3): Lmax × basal_frac heatmaps
# ----------------------------
function fig_H3_sweep2(base::RunConfig; seed::Int=3030, reps::Int=18,
                       L_list = [3,4,5,6],
                       basal_list = [0.10,0.15,0.20,0.25,0.30,0.35])

    G = length(base.geoms)
    nL = length(L_list)
    nb = length(basal_list)

    M_auc = zeros(Float64, G, nL, nb)
    M_ff  = zeros(Float64, G, nL, nb)

    println("\nH3 Sweep2: Lmax × basal_frac (reps=$reps) ...")

    settings = [(iL,ib) for iL in 1:nL for ib in 1:nb]
    @threads for sidx in 1:length(settings)
        iL, ib = settings[sidx]
        cfg = RunConfig(base.n1, base.n2, base.S,
                        basal_list[ib], L_list[iL],
                        base.niche_sigma, base.niche_cut,
                        base.k_prey, base.select_sigma, base.match_sigma,
                        base.Emin_total, base.Epatch,
                        base.fgrid, base.geoms,
                        base.sar_area_fracs, base.sar_samples)

        out = summarize_H3(cfg; seed=seed + 1000*iL + ib, reps=reps)
        for gj in 1:G
            M_auc[gj,iL,ib] = out.auc_err_patch[gj]
            M_ff[gj,iL,ib]  = out.peak_ff[gj]
        end
    end

    fig = Figure(size=(1750, 900))
    for (gj,g) in enumerate(base.geoms)
        ax1 = Axis(fig[1,gj], title="AUC(SAR error) PATCH — geometry=$(String(g))",
                   xlabel="basal_frac", ylabel="Lmax")
        hm1 = heatmap!(ax1, 1:nb, 1:nL, M_auc[gj,:,:])
        ax1.xticks = (1:nb, string.(basal_list))
        ax1.yticks = (1:nL, string.(L_list))
        Colorbar(fig[1,gj+3], hm1)

        ax2 = Axis(fig[2,gj], title="peak frag-fail rate — geometry=$(String(g))",
                   xlabel="basal_frac", ylabel="Lmax")
        hm2 = heatmap!(ax2, 1:nb, 1:nL, M_ff[gj,:,:])
        ax2.xticks = (1:nb, string.(basal_list))
        ax2.yticks = (1:nL, string.(L_list))
        Colorbar(fig[2,gj+3], hm2)
    end

    Label(fig[0,1:6],
          "H3 Sweep2: trophic architecture control (fixed rules).",
          fontsize=16)

    display(fig)
end

# ----------------------------
# MAIN: choose a larger-than-MVP parameterisation (more species)
# ----------------------------
function main()
    # --- Larger system for "definite" results (adjust as needed) ---
    # Keep N moderate so patch BFS stays feasible.
    n1, n2 = 120, 120          # 14,400 cells
    S      = 250               # more species than MVP
    reps_H12 = 45              # baseline robustness
    reps_H3  = 18              # per-setting reps for sweeps (already expensive)

    cfg = RunConfig(
        n1, n2, S,
        0.25,      # basal_frac
        5,         # Lmax
        0.85,      # niche_sigma
        0.45,      # niche_cut
        3,         # k_prey
        0.60,      # select_sigma
        1.00,      # match_sigma
        80,        # Emin_total (scale up with grid; 80 ~ small patch on 120x120)
        80,        # Epatch
        collect(0.0:0.05:0.90),
        [:random, :cluster, :front],
        [0.20, 0.30, 0.40, 0.55, 0.70, 0.85, 1.0],
        20         # SAR samples per area fraction
    )

    println("===============================================")
    println("STATIC FULL PIPELINE (H1–H3), THREADED + SAFE")
    println("grid=$(cfg.n1)x$(cfg.n2) (N=$(cfg.n1*cfg.n2)), S=$(cfg.S)")
    println("threads=$(nthreads())")
    println("===============================================")

    # ----------------------------
    # H1 + H2: baseline aggregation + Figures 1 and 2
    # ----------------------------
    println("\n=== H1/H2 baseline run ===")
    agg = run_many(cfg; seed=1234, reps=reps_H12)

    println("\nFigure 1 (H1): A-only vs AB vs SAR, total-area and patch criteria")
    fig_H1(cfg, agg)

    println("\nFigure 2 (H2): mechanisms + SAR-error vs supported-LCC")
    fig_H2(cfg, agg)

    # ----------------------------
    # H3: sweeps (bigger, but still feasible)
    # ----------------------------
    println("\n=== H3 sweeps (fixed rules; vary system properties) ===")

    println("\nFigure 3a (H3): sweep k_prey × select_sigma")
    fig_H3_sweep1(cfg; seed=777, reps=reps_H3)

    println("\nFigure 3b (H3): sweep Lmax × basal_frac")
    fig_H3_sweep2(cfg; seed=888, reps=reps_H3)

    return nothing
end

main()