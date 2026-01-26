###############################################################
# suitHab_bioletters_static_vNext_FIXED.jl
#
# FIXED version of your "second" script:
#   - Restores the *working mechanics* from the first script:
#       (1) 2D smoothed environmental fields + Gaussian niche threshold (BitMatrix A[S,N])
#       (2) trophic levels + diets with (match_sigma, select_sigma) and a guild-centre prey
#       (3) habitat loss as *nested orders* (monotone in f) for all geometries
#           • random: permutation
#           • cluster: smoothed-noise ranking
#           • front: env-biased ranking (like your first)
#       (4) fast LCC (stamp BFS) and per-species patch BFS only when needed
#       (5) SAR baseline fit from random sub-squares at f=0 (geometry-independent)
#       (6) SAR_eff fit as a mechanistic bridge: AB richness vs effective supported area (and its LCC)
#
# What it produces (DISPLAY ONLY; no saving):
#   Figure 1: Richness vs habitat loss (TOTAL + PATCH) for 3 geometries,
#             with A-only, AB, SAR (A fit), SAR_eff (support-aware) + mechanism row.
#   Figure 2: Bridge scatter: prediction error vs supported-LCC fraction (raw SAR vs SAR_eff),
#             shown for TOTAL and PATCH.
#   Figure 3: Patch-connectivity sensitivity (r=0 vs r=2) for random/cluster (sanity check).
#   Figure 4: Heatmaps (PARAM axes: k_prey × select_sigma) using a LIGHT runner
#             (no per-species patch BFS) for smooth surfaces.
#   Figure 5: Heatmaps (REALISED axes: redundancy × synchrony proxies), binned.
#
# Run:
#   JULIA_NUM_THREADS=12 julia suitHab_bioletters_static_vNext_FIXED.jl
###############################################################

using Random
using Statistics
using Printf
using CairoMakie
using Base.Threads

# ============================================================
# USER KNOBS
# ============================================================

const DEFAULT_SEED = 1234

# --- Main (full) run ---
const N1_MAIN   = 120
const N2_MAIN   = 120
const S_MAIN    = 250
const REPS_MAIN = 20

const FGRID_MAIN = collect(0.0:0.05:0.9)

# survival thresholds (cells)
const EMIN_MAIN   = 1
const EPATCH_MAIN = 10

# niche + metaweb
const BASAL_FRAC = 0.25
const LMAX       = 5
const NICHE_SIGMA = 0.3
const NICHE_CUT   = 0.05
const MATCH_SIGMA = 0.5
const KPREY_MAIN  = 3
const SELECT_MAIN = 0.60

const GEOMS = [:random, :cluster, :front]

# SAR sampling (sub-squares) for baseline fits at f=0
const SAR_AREAFRACS = [0.1, 0.20, 0.30, 0.40, 0.5, 0.6, 0.70, 0.8, 0.9, 1.0]
const SAR_SAMPLES   = 20

# --- Movement sanity (patch connectivity radius) ---
const DO_MOVEMENT_SANITY = true
const R_CONN_SENS   = 2
const REPS_SENS     = 10
const GEOMS_SENS    = [:random, :cluster]

# --- Sweep (LIGHT runner: no per-species patch BFS) ---
const DO_SWEEP      = true
const N1_SWEEP      = 100
const N2_SWEEP      = 100
const S_SWEEP       = 180
const REPS_SWEEP    = 25
const FGRID_SWEEP   = collect(0.0:0.05:0.90)

const KPREY_LIST    = collect(1:12)
const SELECT_LIST   = [0.10, 0.14, 0.20, 0.28, 0.40, 0.55, 0.75, 1.00, 1.35, 1.80]

# LIGHT thresholds scaled roughly with area
const EMIN_SWEEP    = 70
const EPATCH_SWEEP  = 70

# realised-axis binning
const NBX = 10
const NBY = 10

# supported-LCC collapse threshold (optional summary)
const LCC_THRESH = 0.5

# ============================================================
# RNG utilities (thread-schedule independent)
# ============================================================

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

@inline idx(i::Int, j::Int, n2::Int) = (i - 1) * n2 + j

# ============================================================
# Grid helpers
# ============================================================

function make_rowcol(n1::Int, n2::Int)
    N = n1*n2
    row = Vector{Int}(undef, N)
    col = Vector{Int}(undef, N)
    @inbounds for k in 1:N
        r = (k - 1) ÷ n2 + 1
        c = (k - 1) % n2 + 1
        row[k] = r
        col[k] = c
    end
    return row, col
end

function build_neighbors4(n1::Int, n2::Int)
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

# Chebyshev radius (square neighborhood): interpretable "gap-crossing" radius
function radius_offsets(r::Int)
    r <= 0 && return Tuple{Int,Int}[]
    offs = Tuple{Int,Int}[]
    for dr in -r:r, dc in -r:r
        (dr == 0 && dc == 0) && continue
        if max(abs(dr), abs(dc)) <= r
            push!(offs, (dr, dc))
        end
    end
    return offs
end

# ============================================================
# Fast LCC (stamp BFS) on flattened grid
# ============================================================

function lcc_max_size4(present_flag::BitVector, present_idx::Vector{Int},
                       up::Vector{Int}, down::Vector{Int}, left::Vector{Int}, right::Vector{Int},
                       seen::Vector{Int}, queue::Vector{Int}, stamp::Base.RefValue{Int})

    isempty(present_idx) && return 0
    stamp[] += 1
    s = stamp[]
    if s == typemax(Int)
        fill!(seen, 0); stamp[] = 1; s = 1
    end

    maxsz = 0
    @inbounds for start in present_idx
        if seen[start] == s
            continue
        end
        head = 1; tail = 1
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

function lcc_max_size_radius(present_flag::BitVector, present_idx::Vector{Int},
                             row::Vector{Int}, col::Vector{Int}, n1::Int, n2::Int, n2stride::Int,
                             offs::Vector{Tuple{Int,Int}},
                             seen::Vector{Int}, queue::Vector{Int}, stamp::Base.RefValue{Int})

    isempty(present_idx) && return 0
    stamp[] += 1
    s = stamp[]
    if s == typemax(Int)
        fill!(seen, 0); stamp[] = 1; s = 1
    end

    maxsz = 0
    @inbounds for start in present_idx
        if seen[start] == s
            continue
        end
        head = 1; tail = 1
        queue[1] = start
        seen[start] = s
        sz = 0
        while head <= tail
            v = queue[head]; head += 1
            sz += 1
            rv = row[v]; cv = col[v]
            for (dr, dc) in offs
                rr = rv + dr; cc = cv + dc
                if 1 <= rr <= n1 && 1 <= cc <= n2
                    nb = (rr - 1) * n2stride + cc
                    if present_flag[nb] && seen[nb] != s
                        tail += 1
                        queue[tail] = nb
                        seen[nb] = s
                    end
                end
            end
        end
        maxsz = max(maxsz, sz)
    end
    return maxsz
end

# ============================================================
# Environment smoothing + z-score
# ============================================================

function smooth_field!(Z::Matrix{Float64}; iters::Int=30)
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

function make_environment(rng::AbstractRNG; n1::Int, n2::Int, smooth_iters::Int=30)
    env1 = randn(rng, n1, n2)
    env2 = randn(rng, n1, n2)
    smooth_field!(env1; iters=smooth_iters)
    smooth_field!(env2; iters=smooth_iters)
    zscore!(env1); zscore!(env2)
    return env1, env2
end

# ============================================================
# Abiotic maps A[S,N] from 2D Gaussian niche + threshold
# ============================================================

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

# ============================================================
# Trophic levels + diets (match + guild cohesion)
# ============================================================

function make_trophic_levels(rng::AbstractRNG; S::Int, basal_frac::Float64, Lmax::Int)
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
    return TL, basals, is_basal
end

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

function build_diets(rng::AbstractRNG;
                     S::Int, TL::Vector{Int}, basals::Vector{Int},
                     k_prey::Int, match_sigma::Float64, select_sigma::Float64,
                     mu1::Vector{Float64}, mu2::Vector{Float64})

    preylist = [Int[] for _ in 1:S]

    @inbounds for i in 1:S
        TL[i] == 1 && continue
        cands = [j for j in 1:S if TL[j] < TL[i]]
        isempty(cands) && (cands = copy(basals))

        # (1) pick a guild-centre prey based on consumer-prey match
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

        # (2) fill remaining prey with guild cohesion around g + feasibility to consumer
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

    return preylist
end

# ============================================================
# Habitat loss geometries (nested orders)
# ============================================================

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
        # env-biased: remove highest env1 first (like your first script)
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

# ============================================================
# AB fixed-point support map P[S,N] + counts + consumer footprint
# ============================================================

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

    # A-only counts
    @inbounds for i in 1:S
        c = 0
        for k in 1:N
            c += (keep[k] & A[i,k]) ? 1 : 0
        end
        Acount[i] = c
    end

    # AB support bottom-up by TL
    order = sortperm(TL)
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

# ============================================================
# Patch (max connected component) per species
# ============================================================

function clear_present!(present_flag::BitVector, present_idx::Vector{Int})
    @inbounds for k in present_idx
        present_flag[k] = false
    end
    return nothing
end

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

function compute_maxpatches!(Amax::Vector{Int}, ABmax::Vector{Int},
                            A::BitMatrix, P::BitMatrix, keep::BitVector,
                            Acount::Vector{Int}, ABcount::Vector{Int},
                            Epatch::Int, rconn::Int,
                            present_flag::BitVector, present_idx::Vector{Int},
                            up::Vector{Int}, down::Vector{Int}, left::Vector{Int}, right::Vector{Int},
                            row::Vector{Int}, col::Vector{Int}, n1::Int, n2::Int, offs::Vector{Tuple{Int,Int}},
                            seen::Vector{Int}, queue::Vector{Int}, stamp::Base.RefValue{Int},
                            N::Int)

    S = length(Acount)

    @inbounds for i in 1:S
        # A-only
        if Acount[i] < Epatch
            Amax[i] = Acount[i]
        else
            fill_present_from_A!(present_flag, present_idx, A, i, keep, N)
            if rconn == 0
                Amax[i] = lcc_max_size4(present_flag, present_idx, up, down, left, right, seen, queue, stamp)
            else
                Amax[i] = lcc_max_size_radius(present_flag, present_idx, row, col, n1, n2, n2, offs, seen, queue, stamp)
            end
            clear_present!(present_flag, present_idx)
        end

        # AB
        if ABcount[i] < Epatch
            ABmax[i] = ABcount[i]
        else
            fill_present_from_P!(present_flag, present_idx, P, i, N)
            if rconn == 0
                ABmax[i] = lcc_max_size4(present_flag, present_idx, up, down, left, right, seen, queue, stamp)
            else
                ABmax[i] = lcc_max_size_radius(present_flag, present_idx, row, col, n1, n2, n2, offs, seen, queue, stamp)
            end
            clear_present!(present_flag, present_idx)
        end
    end
    return nothing
end

# ============================================================
# SAR power fit helpers
# ============================================================

function sar_fit_power(xs::Vector{Float64}, ys::Vector{Float64})
    x = Float64[]; y = Float64[]
    @inbounds for i in eachindex(xs)
        if xs[i] > 0 && ys[i] > 0 && isfinite(xs[i]) && isfinite(ys[i])
            push!(x, log(xs[i]))
            push!(y, log(ys[i]))
        end
    end
    length(x) < 2 && return (NaN, NaN, NaN)
    x̄ = mean(x); ȳ = mean(y)
    denom = sum((x .- x̄).^2)
    z = denom > 0 ? sum((x .- x̄) .* (y .- ȳ)) / denom : NaN
    b = ȳ - z*x̄
    c = exp(b)
    yhat = b .+ z .* x
    ssr = sum((y .- yhat).^2)
    sst = sum((y .- ȳ).^2)
    r2 = (sst > 0) ? (1 - ssr/sst) : NaN
    return (c, z, r2)
end

# ============================================================
# SAR fits at f=0 using random sub-squares (geometry-independent)
# ============================================================

function sar_fit_total_Aonly(rng::AbstractRNG, A::BitMatrix; n1::Int, n2::Int, Emin_total::Int,
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
        push!(a_eff, frac_eff); push!(Sbar, mean(rich))
    end
    return sar_fit_power(a_eff, Sbar)
end

function sar_fit_patch_Aonly(rng::AbstractRNG, A::BitMatrix; n1::Int, n2::Int, Epatch::Int,
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
                    maxsz = lcc_max_size4(present_flag, present_idx, up, down, left, right, seen, queue, stamp)
                    cntS += (maxsz >= Epatch)
                end
                clear_present!(present_flag, present_idx)
            end
            rich[s] = cntS
        end
        push!(a_eff, frac_eff); push!(Sbar, mean(rich))
    end
    return sar_fit_power(a_eff, Sbar)
end

# ============================================================
# SAR_eff bridge fit at f=0: AB richness vs effective supported area / supported-LCC area
# ============================================================

function sar_fit_eff_AB(rng::AbstractRNG, A::BitMatrix, TL::Vector{Int}, preylist::Vector{Vector{Int}};  # uses A at f=0
                        n1::Int, n2::Int, Emin_total::Int, Epatch::Int,
                        area_fracs::Vector{Float64}, samples::Int,
                        up::Vector{Int}, down::Vector{Int}, left::Vector{Int}, right::Vector{Int},
                        row::Vector{Int}, col::Vector{Int},
                        P::BitMatrix,
                        Acount::Vector{Int}, ABcount::Vector{Int},
                        Amax::Vector{Int}, ABmax::Vector{Int},
                        present_flag::BitVector, present_idx::Vector{Int},
                        cons_present::BitVector, cons_idx::Vector{Int},
                        seen::Vector{Int}, queue::Vector{Int}, stamp::Base.RefValue{Int})

    S, N = size(A)
    keep_sq = falses(N)

    x_eff_tot = Float64[]; y_S_tot = Float64[]
    x_eff_lcc = Float64[]; y_S_pat = Float64[]

    for a in area_fracs
        side = max(3, round(Int, sqrt(a) * min(n1, n2)))
        side = min(side, min(n1, n2))
        for _ in 1:samples
            fill!(keep_sq, false)
            r0 = rand(rng, 1:(n1 - side + 1))
            c0 = rand(rng, 1:(n2 - side + 1))
            @inbounds for r in r0:(r0+side-1), c in c0:(c0+side-1)
                keep_sq[idx(r,c,n2)] = true
            end

            compute_P_and_counts!(P, Acount, ABcount, A, keep_sq, TL, preylist, cons_present, cons_idx, N)

            S_AB_tot = count(>=(Emin_total), ABcount)

            footprint = length(cons_idx)
            if footprint > 0
                push!(x_eff_tot, footprint / N)
                push!(y_S_tot, S_AB_tot)
            end

            compute_maxpatches!(Amax, ABmax, A, P, keep_sq, Acount, ABcount, Epatch, 0,
                                present_flag, present_idx, up, down, left, right,
                                row, col, n1, n2, Tuple{Int,Int}[],
                                seen, queue, stamp, N)
            S_AB_pat = count(>=(Epatch), ABmax)

            if footprint > 0
                maxsz = lcc_max_size4(cons_present, cons_idx, up, down, left, right, seen, queue, stamp)
                push!(x_eff_lcc, maxsz / N)
                push!(y_S_pat, S_AB_pat)
            end
        end
    end

    cT, zT, r2T = sar_fit_power(x_eff_tot, y_S_tot)
    cP, zP, r2P = sar_fit_power(x_eff_lcc, y_S_pat)
    return (cEffTot=cT, zEffTot=zT, r2EffTot=r2T,
            cEffLCC=cP, zEffLCC=zP, r2EffLCC=r2P)
end

# ============================================================
# realised proxies at f=0 (redundancy × synchrony)
# ============================================================

function realised_proxies!(A::BitMatrix, TL::Vector{Int}, preylist::Vector{Vector{Int}},
                           P::BitMatrix, Acount::Vector{Int}, ABcount::Vector{Int})
    S, N = size(A)
    φs = Float64[]
    reds = Float64[]
    @inbounds for i in 1:S
        TL[i] == 1 && continue
        ai = Acount[i]
        ai <= 0 && continue
        asup = ABcount[i]
        push!(φs, asup / ai)

        # redundancy proxy: average prey overlap count per supported cell
        if asup > 0
            so = 0
            for j in preylist[i]
                c = 0
                for k in 1:N
                    c += (A[i,k] & A[j,k]) ? 1 : 0
                end
                so += c
            end
            push!(reds, so / asup)
        else
            push!(reds, 0.0)
        end
    end
    synchrony = isempty(φs) ? 0.0 : mean(φs)
    redundancy = isempty(reds) ? 0.0 : mean(reds)
    return redundancy, synchrony
end

# ============================================================
# Config + buffers
# ============================================================

struct RunConfiguration
    n1::Int
    n2::Int
    S::Int
    basal_frac::Float64
    Lmax::Int
    niche_sigma::Float64
    niche_cut::Float64
    match_sigma::Float64
    k_prey::Int
    select_sigma::Float64
    Emin_total::Int
    Epatch::Int
    fgrid::Vector{Float64}
    geoms::Vector{Symbol}
    sar_area_fracs::Vector{Float64}
    sar_samples::Int
end

mutable struct ThreadBuffers
    Acount::Vector{Int}
    ABcount::Vector{Int}
    Amax::Vector{Int}
    ABmax::Vector{Int}
    P::BitMatrix
    present_flag::BitVector
    present_idx::Vector{Int}
    cons_present::BitVector
    cons_idx::Vector{Int}
    seen::Vector{Int}
    queue::Vector{Int}
    stamp::Base.RefValue{Int}
end

function make_thread_buffers(cfg::RunConfiguration)
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

# ============================================================
# One replicate (FULL): main figures use per-species patch BFS (r=0)
# ============================================================

function run_one_replicate!(cfg::RunConfiguration, rng::AbstractRNG,
                            up::Vector{Int}, down::Vector{Int}, left::Vector{Int}, right::Vector{Int},
                            row::Vector{Int}, col::Vector{Int},
                            buf::ThreadBuffers)

    n1, n2, S = cfg.n1, cfg.n2, cfg.S
    N = n1*n2
    T = length(cfg.fgrid)
    G = length(cfg.geoms)

    env1, env2 = make_environment(rng; n1=n1, n2=n2, smooth_iters=30)
    A, mu1, mu2 = make_abiotic_maps(rng, env1, env2; S=S, niche_sigma=cfg.niche_sigma, niche_cut=cfg.niche_cut)

    TL, basals, _isbasal = make_trophic_levels(rng; S=S, basal_frac=cfg.basal_frac, Lmax=cfg.Lmax)
    preylist = build_diets(rng; S=S, TL=TL, basals=basals,
                           k_prey=cfg.k_prey, match_sigma=cfg.match_sigma, select_sigma=cfg.select_sigma,
                           mu1=mu1, mu2=mu2)

    ords = [make_loss_order(rng, g, env1) for g in cfg.geoms]

    keep0 = trues(N)
    compute_P_and_counts!(buf.P, buf.Acount, buf.ABcount, A, keep0, TL, preylist, buf.cons_present, buf.cons_idx, N)
    compute_maxpatches!(buf.Amax, buf.ABmax, A, buf.P, keep0,
                        buf.Acount, buf.ABcount, cfg.Epatch, 0,
                        buf.present_flag, buf.present_idx, up, down, left, right,
                        row, col, n1, n2, Tuple{Int,Int}[],
                        buf.seen, buf.queue, buf.stamp, N)

    red0, syn0 = realised_proxies!(A, TL, preylist, buf.P, buf.Acount, buf.ABcount)

    S0A_tot  = count(>=(cfg.Emin_total), buf.Acount)
    S0AB_tot = count(>=(cfg.Emin_total), buf.ABcount)
    S0A_pat  = count(>=(cfg.Epatch),     buf.Amax)
    S0AB_pat = count(>=(cfg.Epatch),     buf.ABmax)

    cT, zT, r2T = sar_fit_total_Aonly(rng, A; n1=n1, n2=n2, Emin_total=cfg.Emin_total,
                                      area_fracs=cfg.sar_area_fracs, samples=cfg.sar_samples)

    cP, zP, r2P = sar_fit_patch_Aonly(rng, A; n1=n1, n2=n2, Epatch=cfg.Epatch,
                                      area_fracs=cfg.sar_area_fracs, samples=cfg.sar_samples,
                                      up=up, down=down, left=left, right=right,
                                      present_flag=buf.present_flag, present_idx=buf.present_idx,
                                      seen=buf.seen, queue=buf.queue, stamp=buf.stamp)

    eff = sar_fit_eff_AB(rng, A, TL, preylist;
                         n1=n1, n2=n2, Emin_total=cfg.Emin_total, Epatch=cfg.Epatch,
                         area_fracs=cfg.sar_area_fracs, samples=max(10, cfg.sar_samples ÷ 2),
                         up=up, down=down, left=left, right=right,
                         row=row, col=col,
                         P=buf.P,
                         Acount=buf.Acount, ABcount=buf.ABcount,
                         Amax=buf.Amax, ABmax=buf.ABmax,
                         present_flag=buf.present_flag, present_idx=buf.present_idx,
                         cons_present=buf.cons_present, cons_idx=buf.cons_idx,
                         seen=buf.seen, queue=buf.queue, stamp=buf.stamp)

    SA_tot   = zeros(Float64, G, T)
    SAB_tot  = zeros(Float64, G, T)
    SA_pat   = zeros(Float64, G, T)
    SAB_pat  = zeros(Float64, G, T)

    Aeff_tot = zeros(Float64, G, T)
    Aeff_lcc = zeros(Float64, G, T)

    phi      = zeros(Float64, G, T)
    lccFrac  = zeros(Float64, G, T)
    fragfail = zeros(Float64, G, T)

    for (gj, g) in enumerate(cfg.geoms)
        ord = ords[gj]
        for (t, f) in enumerate(cfg.fgrid)
            keep = (f == 0.0) ? keep0 : keep_from_order(ord, f, N)

            compute_P_and_counts!(buf.P, buf.Acount, buf.ABcount, A, keep, TL, preylist, buf.cons_present, buf.cons_idx, N)
            compute_maxpatches!(buf.Amax, buf.ABmax, A, buf.P, keep,
                                buf.Acount, buf.ABcount, cfg.Epatch, 0,
                                buf.present_flag, buf.present_idx, up, down, left, right,
                                row, col, n1, n2, Tuple{Int,Int}[],
                                buf.seen, buf.queue, buf.stamp, N)

            SA_tot[gj,t]  = count(>=(cfg.Emin_total), buf.Acount)
            SAB_tot[gj,t] = count(>=(cfg.Emin_total), buf.ABcount)
            SA_pat[gj,t]  = count(>=(cfg.Epatch),     buf.Amax)
            SAB_pat[gj,t] = count(>=(cfg.Epatch),     buf.ABmax)

            # mean φ among consumers
            sφ = 0.0; nφ = 0
            @inbounds for i in 1:S
                if TL[i] > 1 && buf.Acount[i] > 0
                    sφ += buf.ABcount[i] / buf.Acount[i]
                    nφ += 1
                end
            end
            phi[gj,t] = (nφ > 0 ? sφ/nφ : NaN)

            # consumer footprint + its LCC
            footprint = length(buf.cons_idx)
            Aeff_tot[gj,t] = footprint / N
            if footprint == 0
                lccFrac[gj,t] = 0.0
                Aeff_lcc[gj,t] = 0.0
            else
                maxsz = lcc_max_size4(buf.cons_present, buf.cons_idx, up, down, left, right, buf.seen, buf.queue, buf.stamp)
                lccFrac[gj,t] = maxsz / footprint
                Aeff_lcc[gj,t] = maxsz / N
            end

            # frag-fail rate among consumers: total survives but patch fails
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

    fits = (cT=cT, zT=zT, r2T=r2T, cP=cP, zP=zP, r2P=r2P,
            cEffTot=eff.cEffTot, zEffTot=eff.zEffTot, r2EffTot=eff.r2EffTot,
            cEffLCC=eff.cEffLCC, zEffLCC=eff.zEffLCC, r2EffLCC=eff.r2EffLCC)

    base = (S0A_tot=S0A_tot, S0AB_tot=S0AB_tot, S0A_pat=S0A_pat, S0AB_pat=S0AB_pat)

    return (SA_tot=SA_tot, SAB_tot=SAB_tot, SA_pat=SA_pat, SAB_pat=SAB_pat,
            Aeff_tot=Aeff_tot, Aeff_lcc=Aeff_lcc,
            phi=phi, lcc=lccFrac, fragfail=fragfail,
            fits=fits, base=base,
            red0=red0, syn0=syn0)
end

# ============================================================
# Aggregation (FULL)
# ============================================================

function run_many(cfg::RunConfiguration; seed::Int, reps::Int)
    println("\nRunning FULL reps=$reps (threads=$(nthreads())) ...")
    n1,n2 = cfg.n1, cfg.n2
    up, down, left, right = build_neighbors4(n1,n2)
    row, col = make_rowcol(n1,n2)

    G = length(cfg.geoms)
    T = length(cfg.fgrid)

    tbufs = [make_thread_buffers(cfg) for _ in 1:nthreads()]

    SA_tot_sum   = [zeros(Float64, G, T) for _ in 1:nthreads()]
    SAB_tot_sum  = [zeros(Float64, G, T) for _ in 1:nthreads()]
    SA_pat_sum   = [zeros(Float64, G, T) for _ in 1:nthreads()]
    SAB_pat_sum  = [zeros(Float64, G, T) for _ in 1:nthreads()]
    Aeff_tot_sum = [zeros(Float64, G, T) for _ in 1:nthreads()]
    Aeff_lcc_sum = [zeros(Float64, G, T) for _ in 1:nthreads()]
    phi_sum      = [zeros(Float64, G, T) for _ in 1:nthreads()]
    lcc_sum      = [zeros(Float64, G, T) for _ in 1:nthreads()]
    ff_sum       = [zeros(Float64, G, T) for _ in 1:nthreads()]

    fits_acc = [zeros(Float64, 12) for _ in 1:nthreads()] # cT zT cP zP cEffTot zEffTot cEffLCC zEffLCC r2T r2P r2EffTot r2EffLCC
    base_acc = [zeros(Float64, 4)  for _ in 1:nthreads()]
    red_acc  = zeros(Float64, nthreads())
    syn_acc  = zeros(Float64, nthreads())
    nfits    = zeros(Int, nthreads())

    @threads for rep in 1:reps
        tid = threadid()
        rng = rep_rng(seed, rep)
        out = run_one_replicate!(cfg, rng, up, down, left, right, row, col, tbufs[tid])

        SA_tot_sum[tid]   .+= out.SA_tot
        SAB_tot_sum[tid]  .+= out.SAB_tot
        SA_pat_sum[tid]   .+= out.SA_pat
        SAB_pat_sum[tid]  .+= out.SAB_pat
        Aeff_tot_sum[tid] .+= out.Aeff_tot
        Aeff_lcc_sum[tid] .+= out.Aeff_lcc
        phi_sum[tid]      .+= out.phi
        lcc_sum[tid]      .+= out.lcc
        ff_sum[tid]       .+= out.fragfail

        red_acc[tid] += out.red0
        syn_acc[tid] += out.syn0

        f = out.fits
        if all(isfinite, (f.cT,f.zT,f.cP,f.zP,f.cEffTot,f.zEffTot,f.cEffLCC,f.zEffLCC))
            v = fits_acc[tid]
            v[1]  += f.cT;       v[2]  += f.zT
            v[3]  += f.cP;       v[4]  += f.zP
            v[5]  += f.cEffTot;  v[6]  += f.zEffTot
            v[7]  += f.cEffLCC;  v[8]  += f.zEffLCC
            v[9]  += f.r2T;      v[10] += f.r2P
            v[11] += f.r2EffTot; v[12] += f.r2EffLCC

            b = base_acc[tid]
            b[1] += out.base.S0A_tot
            b[2] += out.base.S0AB_tot
            b[3] += out.base.S0A_pat
            b[4] += out.base.S0AB_pat

            nfits[tid] += 1
        end
    end

    SA_tot   = zeros(Float64, G, T)
    SAB_tot  = zeros(Float64, G, T)
    SA_pat   = zeros(Float64, G, T)
    SAB_pat  = zeros(Float64, G, T)
    Aeff_tot = zeros(Float64, G, T)
    Aeff_lcc = zeros(Float64, G, T)
    phi      = zeros(Float64, G, T)
    lcc      = zeros(Float64, G, T)
    ff       = zeros(Float64, G, T)

    fits = zeros(Float64, 12)
    base = zeros(Float64, 4)
    nf   = sum(nfits)

    red0 = 0.0
    syn0 = 0.0

    for tid in 1:nthreads()
        SA_tot   .+= SA_tot_sum[tid]
        SAB_tot  .+= SAB_tot_sum[tid]
        SA_pat   .+= SA_pat_sum[tid]
        SAB_pat  .+= SAB_pat_sum[tid]
        Aeff_tot .+= Aeff_tot_sum[tid]
        Aeff_lcc .+= Aeff_lcc_sum[tid]
        phi      .+= phi_sum[tid]
        lcc      .+= lcc_sum[tid]
        ff       .+= ff_sum[tid]
        fits     .+= fits_acc[tid]
        base     .+= base_acc[tid]
        red0     += red_acc[tid]
        syn0     += syn_acc[tid]
    end

    SA_tot   ./= reps
    SAB_tot  ./= reps
    SA_pat   ./= reps
    SAB_pat  ./= reps
    Aeff_tot ./= reps
    Aeff_lcc ./= reps
    phi      ./= reps
    lcc      ./= reps
    ff       ./= reps

    red0 /= reps
    syn0 /= reps

    if nf > 0
        fits ./= nf
        base ./= nf
    else
        fits .= NaN
        base .= NaN
    end

    fitpack = (cT=fits[1], zT=fits[2], cP=fits[3], zP=fits[4],
               cEffTot=fits[5], zEffTot=fits[6], cEffLCC=fits[7], zEffLCC=fits[8],
               r2T=fits[9], r2P=fits[10], r2EffTot=fits[11], r2EffLCC=fits[12])
    basepack = (S0A_tot=base[1], S0AB_tot=base[2], S0A_pat=base[3], S0AB_pat=base[4])

    return (SA_tot=SA_tot, SAB_tot=SAB_tot, SA_pat=SA_pat, SAB_pat=SAB_pat,
            Aeff_tot=Aeff_tot, Aeff_lcc=Aeff_lcc,
            phi=phi, lcc=lcc, fragfail=ff,
            fits=fitpack, base=basepack,
            red0=red0, syn0=syn0)
end

# ============================================================
# Predictions for plotting
# ============================================================

function predictions(cfg::RunConfiguration, agg)
    T = length(cfg.fgrid)
    G = length(cfg.geoms)
    h = 1.0 .- cfg.fgrid

    Ssar_tot = isfinite(agg.fits.cT) && isfinite(agg.fits.zT) ? agg.fits.cT .* (h .^ agg.fits.zT) : fill(NaN, T)
    Ssar_pat = isfinite(agg.fits.cP) && isfinite(agg.fits.zP) ? agg.fits.cP .* (h .^ agg.fits.zP) : fill(NaN, T)

    Seff_tot = fill(NaN, G, T)
    Seff_pat = fill(NaN, G, T)

    if isfinite(agg.fits.cEffTot) && isfinite(agg.fits.zEffTot)
        @inbounds for gj in 1:G, t in 1:T
            Seff_tot[gj,t] = agg.fits.cEffTot * (max(1e-12, agg.Aeff_tot[gj,t]) ^ agg.fits.zEffTot)
        end
    end
    if isfinite(agg.fits.cEffLCC) && isfinite(agg.fits.zEffLCC)
        @inbounds for gj in 1:G, t in 1:T
            Seff_pat[gj,t] = agg.fits.cEffLCC * (max(1e-12, agg.Aeff_lcc[gj,t]) ^ agg.fits.zEffLCC)
        end
    end

    # clamp for nicer plots
    @inbounds for t in 1:T
        Ssar_tot[t] = clamp(Ssar_tot[t], 0.0, agg.base.S0A_tot)
        Ssar_pat[t] = clamp(Ssar_pat[t], 0.0, agg.base.S0A_pat)
    end
    @inbounds for gj in 1:G, t in 1:T
        Seff_tot[gj,t] = clamp(Seff_tot[gj,t], 0.0, agg.base.S0AB_tot)
        Seff_pat[gj,t] = clamp(Seff_pat[gj,t], 0.0, agg.base.S0AB_pat)
    end

    return (Ssar_tot=Ssar_tot, Ssar_pat=Ssar_pat, Seff_tot=Seff_tot, Seff_pat=Seff_pat)
end

# ============================================================
# FIGURE 1 (main)
# ============================================================

function fig_main(cfg::RunConfiguration, agg)
    pred = predictions(cfg, agg)
    G = length(cfg.geoms)

    fig = Figure(size=(1750, 980))

    for (j,g) in enumerate(cfg.geoms)
        ax = Axis(fig[1,j], xlabel="habitat loss f",
                  ylabel=(j==1 ? "richness (species)" : ""),
                  title="TOTAL-AREA   geometry=$(String(g))")
        lines!(ax, cfg.fgrid, vec(agg.SA_tot[j,:]);  label="A-only observed", linewidth=3)
        lines!(ax, cfg.fgrid, vec(agg.SAB_tot[j,:]); label="AB observed (biotic)", linewidth=3)
        lines!(ax, cfg.fgrid, pred.Ssar_tot;         label="SAR (area-only; A fit)", linestyle=:dash, linewidth=3)
        lines!(ax, cfg.fgrid, vec(pred.Seff_tot[j,:]); label="SAR_eff (supported area)", linestyle=:dot, linewidth=3)
        axislegend(ax; position=:lb)
    end

    for (j,g) in enumerate(cfg.geoms)
        ax = Axis(fig[2,j], xlabel="habitat loss f",
                  ylabel=(j==1 ? "richness (species)" : ""),
                  title="PATCH (r=0)   geometry=$(String(g))")
        lines!(ax, cfg.fgrid, vec(agg.SA_pat[j,:]);  label="A-only observed", linewidth=3)
        lines!(ax, cfg.fgrid, vec(agg.SAB_pat[j,:]); label="AB observed (biotic)", linewidth=3)
        lines!(ax, cfg.fgrid, pred.Ssar_pat;         label="SAR (area-only; A fit)", linestyle=:dash, linewidth=3)
        lines!(ax, cfg.fgrid, vec(pred.Seff_pat[j,:]); label="SAR_eff (supported LCC area)", linestyle=:dot, linewidth=3)
        axislegend(ax; position=:lb)
    end

    for (j,g) in enumerate(cfg.geoms)
        axL = Axis(fig[3,j], xlabel="habitat loss f",
                   ylabel=(j==1 ? "mean φ (consumers)" : ""),
                   title="Mechanism: support amount vs support connectivity")
        axR = Axis(fig[3,j], yaxisposition=:right, ylabel=(j==G ? "supported-LCC / frag-fail" : ""))
        linkxaxes!(axL, axR)
        hidespines!(axR); hidexdecorations!(axR)

        ylims!(axL, (0.0, 1.0))
        ylims!(axR, (0.0, 1.0))

        lines!(axL, cfg.fgrid, vec(agg.phi[j,:]); linewidth=3, label="mean φ")
        lines!(axR, cfg.fgrid, vec(agg.lcc[j,:]); linewidth=3, linestyle=:dash, label="supported-LCC fraction")
        lines!(axR, cfg.fgrid, vec(agg.fragfail[j,:]); linewidth=3, linestyle=:dot, label="frag-fail rate")
        axislegend(axL; position=:lc)
    end

    Label(fig[0,1:3],
          "Static synthetic test: SAR misses losses when trophically supported habitat fragments.  SAR_eff is a mechanistic bridge.",
          fontsize=16)

    display(fig)
end

# ============================================================
# FIGURE 2 (bridge scatter)
# ============================================================

function fig_bridge(cfg::RunConfiguration, agg)
    pred = predictions(cfg, agg)

    fig = Figure(size=(1750, 700))
    ax = Axis(fig[1,1:3],
        xlabel="supported-LCC fraction (consumer footprint)",
        ylabel="error in predicted losses (observed AB loss − predicted loss)",
        title="Mechanistic link: fragmentation of supported habitat organizes SAR prediction error"
    )

    for (j,g) in enumerate(cfg.geoms)
        x = vec(agg.lcc[j,:])

        # TOTAL losses
        LobsT = agg.base.S0AB_tot .- vec(agg.SAB_tot[j,:])
        LrawT = agg.base.S0A_tot  .- pred.Ssar_tot
        LeffT = agg.base.S0AB_tot .- vec(pred.Seff_tot[j,:])
        err_raw_T = LobsT .- LrawT
        err_eff_T = LobsT .- LeffT

        scatter!(ax, x, err_raw_T; markersize=10, label="$(String(g)) raw SAR (TOTAL)")
        scatter!(ax, x, err_eff_T; markersize=10, marker=:utriangle, label="$(String(g)) SAR_eff (TOTAL)")
        lines!(ax, x, err_raw_T; linewidth=2)
        lines!(ax, x, err_eff_T; linewidth=2, linestyle=:dash)

        # PATCH losses
        LobsP = agg.base.S0AB_pat .- vec(agg.SAB_pat[j,:])
        LrawP = agg.base.S0A_pat  .- pred.Ssar_pat
        LeffP = agg.base.S0AB_pat .- vec(pred.Seff_pat[j,:])
        err_raw_P = LobsP .- LrawP
        err_eff_P = LobsP .- LeffP

        scatter!(ax, x, err_raw_P; markersize=8, marker=:diamond, label="$(String(g)) raw SAR (PATCH)")
        scatter!(ax, x, err_eff_P; markersize=8, marker=:dtriangle, label="$(String(g)) SAR_eff (PATCH)")
        lines!(ax, x, err_raw_P; linewidth=2, linestyle=:dot)
        lines!(ax, x, err_eff_P; linewidth=2, linestyle=:dashdot)
    end

    axislegend(ax; position=:rt, nbanks=2)
    display(fig)
end

# ============================================================
# FIGURE 3 (movement sanity): r=0 vs r=2 patch connectivity
# ============================================================

function fig_patch_connectivity_sensitivity(cfg::RunConfiguration; seed::Int, reps::Int, rconn::Int, geoms::Vector{Symbol})
    println("\n=== Movement sanity check: patch connectivity radius r=$rconn (reps=$reps) ===")
    cfg2 = RunConfiguration(cfg.n1, cfg.n2, cfg.S, cfg.basal_frac, cfg.Lmax,
                            cfg.niche_sigma, cfg.niche_cut, cfg.match_sigma,
                            cfg.k_prey, cfg.select_sigma,
                            cfg.Emin_total, cfg.Epatch,
                            cfg.fgrid, geoms, cfg.sar_area_fracs, cfg.sar_samples)

    n1,n2 = cfg2.n1,cfg2.n2
    N = n1*n2
    up, down, left, right = build_neighbors4(n1,n2)
    row, col = make_rowcol(n1,n2)
    offs = radius_offsets(rconn)

    tbufs = [make_thread_buffers(cfg2) for _ in 1:nthreads()]
    G = length(cfg2.geoms)
    T = length(cfg2.fgrid)

    SA_pat_sum  = [zeros(Float64, G, T) for _ in 1:nthreads()]
    SAB_pat_sum = [zeros(Float64, G, T) for _ in 1:nthreads()]

    @threads for rep in 1:reps
        tid = threadid()
        rng = rep_rng(seed, rep)
        buf = tbufs[tid]

        env1, env2 = make_environment(rng; n1=n1, n2=n2, smooth_iters=30)
        A, mu1, mu2 = make_abiotic_maps(rng, env1, env2; S=cfg2.S, niche_sigma=cfg2.niche_sigma, niche_cut=cfg2.niche_cut)
        TL, basals, _ = make_trophic_levels(rng; S=cfg2.S, basal_frac=cfg2.basal_frac, Lmax=cfg2.Lmax)
        preylist = build_diets(rng; S=cfg2.S, TL=TL, basals=basals,
                               k_prey=cfg2.k_prey, match_sigma=cfg2.match_sigma, select_sigma=cfg2.select_sigma,
                               mu1=mu1, mu2=mu2)
        ords = [make_loss_order(rng, g, env1) for g in cfg2.geoms]

        keep0 = trues(N)

        for (gj,g) in enumerate(cfg2.geoms)
            ord = ords[gj]
            for (t,f) in enumerate(cfg2.fgrid)
                keep = (f == 0.0) ? keep0 : keep_from_order(ord, f, N)
                compute_P_and_counts!(buf.P, buf.Acount, buf.ABcount, A, keep, TL, preylist, buf.cons_present, buf.cons_idx, N)
                compute_maxpatches!(buf.Amax, buf.ABmax, A, buf.P, keep,
                                    buf.Acount, buf.ABcount, cfg2.Epatch, rconn,
                                    buf.present_flag, buf.present_idx, up, down, left, right,
                                    row, col, n1, n2, offs,
                                    buf.seen, buf.queue, buf.stamp, N)
                SA_pat_sum[tid][gj,t]  += count(>=(cfg2.Epatch), buf.Amax)
                SAB_pat_sum[tid][gj,t] += count(>=(cfg2.Epatch), buf.ABmax)
            end
        end
    end

    SA_pat  = zeros(Float64, G, T)
    SAB_pat = zeros(Float64, G, T)
    for tid in 1:nthreads()
        SA_pat  .+= SA_pat_sum[tid]
        SAB_pat .+= SAB_pat_sum[tid]
    end
    SA_pat  ./= reps
    SAB_pat ./= reps

    fig = Figure(size=(1400, 450))
    for (j,g) in enumerate(cfg2.geoms)
        ax = Axis(fig[1,j],
            xlabel="habitat loss f",
            ylabel=(j==1 ? "patch richness (r=$(rconn))" : ""),
            title="PATCH connectivity radius r=$(rconn) — geometry=$(String(g))"
        )
        lines!(ax, cfg2.fgrid, vec(SA_pat[j,:]);  label="A-only", linewidth=3)
        lines!(ax, cfg2.fgrid, vec(SAB_pat[j,:]); label="AB", linewidth=3)
        axislegend(ax; position=:lb)
    end
    display(fig)
end

# ============================================================
# LIGHT runner for sweeps (no per-species patch BFS):
#   returns total richness curves + footprint LCC fraction + SAR fit on A-only (total)
# ============================================================

function run_many_light(cfg::RunConfiguration; seed::Int, reps::Int)
    println("  LIGHT reps=$reps (threads=$(nthreads())) ...")
    n1,n2 = cfg.n1, cfg.n2
    up, down, left, right = build_neighbors4(n1,n2)
    row, col = make_rowcol(n1,n2)
    N = n1*n2

    G = length(cfg.geoms)
    T = length(cfg.fgrid)

    tbufs = [make_thread_buffers(cfg) for _ in 1:nthreads()]

    SA_tot_sum   = [zeros(Float64, G, T) for _ in 1:nthreads()]
    SAB_tot_sum  = [zeros(Float64, G, T) for _ in 1:nthreads()]
    lcc_sum      = [zeros(Float64, G, T) for _ in 1:nthreads()]

    fits_acc     = [zeros(Float64, 2) for _ in 1:nthreads()]  # cT,zT
    base_acc     = [zeros(Float64, 2) for _ in 1:nthreads()]  # S0A_tot, S0AB_tot
    red_acc      = zeros(Float64, nthreads())
    syn_acc      = zeros(Float64, nthreads())
    nfits        = zeros(Int, nthreads())

    @threads for rep in 1:reps
        tid = threadid()
        rng = rep_rng(seed, rep)
        buf = tbufs[tid]

        env1, env2 = make_environment(rng; n1=n1, n2=n2, smooth_iters=30)
        A, mu1, mu2 = make_abiotic_maps(rng, env1, env2; S=cfg.S, niche_sigma=cfg.niche_sigma, niche_cut=cfg.niche_cut)

        TL, basals, _ = make_trophic_levels(rng; S=cfg.S, basal_frac=cfg.basal_frac, Lmax=cfg.Lmax)
        preylist = build_diets(rng; S=cfg.S, TL=TL, basals=basals,
                               k_prey=cfg.k_prey, match_sigma=cfg.match_sigma, select_sigma=cfg.select_sigma,
                               mu1=mu1, mu2=mu2)

        ords = [make_loss_order(rng, g, env1) for g in cfg.geoms]

        keep0 = trues(N)
        compute_P_and_counts!(buf.P, buf.Acount, buf.ABcount, A, keep0, TL, preylist, buf.cons_present, buf.cons_idx, N)

        red0, syn0 = realised_proxies!(A, TL, preylist, buf.P, buf.Acount, buf.ABcount)
        red_acc[tid] += red0
        syn_acc[tid] += syn0

        S0A_tot  = count(>=(cfg.Emin_total), buf.Acount)
        S0AB_tot = count(>=(cfg.Emin_total), buf.ABcount)

        cT, zT, _ = sar_fit_total_Aonly(rng, A; n1=n1, n2=n2, Emin_total=cfg.Emin_total,
                                        area_fracs=cfg.sar_area_fracs, samples=cfg.sar_samples)

        if isfinite(cT) && isfinite(zT)
            fits_acc[tid][1] += cT
            fits_acc[tid][2] += zT
            base_acc[tid][1] += S0A_tot
            base_acc[tid][2] += S0AB_tot
            nfits[tid] += 1
        end

        for (gj,g) in enumerate(cfg.geoms)
            ord = ords[gj]
            for (t,f) in enumerate(cfg.fgrid)
                keep = (f == 0.0) ? keep0 : keep_from_order(ord, f, N)
                compute_P_and_counts!(buf.P, buf.Acount, buf.ABcount, A, keep, TL, preylist, buf.cons_present, buf.cons_idx, N)

                SA_tot_sum[tid][gj,t]  += count(>=(cfg.Emin_total), buf.Acount)
                SAB_tot_sum[tid][gj,t] += count(>=(cfg.Emin_total), buf.ABcount)

                footprint = length(buf.cons_idx)
                if footprint == 0
                    lcc_sum[tid][gj,t] += 0.0
                else
                    maxsz = lcc_max_size4(buf.cons_present, buf.cons_idx, up, down, left, right, buf.seen, buf.queue, buf.stamp)
                    lcc_sum[tid][gj,t] += maxsz / footprint
                end
            end
        end
    end

    SA_tot  = zeros(Float64, G, T)
    SAB_tot = zeros(Float64, G, T)
    lcc     = zeros(Float64, G, T)
    fits    = zeros(Float64, 2)
    base    = zeros(Float64, 2)
    nf      = sum(nfits)

    red0 = 0.0
    syn0 = 0.0

    for tid in 1:nthreads()
        SA_tot  .+= SA_tot_sum[tid]
        SAB_tot .+= SAB_tot_sum[tid]
        lcc     .+= lcc_sum[tid]
        fits    .+= fits_acc[tid]
        base    .+= base_acc[tid]
        red0    += red_acc[tid]
        syn0    += syn_acc[tid]
    end

    SA_tot  ./= reps
    SAB_tot ./= reps
    lcc     ./= reps
    red0    /= reps
    syn0    /= reps

    if nf > 0
        fits ./= nf
        base ./= nf
    else
        fits .= NaN
        base .= NaN
    end

    return (SA_tot=SA_tot, SAB_tot=SAB_tot, lcc=lcc,
            cT=fits[1], zT=fits[2], S0A_tot=base[1], S0AB_tot=base[2],
            red0=red0, syn0=syn0)
end

# ============================================================
# Sweep metrics + heatmaps
# ============================================================

@inline function trapz(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    @assert length(y) == n
    s = 0.0
    @inbounds for i in 1:n-1
        dx = float(x[i+1] - x[i])
        s += dx * (float(y[i]) + float(y[i+1])) / 2
    end
    return s
end

function bin2d_mean(x::AbstractArray{<:Real}, y::AbstractArray{<:Real}, v::AbstractArray{<:Real}; nbx::Int=NBX, nby::Int=NBY)
    xs = vec(Float64.(x)); ys = vec(Float64.(y)); vs = vec(Float64.(v))
    keep = findall(i -> isfinite(xs[i]) && isfinite(ys[i]) && isfinite(vs[i]), eachindex(vs))
    isempty(keep) && return (fill(NaN, nby, nbx), collect(1:nbx), collect(1:nby))

    xmin, xmax = minimum(xs[keep]), maximum(xs[keep])
    ymin, ymax = minimum(ys[keep]), maximum(ys[keep])
    if xmin == xmax; xmin -= 1e-6; xmax += 1e-6; end
    if ymin == ymax; ymin -= 1e-6; ymax += 1e-6; end

    xedges = range(xmin, xmax; length=nbx+1)
    yedges = range(ymin, ymax; length=nby+1)
    M = fill(NaN, nby, nbx)
    C = fill(0,   nby, nbx)

    for i in keep
        xi, yi, vi = xs[i], ys[i], vs[i]
        ix = clamp(searchsortedlast(xedges, xi), 1, nbx)
        iy = clamp(searchsortedlast(yedges, yi), 1, nby)
        if isnan(M[iy, ix])
            M[iy, ix] = vi
            C[iy, ix] = 1
        else
            M[iy, ix] += vi
            C[iy, ix] += 1
        end
    end

    for iy in 1:nby, ix in 1:nbx
        if C[iy, ix] > 0
            M[iy, ix] /= C[iy, ix]
        end
    end

    xcent = (xedges[1:end-1] .+ xedges[2:end]) ./ 2
    ycent = (yedges[1:end-1] .+ yedges[2:end]) ./ 2
    return M, xcent, ycent
end

function do_sweep_heatmaps()
    println("\n=== SWEEP (LIGHT): k_prey × select_sigma ===")
    cfg0 = RunConfiguration(
        N1_SWEEP, N2_SWEEP, S_SWEEP,
        BASAL_FRAC, LMAX, NICHE_SIGMA, NICHE_CUT,
        MATCH_SIGMA, KPREY_MAIN, SELECT_MAIN,
        EMIN_SWEEP, EPATCH_SWEEP,
        FGRID_SWEEP, GEOMS,
        SAR_AREAFRACS, max(12, SAR_SAMPLES ÷ 2)
    )

    nk = length(KPREY_LIST)
    ns = length(SELECT_LIST)
    G  = length(GEOMS)

    auc_frag = fill(NaN, G, nk, ns)
    peak_err = fill(NaN, G, nk, ns)
    wmiss    = fill(NaN, G, nk, ns)
    fcol     = fill(NaN, G, nk, ns)
    red0     = fill(NaN, nk, ns)
    syn0     = fill(NaN, nk, ns)

    settings = [(ik,is) for ik in 1:nk for is in 1:ns]

    @threads for sidx in 1:length(settings)
        ik, is = settings[sidx]

        cfg = RunConfiguration(
            cfg0.n1, cfg0.n2, cfg0.S,
            cfg0.basal_frac, cfg0.Lmax, cfg0.niche_sigma, cfg0.niche_cut,
            cfg0.match_sigma, KPREY_LIST[ik], SELECT_LIST[is],
            cfg0.Emin_total, cfg0.Epatch,
            cfg0.fgrid, cfg0.geoms,
            cfg0.sar_area_fracs, cfg0.sar_samples
        )

        out = run_many_light(cfg; seed=DEFAULT_SEED + 10_000*ik + 37*is, reps=REPS_SWEEP)

        # raw SAR prediction (TOTAL) from A-only fit at f=0 (sub-squares)
        h = 1.0 .- cfg.fgrid
        Ssar = out.cT .* (h .^ out.zT)
        Ssar = clamp.(Ssar, 0.0, out.S0A_tot)
        Lraw = out.S0A_tot .- Ssar

        red0[ik,is] = out.red0
        syn0[ik,is] = out.syn0

        for gj in 1:G
            frag = 1 .- vec(out.lcc[gj,:])
            auc_frag[gj,ik,is] = trapz(cfg.fgrid, frag) / (cfg.fgrid[end] - cfg.fgrid[1] + eps())

            Lobs = out.S0AB_tot .- vec(out.SAB_tot[gj,:])
            err  = Lobs .- Lraw
            peak_err[gj,ik,is] = maximum(err)
            wmiss[gj,ik,is]    = trapz(cfg.fgrid, err .* frag) / (cfg.fgrid[end] - cfg.fgrid[1] + eps())

            tcol = findfirst(t -> out.lcc[gj,t] < LCC_THRESH, 1:length(cfg.fgrid))
            fcol[gj,ik,is] = isnothing(tcol) ? cfg.fgrid[end] : cfg.fgrid[tcol]
        end
    end

    # --- Figure 4: PARAM axes ---
    fig = Figure(size=(1750, 980))
    Label(fig[0,1:6], "Sweep heatmaps (LIGHT): k_prey × select_sigma", fontsize=16)

    metric_blocks = [
        (:auc,  auc_frag, "AUC frag (mean over f of 1−supportedLCC)"),
        (:peak, peak_err, "peak raw SAR error (TOTAL losses)"),
        (:wm,   wmiss,    "weighted miss: AUC[err_raw(f) × frag(f)]"),
    ]

    for (r, (_, M, titlebase)) in enumerate(metric_blocks)
        for (gj,g) in enumerate(GEOMS)
            ax = Axis(fig[r,gj], title="$titlebase — geometry=$(String(g))",
                      xlabel="select_sigma", ylabel="k_prey")
            hm = heatmap!(ax, 1:ns, 1:nk, M[gj,:,:])
            ax.xticks = (1:ns, string.(SELECT_LIST))
            ax.yticks = (1:nk, string.(KPREY_LIST))
            Colorbar(fig[r,gj+3], hm)
        end
    end
    display(fig)

    # --- Figure 5: REALISED axes (binned) ---
    fig2 = Figure(size=(1750, 980))
    Label(fig2[0,1:6], "Sweep heatmaps (REALISED axes): redundancy × synchrony (binned)", fontsize=16)

    realised_metrics = [
        (auc_frag, "AUC frag"),
        (peak_err, "peak raw SAR error"),
        (wmiss,    "weighted miss"),
    ]

    for (r, (M, titlebase)) in enumerate(realised_metrics)
        for (gj,g) in enumerate(GEOMS)
            B, xcent, ycent = bin2d_mean(red0, syn0, M[gj,:,:]; nbx=NBX, nby=NBY)
            ax = Axis(fig2[r,gj], title="$titlebase — geometry=$(String(g))",
                      xlabel="realised redundancy (avg prey overlap / supported cells)",
                      ylabel="realised synchrony (mean φ at f=0)")
            hm = heatmap!(ax, xcent, ycent, B; interpolate=false)
            Colorbar(fig2[r,gj+3], hm)
        end
    end
    display(fig2)
end

# ============================================================
# MAIN
# ============================================================

function main()
    cfg = RunConfiguration(
        N1_MAIN, N2_MAIN, S_MAIN,
        BASAL_FRAC, LMAX,
        NICHE_SIGMA, NICHE_CUT,
        MATCH_SIGMA,
        KPREY_MAIN, SELECT_MAIN,
        EMIN_MAIN, EPATCH_MAIN,
        FGRID_MAIN,
        GEOMS,
        SAR_AREAFRACS,
        SAR_SAMPLES
    )

    println("==========================================================")
    println("BioLetters-oriented static synthetic pipeline (vNext FIXED)")
    println("grid=$(cfg.n1)x$(cfg.n2)  N=$(cfg.n1*cfg.n2)   S=$(cfg.S)")
    println("threads=$(nthreads())")
    println("==========================================================")

    agg = run_many(cfg; seed=DEFAULT_SEED, reps=REPS_MAIN)
    fig_main(cfg, agg)
    fig_bridge(cfg, agg)

    if DO_MOVEMENT_SANITY
        fig_patch_connectivity_sensitivity(cfg; seed=DEFAULT_SEED+999, reps=REPS_SENS, rconn=R_CONN_SENS, geoms=GEOMS_SENS)
    end

    if DO_SWEEP
        do_sweep_heatmaps()
    end

    return nothing
end

main()
