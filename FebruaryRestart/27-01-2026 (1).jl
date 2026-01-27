###############################################################
# suitHab_PATCH_robust_niches_viability.jl
#
# Goal:
#   Explore ONLY the PATCH richness curves (r=0 connectivity),
#   keeping your Emin thresholds but adding an explicit viability
#   definition (relative-to-baseline).
#
# Key robustness changes:
#   1) Variable niche breadths across species (sigma_i lognormal)
#   2) Occupancy-calibrated niches: each species is assigned a
#      target occupancy p_i (very variable across species), and
#      we set its niche cutoff so that fraction of suitable cells
#      at f=0 matches p_i (approximately exact).
#      => "niche construction" becomes robust: sigma/cut knobs
#         do NOT arbitrarily change baseline occupancy.
#
# What it plots (ONLY PATCH):
#   For each niche parameterization:
#     3 columns: geometry = random / cluster / front
#     2 rows:
#       Row 1: Patch richness using Emin_patch (absolute cells)
#       Row 2: Patch richness using Viability (fraction of baseline patch)
#   In each panel: A-only and AB curves.
#
# Run:
#   JULIA_NUM_THREADS=12 julia suitHab_PATCH_robust_niches_viability.jl
###############################################################

using Random
using Statistics
using CairoMakie
using Base.Threads
using Printf

# ============================================================
# USER KNOBS
# ============================================================

const DEFAULT_SEED = 1234

# Grid / community size
const N1 = 120
const N2 = 120
const S  = 250

# Replicates (increase for smoother curves)
const REPS = 12

# Habitat-loss grid (include 1.0 if you want the end behavior)
const FGRID = collect(0.0:0.05:1.0)

# Patch threshold (your original meaning)
const EMIN_PATCH = 80

# Viability: require a fraction of baseline patch size
# (still keeps EMIN_PATCH as a floor so very tiny baselines don't dominate)
const VFRAC_PATCH = 0.25

# Trophic construction
const BASAL_FRAC  = 0.25
const LMAX        = 5
const MATCH_SIGMA = 1.00
const KPREY       = 3
const SELECT_SIGMA = 0.60

# Habitat loss geometries
const GEOMS = [:random, :cluster, :front]

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
    s = splitmix64(
        UInt64(master_seed) ⊻
        splitmix64(reinterpret(UInt64, Int64(rep)))
    )
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
# Robust niche construction:
#   - sigma_i varies strongly across species (lognormal)
#   - each species gets a target occupancy p_i (very variable)
#   - cutoff is set per species to hit that occupancy at f=0
# ============================================================

# simple Beta sampler via Gamma (no external packages)
@inline function rand_beta(rng::AbstractRNG, a::Float64, b::Float64)
    x = rand(rng, Gamma(a, 1.0))
    y = rand(rng, Gamma(b, 1.0))
    return x / (x + y)
end

# fallback Gamma if not imported: implement using Random's randexp? -> better: use Distributions
# But we promised no extra deps. So implement Marsaglia-Tsang for Gamma(k,1).
# Works for k>0.
function rand_gamma_mt(rng::AbstractRNG, k::Float64)
    if k < 1.0
        # Johnk's generator via boosting
        u = rand(rng)
        return rand_gamma_mt(rng, k + 1.0) * u^(1.0/k)
    end
    d = k - 1.0/3.0
    c = 1.0 / sqrt(9.0d)
    while true
        x = randn(rng)
        v = (1.0 + c*x)
        v <= 0 && continue
        v = v^3
        u = rand(rng)
        if u < 1.0 - 0.331*(x^4)
            return d*v
        end
        if log(u) < 0.5*x^2 + d*(1.0 - v + log(v))
            return d*v
        end
    end
end

@inline function rand_beta_mt(rng::AbstractRNG, a::Float64, b::Float64)
    x = rand_gamma_mt(rng, a)
    y = rand_gamma_mt(rng, b)
    return x / (x + y)
end

struct NicheParam
    name::String
    sigma_base::Float64      # typical sigma in env space
    sigma_logsd::Float64     # variability of sigma across species (lognormal sd)
    occ_a::Float64           # Beta(a,b) for target occupancy p_i
    occ_b::Float64
    pmin::Float64
    pmax::Float64
end

function make_abiotic_maps_robust(rng::AbstractRNG, env1::Matrix{Float64}, env2::Matrix{Float64};
                                 S::Int, par::NicheParam)

    n1, n2 = size(env1)
    N = n1*n2

    # species niche centers in env space
    mu1 = randn(rng, S)
    mu2 = randn(rng, S)

    # strong interspecific variability in breadth
    # sigma_i ~ LogNormal(log(sigma_base), sigma_logsd)
    sigma_i = exp.(log(par.sigma_base) .+ par.sigma_logsd .* randn(rng, S))

    # target occupancy fractions (variable across species)
    p_i = Vector{Float64}(undef, S)
    @inbounds for i in 1:S
        p = rand_beta_mt(rng, par.occ_a, par.occ_b)
        p_i[i] = clamp(par.pmin + (par.pmax - par.pmin)*p, 1e-4, 0.9999)
    end

    A = BitMatrix(undef, S, N)

    # For each species:
    # 1) compute suitability over all cells
    # 2) choose cutoff as quantile so that fraction >= cutoff equals p_i
    suit = Vector{Float64}(undef, N)

    @inbounds for i in 1:S
        si = sigma_i[i]
        inv2s2 = 1.0 / (2.0*si*si)

        for r in 1:n1, c in 1:n2
            d1 = env1[r,c] - mu1[i]
            d2 = env2[r,c] - mu2[i]
            suit[idx(r,c,n2)] = exp(-(d1*d1 + d2*d2) * inv2s2)
        end

        # Find cutoff to match occupancy fraction p_i:
        # want top p_i fraction of cells to be TRUE
        kkeep = clamp(round(Int, p_i[i]*N), 1, N)
        # partial selection: get kth largest by sorting indices (N is 14400 -> OK)
        ord = sortperm(suit; rev=true)
        cutoff = suit[ord[kkeep]]

        for k in 1:N
            A[i,k] = (suit[k] >= cutoff)
        end
    end

    return A, mu1, mu2, sigma_i, p_i
end

# ============================================================
# Trophic levels + diets
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

        # (1) pick guild-centre prey based on consumer-prey match
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
# AB support map P[S,N] + counts
# ============================================================

function compute_P_and_counts!(P::BitMatrix,
                              Acount::Vector{Int}, ABcount::Vector{Int},
                              A::BitMatrix, keep::BitVector,
                              TL::Vector{Int}, preylist::Vector{Vector{Int}},
                              N::Int)

    S = size(A,1)
    fill!(Acount, 0)
    fill!(ABcount, 0)
    P .= false

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
                c += v ? 1 : 0
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
                            Epatch::Int,
                            present_flag::BitVector, present_idx::Vector{Int},
                            up::Vector{Int}, down::Vector{Int}, left::Vector{Int}, right::Vector{Int},
                            seen::Vector{Int}, queue::Vector{Int}, stamp::Base.RefValue{Int},
                            N::Int)

    S = length(Acount)
    @inbounds for i in 1:S
        # A-only
        if Acount[i] < Epatch
            Amax[i] = Acount[i]
        else
            fill_present_from_A!(present_flag, present_idx, A, i, keep, N)
            Amax[i] = lcc_max_size4(present_flag, present_idx, up, down, left, right, seen, queue, stamp)
            clear_present!(present_flag, present_idx)
        end

        # AB
        if ABcount[i] < Epatch
            ABmax[i] = ABcount[i]
        else
            fill_present_from_P!(present_flag, present_idx, P, i, N)
            ABmax[i] = lcc_max_size4(present_flag, present_idx, up, down, left, right, seen, queue, stamp)
            clear_present!(present_flag, present_idx)
        end
    end
    return nothing
end

# ============================================================
# Thread buffers
# ============================================================

mutable struct ThreadBuffers1
    Acount::Vector{Int}
    ABcount::Vector{Int}
    Amax::Vector{Int}
    ABmax::Vector{Int}
    P::BitMatrix
    present_flag::BitVector
    present_idx::Vector{Int}
    seen::Vector{Int}
    queue::Vector{Int}
    stamp::Base.RefValue{Int}
end

function make_thread_buffers(S::Int, N::Int)
    return ThreadBuffers1(
        zeros(Int, S),
        zeros(Int, S),
        zeros(Int, S),
        zeros(Int, S),
        falses(S, N),
        falses(N),
        Int[],
        zeros(Int, N),
        Vector{Int}(undef, N),
        Ref(0)
    )
end

# ============================================================
# One replicate: PATCH curves (Epatch + viability)
# ============================================================

function run_one_rep_patch!(rng::AbstractRNG, par::NicheParam,
                            up, down, left, right,
                            buf::ThreadBuffers1)

    n1, n2 = N1, N2
    N = n1*n2
    T = length(FGRID)
    G = length(GEOMS)

    env1, env2 = make_environment(rng; n1=n1, n2=n2, smooth_iters=30)
    A, mu1, mu2, sigma_i, p_i = make_abiotic_maps_robust(rng, env1, env2; S=S, par=par)

    TL, basals, _ = make_trophic_levels(rng; S=S, basal_frac=BASAL_FRAC, Lmax=LMAX)
    preylist = build_diets(rng; S=S, TL=TL, basals=basals,
                           k_prey=KPREY, match_sigma=MATCH_SIGMA, select_sigma=SELECT_SIGMA,
                           mu1=mu1, mu2=mu2)

    ords = [make_loss_order(rng, g, env1) for g in GEOMS]
    keep0 = trues(N)

    # baseline f=0 to set viability thresholds per species
    compute_P_and_counts!(buf.P, buf.Acount, buf.ABcount, A, keep0, TL, preylist, N)
    compute_maxpatches!(buf.Amax, buf.ABmax, A, buf.P, keep0, buf.Acount, buf.ABcount,
                        EMIN_PATCH, buf.present_flag, buf.present_idx,
                        up, down, left, right,
                        buf.seen, buf.queue, buf.stamp, N)

    Amax0  = copy(buf.Amax)
    ABmax0 = copy(buf.ABmax)

    # Viability threshold per species (keep Emin as floor)
    vthA  = Vector{Int}(undef, S)
    vthAB = Vector{Int}(undef, S)
    @inbounds for i in 1:S
        vthA[i]  = max(EMIN_PATCH, ceil(Int, VFRAC_PATCH * max(1, Amax0[i])))
        vthAB[i] = max(EMIN_PATCH, ceil(Int, VFRAC_PATCH * max(1, ABmax0[i])))
    end

    # outputs: richness for A-only vs AB under (Epatch) and (viability)
    SA_E   = zeros(Float64, G, T)
    SAB_E  = zeros(Float64, G, T)
    SA_V   = zeros(Float64, G, T)
    SAB_V  = zeros(Float64, G, T)

    for (gj,g) in enumerate(GEOMS)
        ord = ords[gj]
        for (t,f) in enumerate(FGRID)
            keep = (f == 0.0) ? keep0 : keep_from_order(ord, f, N)

            compute_P_and_counts!(buf.P, buf.Acount, buf.ABcount, A, keep, TL, preylist, N)
            compute_maxpatches!(buf.Amax, buf.ABmax, A, buf.P, keep, buf.Acount, buf.ABcount,
                                EMIN_PATCH, buf.present_flag, buf.present_idx,
                                up, down, left, right,
                                buf.seen, buf.queue, buf.stamp, N)

            # Emin-only
            SA_E[gj,t]  = count(>=(EMIN_PATCH), buf.Amax)
            SAB_E[gj,t] = count(>=(EMIN_PATCH), buf.ABmax)

            # Viability (relative-to-baseline thresholds)
            cA = 0; cAB = 0
            @inbounds for i in 1:S
                cA  += (buf.Amax[i]  >= vthA[i])
                cAB += (buf.ABmax[i] >= vthAB[i])
            end
            SA_V[gj,t]  = cA
            SAB_V[gj,t] = cAB
        end
    end

    # quick diagnostics to confirm niche variability (returned, optional)
    return (SA_E=SA_E, SAB_E=SAB_E, SA_V=SA_V, SAB_V=SAB_V,
            occ_mean=mean(p_i), occ_sd=std(p_i), sigma_mean=mean(sigma_i), sigma_sd=std(sigma_i))
end

# ============================================================
# Aggregate over reps (threads)
# ============================================================

function run_many_patch(par::NicheParam; seed::Int=DEFAULT_SEED, reps::Int=REPS)
    println("\nRunning PATCH-only (reps=$reps, threads=$(nthreads())) — niche: $(par.name)")
    n1,n2 = N1,N2
    N = n1*n2
    up, down, left, right = build_neighbors4(n1,n2)

    G = length(GEOMS)
    T = length(FGRID)

    tbufs = [make_thread_buffers(S, N) for _ in 1:nthreads()]

    SA_E_sum  = [zeros(Float64, G, T) for _ in 1:nthreads()]
    SAB_E_sum = [zeros(Float64, G, T) for _ in 1:nthreads()]
    SA_V_sum  = [zeros(Float64, G, T) for _ in 1:nthreads()]
    SAB_V_sum = [zeros(Float64, G, T) for _ in 1:nthreads()]

    occ_m = zeros(Float64, nthreads())
    occ_s = zeros(Float64, nthreads())
    sig_m = zeros(Float64, nthreads())
    sig_s = zeros(Float64, nthreads())

    @threads for rep in 1:reps
        tid = threadid()
        h = Int(mod(hash(par.name), UInt64(typemax(Int))))
        rng = rep_rng(seed, rep + 10_000 * h)

        out = run_one_rep_patch!(rng, par, up, down, left, right, tbufs[tid])

        SA_E_sum[tid]  .+= out.SA_E
        SAB_E_sum[tid] .+= out.SAB_E
        SA_V_sum[tid]  .+= out.SA_V
        SAB_V_sum[tid] .+= out.SAB_V

        occ_m[tid] += out.occ_mean
        occ_s[tid] += out.occ_sd
        sig_m[tid] += out.sigma_mean
        sig_s[tid] += out.sigma_sd
    end

    SA_E  = zeros(Float64, G, T)
    SAB_E = zeros(Float64, G, T)
    SA_V  = zeros(Float64, G, T)
    SAB_V = zeros(Float64, G, T)

    for tid in 1:nthreads()
        SA_E  .+= SA_E_sum[tid]
        SAB_E .+= SAB_E_sum[tid]
        SA_V  .+= SA_V_sum[tid]
        SAB_V .+= SAB_V_sum[tid]
    end

    SA_E  ./= reps
    SAB_E ./= reps
    SA_V  ./= reps
    SAB_V ./= reps

    occ_mean = sum(occ_m) / reps
    occ_sd   = sum(occ_s) / reps
    sigma_mean = sum(sig_m) / reps
    sigma_sd   = sum(sig_s) / reps

    println(@sprintf("  niche diagnostics (avg over reps):  occ mean=%.3f  occ sd=%.3f   sigma mean=%.3f  sigma sd=%.3f",
                     occ_mean, occ_sd, sigma_mean, sigma_sd))

    return (SA_E=SA_E, SAB_E=SAB_E, SA_V=SA_V, SAB_V=SAB_V,
            occ_mean=occ_mean, occ_sd=occ_sd, sigma_mean=sigma_mean, sigma_sd=sigma_sd)
end

# ============================================================
# Plotting: only PATCH figure, 2 rows × 3 columns
# ============================================================
function fig_patch_only(par::NicheParam, agg)
    fig = Figure(size=(1750, 820))
    Label(fig[0, 1:3],
          "PATCH-only — Emin_patch=$(EMIN_PATCH) and Viability (≥ $(VFRAC_PATCH)×baseline patch, with Emin floor).  Niche: $(par.name)",
          fontsize=16)

    # Row 1: Emin_patch
    for (j,g) in enumerate(GEOMS)
        ax = Axis(fig[1, j],
            xlabel="habitat loss f",
            ylabel=(j==1 ? "patch richness" : ""),
            title="PATCH (Emin) — geometry=$(String(g))"
        )
        lines!(ax, FGRID, vec(agg.SA_E[j,:]);  linewidth=3, label="A-only")
        lines!(ax, FGRID, vec(agg.SAB_E[j,:]); linewidth=3, label="AB")
        axislegend(ax; position=:lb)
    end

    # Row 2: viability
    for (j,g) in enumerate(GEOMS)
        ax = Axis(fig[2, j],
            xlabel="habitat loss f",
            ylabel=(j==1 ? "patch richness" : ""),
            title="PATCH (viability) — geometry=$(String(g))"
        )
        lines!(ax, FGRID, vec(agg.SA_V[j,:]);  linewidth=3, label="A-only")
        lines!(ax, FGRID, vec(agg.SAB_V[j,:]); linewidth=3, label="AB")
        axislegend(ax; position=:lb)
    end

    display(fig)
    return fig
end

# ============================================================
# MAIN: try a few niche parameterizations to see shape robustness
# ============================================================

function main()
    println("==========================================================")
    println("PATCH robustness explorer (niches + viability)")
    println("grid=$(N1)x$(N2)  N=$(N1*N2)   S=$(S)   reps=$(REPS)")
    println("EMIN_PATCH=$(EMIN_PATCH)   VFRAC_PATCH=$(VFRAC_PATCH)")
    println("==========================================================")

    # A small set of niche parameterizations:
    # - All have *highly variable* occupancy via Beta,
    # - and also *highly variable* sigma via lognormal,
    # but differ in how extreme variability is.
    pars = NicheParam[
        NicheParam("Very variable occupancy + very variable sigma",
                   0.85, 0.80,   0.45, 0.45,   0.02, 0.85),
        NicheParam("Same occupancy range, slightly less extreme sigma variability",
                   0.85, 0.45,   0.45, 0.45,   0.02, 0.85),
        NicheParam("Skewed to narrower niches (lower mean occupancy), still variable",
                   0.85, 0.80,   0.35, 1.20,   0.01, 0.60),
    ]

    for (k,par) in enumerate(pars)
        agg = run_many_patch(par; seed=DEFAULT_SEED + 1000*k, reps=REPS)
        fig_patch_only(par, agg)
    end

    println("\nDone. Increase REPS for smoother averages.")
    return nothing
end

main()
