###############################################################
# postpass_auc_patchonly_36.jl
#
# Standalone sweep (same design as the 36-panel script) BUT:
#   - does NOT plot
#   - computes trapezoidal AUC of (A-only − AB) vs habitat loss f
#     under:
#       (i) Emin_patch criterion
#       (ii) viability criterion
#   - saves a tidy CSV:
#       ./Figures/patchonly_36/auc_summary.csv
#
# Run:
#   JULIA_NUM_THREADS=12 julia postpass_auc_patchonly_36.jl
###############################################################

using Random
using Statistics
using Printf
using Base.Threads

# =========================
# USER KNOBS (match your run)
# =========================
const DEFAULT_SEED = 1234

const N1      = 100
const N2      = 100
const S       = 200
const REPS    = 10

const FGRID   = collect(0.0:0.05:1.00)

const EMIN_PATCH = 80
const VFRAC      = 0.25

const BASAL_FRAC   = 0.25
const LMAX         = 5
const MATCH_SIGMA  = 1.0
const SELECT_SIGMA = 0.60

const KPREY_LIST = [1, 4, 9]
const CORR_LIST  = [0.0, 0.5, 1.0]
const GEOMS      = [:random, :cluster, :front]

const OUTDIR = joinpath(pwd(), "Figures", "patchonly_36")
const OUTCSV = joinpath(OUTDIR, "auc_summary.csv")

# =========================
# Utilities
# =========================
@inline idx(i::Int, j::Int, n2::Int) = (i - 1) * n2 + j

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

function ensure_dir(path::String)
    isdir(path) || mkpath(path)
    return path
end

# trapezoidal AUC over (x,y)
function auc_trapz(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n == length(y) || error("x/y length mismatch")
    a = 0.0
    @inbounds for i in 1:(n-1)
        dx = x[i+1] - x[i]
        a += 0.5 * (y[i] + y[i+1]) * dx
    end
    return a
end

# =========================
# Grid helpers + BFS LCC
# =========================
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

function lcc_max_size4(present_flag::BitVector, present_idx::Vector{Int},
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

# =========================
# Environment smoothing + z-score
# =========================
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

# =========================
# Niche scenarios (same 4)
# =========================
struct NicheScenario
    name::String
    sig_μlog::Float64
    sig_σlog::Float64
    sig_min::Float64
    sig_max::Float64
    occ_a::Float64
    occ_b::Float64
    occ_min::Float64
    occ_max::Float64
end

const NICHE_SCENARIOS = [
    NicheScenario("Very variable occupancy + very variable sigma",
                  log(0.85), 0.60, 0.20, 2.50,
                  0.7, 0.7, 0.05, 0.95),
    NicheScenario("Same occupancy range, less extreme sigma variability",
                  log(0.85), 0.35, 0.30, 2.00,
                  0.7, 0.7, 0.05, 0.95),
    NicheScenario("Skewed to narrower niches (lower mean occupancy), still variable",
                  log(0.60), 0.45, 0.15, 1.50,
                  1.2, 3.0, 0.02, 0.60),
    NicheScenario("Skewed to broader niches (higher mean occupancy), still variable",
                  log(1.10), 0.45, 0.25, 3.00,
                  3.0, 1.2, 0.40, 0.95)
]

function rand_beta(rng::AbstractRNG, a::Float64, b::Float64)
    x = randgamma(rng, a)
    y = randgamma(rng, b)
    return x / (x + y)
end

@inline clamp01(x::Float64) = x < 0 ? 0.0 : (x > 1 ? 1.0 : x)

function make_abiotic_maps_quantile!(A::BitMatrix,
                                    mu1::Vector{Float64}, mu2::Vector{Float64},
                                    sig::Vector{Float64}, occ::Vector{Float64},
                                    rng::AbstractRNG, env1::Matrix{Float64}, env2::Matrix{Float64},
                                    scen::NicheScenario)

    n1, n2 = size(env1)
    N = n1*n2
    S = size(A, 1)

    mu1 .= randn(rng, S)
    mu2 .= randn(rng, S)

    @inbounds for i in 1:S
        si = exp(scen.sig_μlog + scen.sig_σlog * randn(rng))
        sig[i] = clamp(si, scen.sig_min, scen.sig_max)

        oi = rand_beta(rng, scen.occ_a, scen.occ_b)
        occ[i] = clamp(oi, scen.occ_min, scen.occ_max)
    end

    suits  = Vector{Float64}(undef, N)
    sorted = Vector{Float64}(undef, N)

    @inbounds for i in 1:S
        si = sig[i]
        denom = 2.0 * si * si

        k = 0
        for r in 1:n1, c in 1:n2
            k += 1
            d1 = env1[r,c] - mu1[i]
            d2 = env2[r,c] - mu2[i]
            suits[k] = exp(-(d1*d1 + d2*d2) / denom)
        end

        copyto!(sorted, suits)
        sort!(sorted)

        q = clamp01(1.0 - occ[i])
        j = clamp(Int(floor(q*(N-1))) + 1, 1, N)
        thr = sorted[j]

        for k in 1:N
            A[i,k] = (suits[k] >= thr)
        end
    end
    return nothing
end

# =========================
# Trophic + diets w/ correlation knob
# =========================
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
                     mu1::Vector{Float64}, mu2::Vector{Float64},
                     corr_strength::Float64)

    preylist = [Int[] for _ in 1:S]
    α = clamp(corr_strength, 0.0, 1.0)

    @inbounds for i in 1:S
        TL[i] == 1 && continue
        cands = [j for j in 1:S if TL[j] < TL[i]]
        isempty(cands) && (cands = copy(basals))

        wmatch = Vector{Float64}(undef, length(cands))
        for (t, j) in enumerate(cands)
            d1 = mu1[i]-mu1[j]; d2 = mu2[i]-mu2[j]
            w0 = exp(-(d1*d1 + d2*d2) / (2.0*match_sigma^2))
            wmatch[t] = w0^α  # α=0 => 1 (no matching), α=1 => full matching
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
            w0 = exp(-(d1i*d1i + d2i*d2i) / (2.0*match_sigma^2))
            wfeas = w0^α

            w[t] = wguild * wfeas
        end

        others = sample_weighted_no_replace(rng, cands2, w, k_prey-1)
        preylist[i] = vcat([g], others)
    end

    return preylist
end

# =========================
# Habitat-loss orders
# =========================
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

# =========================
# AB support + patches
# =========================
mutable struct Buffers
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

function make_buffers(S::Int, N::Int)
    return Buffers(
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

function compute_P_and_counts!(buf::Buffers,
                              A::BitMatrix, keep::BitVector,
                              TL::Vector{Int}, preylist::Vector{Vector{Int}},
                              N::Int)

    S = size(A, 1)
    fill!(buf.Acount, 0)
    fill!(buf.ABcount, 0)
    buf.P .= false

    @inbounds for i in 1:S
        c = 0
        for k in 1:N
            c += (keep[k] & A[i,k]) ? 1 : 0
        end
        buf.Acount[i] = c
    end

    order = sortperm(TL)
    @inbounds for i in order
        if TL[i] == 1
            c = 0
            for k in 1:N
                v = keep[k] & A[i,k]
                buf.P[i,k] = v
                c += v ? 1 : 0
            end
            buf.ABcount[i] = c
        else
            prey = preylist[i]
            c = 0
            for k in 1:N
                sup = false
                for pj in prey
                    sup |= buf.P[pj,k]
                    sup && break
                end
                v = keep[k] & A[i,k] & sup
                buf.P[i,k] = v
                c += v ? 1 : 0
            end
            buf.ABcount[i] = c
        end
    end
    return nothing
end

function compute_maxpatches!(buf::Buffers,
                            A::BitMatrix, keep::BitVector,
                            up::Vector{Int}, down::Vector{Int}, left::Vector{Int}, right::Vector{Int},
                            N::Int)

    S = size(A, 1)

    @inbounds for i in 1:S
        if buf.Acount[i] == 0
            buf.Amax[i] = 0
        else
            fill_present_from_A!(buf.present_flag, buf.present_idx, A, i, keep, N)
            buf.Amax[i] = lcc_max_size4(buf.present_flag, buf.present_idx, up, down, left, right,
                                        buf.seen, buf.queue, buf.stamp)
            clear_present!(buf.present_flag, buf.present_idx)
        end

        if buf.ABcount[i] == 0
            buf.ABmax[i] = 0
        else
            fill_present_from_P!(buf.present_flag, buf.present_idx, buf.P, i, N)
            buf.ABmax[i] = lcc_max_size4(buf.present_flag, buf.present_idx, up, down, left, right,
                                         buf.seen, buf.queue, buf.stamp)
            clear_present!(buf.present_flag, buf.present_idx)
        end
    end
    return nothing
end

# =========================
# One rep + combo aggregation
# =========================
function run_one_rep!(rng::AbstractRNG,
                      scen::NicheScenario, k_prey::Int, corr_strength::Float64,
                      up::Vector{Int}, down::Vector{Int}, left::Vector{Int}, right::Vector{Int},
                      buf::Buffers)

    n1, n2 = N1, N2
    N = n1*n2
    T = length(FGRID)
    G = length(GEOMS)

    env1, env2 = make_environment(rng; n1=n1, n2=n2, smooth_iters=30)

    A = BitMatrix(undef, S, N)
    mu1 = Vector{Float64}(undef, S)
    mu2 = Vector{Float64}(undef, S)
    sig = Vector{Float64}(undef, S)
    occ = Vector{Float64}(undef, S)
    make_abiotic_maps_quantile!(A, mu1, mu2, sig, occ, rng, env1, env2, scen)

    TL, basals, _ = make_trophic_levels(rng; S=S, basal_frac=BASAL_FRAC, Lmax=LMAX)
    preylist = build_diets(rng; S=S, TL=TL, basals=basals,
                           k_prey=k_prey, match_sigma=MATCH_SIGMA, select_sigma=SELECT_SIGMA,
                           mu1=mu1, mu2=mu2,
                           corr_strength=corr_strength)

    ords = [make_loss_order(rng, g, env1) for g in GEOMS]

    keep0 = trues(N)
    compute_P_and_counts!(buf, A, keep0, TL, preylist, N)
    compute_maxpatches!(buf, A, keep0, up, down, left, right, N)

    baseA  = copy(buf.Amax)
    baseAB = copy(buf.ABmax)

    thrA  = Vector{Int}(undef, S)
    thrAB = Vector{Int}(undef, S)
    @inbounds for i in 1:S
        thrA[i]  = max(EMIN_PATCH, round(Int, VFRAC * baseA[i]))
        thrAB[i] = max(EMIN_PATCH, round(Int, VFRAC * baseAB[i]))
    end

    A_Emin   = zeros(Float64, G, T)
    AB_Emin  = zeros(Float64, G, T)
    A_Viab   = zeros(Float64, G, T)
    AB_Viab  = zeros(Float64, G, T)

    for (gj, g) in enumerate(GEOMS)
        ord = ords[gj]
        for (t, f) in enumerate(FGRID)
            keep = (f == 0.0) ? keep0 : keep_from_order(ord, f, N)

            compute_P_and_counts!(buf, A, keep, TL, preylist, N)
            compute_maxpatches!(buf, A, keep, up, down, left, right, N)

            A_Emin[gj,t]  = count(>=(EMIN_PATCH), buf.Amax)
            AB_Emin[gj,t] = count(>=(EMIN_PATCH), buf.ABmax)

            vA = 0; vB = 0
            @inbounds for i in 1:S
                vA += (buf.Amax[i]  >= thrA[i])
                vB += (buf.ABmax[i] >= thrAB[i])
            end
            A_Viab[gj,t]  = vA
            AB_Viab[gj,t] = vB
        end
    end

    return (A_Emin=A_Emin, AB_Emin=AB_Emin, A_Viab=A_Viab, AB_Viab=AB_Viab)
end

function run_combo(scen::NicheScenario, k_prey::Int, corr_strength::Float64; seed::Int)
    N = N1*N2
    up, down, left, right = build_neighbors4(N1, N2)

    tbufs = [make_buffers(S, N) for _ in 1:nthreads()]

    G = length(GEOMS)
    T = length(FGRID)

    A_Emin_sum  = [zeros(Float64, G, T) for _ in 1:nthreads()]
    AB_Emin_sum = [zeros(Float64, G, T) for _ in 1:nthreads()]
    A_Viab_sum  = [zeros(Float64, G, T) for _ in 1:nthreads()]
    AB_Viab_sum = [zeros(Float64, G, T) for _ in 1:nthreads()]

    @threads for rep in 1:REPS
        tid = threadid()
        rng = rep_rng(seed, rep)
        out = run_one_rep!(rng, scen, k_prey, corr_strength, up, down, left, right, tbufs[tid])
        A_Emin_sum[tid]  .+= out.A_Emin
        AB_Emin_sum[tid] .+= out.AB_Emin
        A_Viab_sum[tid]  .+= out.A_Viab
        AB_Viab_sum[tid] .+= out.AB_Viab
    end

    A_Emin  = zeros(Float64, G, T)
    AB_Emin = zeros(Float64, G, T)
    A_Viab  = zeros(Float64, G, T)
    AB_Viab = zeros(Float64, G, T)

    for tid in 1:nthreads()
        A_Emin  .+= A_Emin_sum[tid]
        AB_Emin .+= AB_Emin_sum[tid]
        A_Viab  .+= A_Viab_sum[tid]
        AB_Viab .+= AB_Viab_sum[tid]
    end

    A_Emin  ./= REPS
    AB_Emin ./= REPS
    A_Viab  ./= REPS
    AB_Viab ./= REPS

    return (A_Emin=A_Emin, AB_Emin=AB_Emin, A_Viab=A_Viab, AB_Viab=AB_Viab)
end

# =========================
# MAIN: compute AUCs + write CSV
# =========================
function main()
    ensure_dir(OUTDIR)

    open(OUTCSV, "w") do io
        println(io, "scenario_id,scenario_name,k_prey,corr,geometry,auc_emin,auc_viability,aucmean_emin_over_geoms,aucmean_viab_over_geoms")

        combo_i = 0
        for (si, scen) in enumerate(NICHE_SCENARIOS)
            for k_prey in KPREY_LIST
                for α in CORR_LIST
                    combo_i += 1
                    seed = DEFAULT_SEED + 100_000*si + 1_000*k_prey + round(Int, 100*α)

                    println(@sprintf("[%2d/36] AUCs: scenario=%d | k_prey=%d | corr=%.2f", combo_i, si, k_prey, α))
                    out = run_combo(scen, k_prey, α; seed=seed)

                    # AUC per geometry
                    aucE = zeros(Float64, length(GEOMS))
                    aucV = zeros(Float64, length(GEOMS))
                    for (gj, g) in enumerate(GEOMS)
                        yE = vec(out.A_Emin[gj,:] .- out.AB_Emin[gj,:])
                        yV = vec(out.A_Viab[gj,:] .- out.AB_Viab[gj,:])
                        aucE[gj] = auc_trapz(FGRID, yE)
                        aucV[gj] = auc_trapz(FGRID, yV)
                    end

                    meanE = mean(aucE)
                    meanV = mean(aucV)

                    for (gj, g) in enumerate(GEOMS)
                        # quote scenario_name safely (no commas assumption is ok; if you add commas, switch to TSV)
                        sname = replace(scen.name, "\"" => "'")
                        println(io, @sprintf("%d,\"%s\",%d,%.2f,%s,%.6f,%.6f,%.6f,%.6f",
                                            si, sname, k_prey, α, String(g),
                                            aucE[gj], aucV[gj], meanE, meanV))
                    end
                end
            end
        end
    end

    println("\nSaved AUC summary to:\n  $OUTCSV")
    return nothing
end

main()

