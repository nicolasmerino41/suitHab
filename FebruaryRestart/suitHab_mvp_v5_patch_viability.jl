# suitHab_MVP_static_v5_patch_viability.jl
#
# vNext = Static MVP + patch-viability criterion (in addition to total-area criterion)
#
# What it does (in one run):
#   - Build synthetic env + abiotic maps A_i(x)
#   - Build synthetic metaweb (TL + diets) with fixed local prey rule: ">= 1 prey present in cell"
#   - Apply habitat loss geometries (random / cluster / front) along a nested removal order
#   - Compute BOTH extinction definitions, side-by-side (no switches):
#       (1) Total-area viability:  area >= Emin
#       (2) Patch viability:      largest connected patch >= Epatch
#     for BOTH A-only and AB (biotic-supported) outcomes.
#   - Fit two SAR baselines from pre-loss A-only:
#       SAR_total  fitted with total-area presence rule
#       SAR_patch  fitted with patch presence rule
#     and compare SAR predictions against AB "reality" under each criterion.
#   - Mechanistic panels:
#       row3: mean φ (consumers), LCC fraction (consumer-supported overlap), and fragmentation-failure rate
#       row4: SAR error vs LCC, for BOTH criteria and all geometries
#

using Random
using Statistics
using Printf
using CairoMakie

# ----------------------------
# Utilities
# ----------------------------
@inline idx(i::Int, j::Int, n2::Int) = (i - 1) * n2 + j

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

function weighted_pick_index(rng::AbstractRNG, w::Vector{Float64})
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

function auc_trapz(x::Vector{Float64}, y::Vector{Float64})
    s = 0.0
    @inbounds for i in 1:length(x)-1
        dx = x[i+1] - x[i]
        s += dx * (y[i] + y[i+1]) / 2.0
    end
    return s
end

# ----------------------------
# Neighbors for patch connectivity (4-neighborhood on the full grid)
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

# Compute largest component size given:
#   present_flag[k] = true if cell k is present
#   present_idx = list of indices where present_flag is true
# Uses stamp-based visited marking to avoid clearing seen each call.
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
        # BFS
        head = 1
        tail = 1
        queue[1] = start
        seen[start] = s
        sz = 0
        while head <= tail
            v = queue[head]
            head += 1
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

# Fill present_flag and present_idx for a species i using:
#   present if keep[k] & A[i,k]  (A-only) OR P[i,k] (biotic)
# Returns count.
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

# ----------------------------
# Synthetic environment + abiotic maps
# ----------------------------
function make_environment(rng::AbstractRNG; n1::Int, n2::Int, smooth_iters::Int=25)
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
# Metaweb
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
    for i in 1:S
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

function realized_prey_diversity(mu1::Vector{Float64}, mu2::Vector{Float64},
                                 TL::Vector{Int}, preylist::Vector{Vector{Int}})
    vals = Float64[]
    @inbounds for i in eachindex(TL)
        TL[i] == 1 && continue
        prey = preylist[i]
        m = length(prey)
        m < 2 && continue
        s = 0.0
        np = 0
        for a in 1:m-1
            ja = prey[a]
            for b in a+1:m
                jb = prey[b]
                d1 = mu1[ja]-mu1[jb]; d2 = mu2[ja]-mu2[jb]
                s += sqrt(d1*d1 + d2*d2)
                np += 1
            end
        end
        np > 0 && push!(vals, s/np)
    end
    return isempty(vals) ? NaN : mean(vals)
end

# ----------------------------
# Habitat loss geometries (nested removal)
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
        smooth_field!(Z; iters=30)
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
# Compute AB presence P(i,k) under fixed prey rule: ">= 1 prey present in cell"
# Also returns:
#   - Acount[i], ABcount[i]
#   - cons_present (cells where any consumer is present) and its index list
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

    # A counts (keep & A)
    @inbounds for i in 1:S
        c = 0
        for k in 1:N
            c += (keep[k] & A[i,k]) ? 1 : 0
        end
        Acount[i] = c
    end

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

# ----------------------------
# Patch viability arrays: max patch size for each species (A-only and AB)
# Uses shortcuts: if total count < Epatch, maxpatch = count (no BFS).
# ----------------------------
function compute_maxpatches!(Amax::Vector{Int}, ABmax::Vector{Int},
                            A::BitMatrix, P::BitMatrix, keep::BitVector,
                            Acount::Vector{Int}, ABcount::Vector{Int},
                            Emin_patch::Int,
                            present_flag::BitVector, present_idx::Vector{Int},
                            up::Vector{Int}, down::Vector{Int}, left::Vector{Int}, right::Vector{Int},
                            seen::Vector{Int}, queue::Vector{Int}, stamp::Base.RefValue{Int},
                            N::Int)

    S = length(Acount)

    @inbounds for i in 1:S
        # A-only
        if Acount[i] < Emin_patch
            Amax[i] = Acount[i]
        else
            cnt = fill_present_from_A!(present_flag, present_idx, A, i, keep, N)
            # cnt == Acount[i] but computed here
            Amax[i] = lcc_max_size(present_flag, present_idx, up, down, left, right, seen, queue, stamp)
            clear_present!(present_flag, present_idx)
        end

        # AB
        if ABcount[i] < Emin_patch
            ABmax[i] = ABcount[i]
        else
            cnt = fill_present_from_P!(present_flag, present_idx, P, i, N)
            ABmax[i] = lcc_max_size(present_flag, present_idx, up, down, left, right, seen, queue, stamp)
            clear_present!(present_flag, present_idx)
        end
    end

    return nothing
end

# ----------------------------
# SAR fits from pre-loss A-only (two definitions):
#   - total area presence: count >= Emin_total
#   - patch presence:      maxpatch >= Emin_patch
# Squares are contiguous; we sample many squares per area fraction.
# ----------------------------
function sar_fit_total(rng::AbstractRNG, A::BitMatrix;
                       n1::Int, n2::Int, Emin_total::Int,
                       area_fracs::Vector{Float64}, samples::Int)

    S, N = size(A)
    keep_sq = falses(N)

    a_eff = Float64[]
    Sbar  = Float64[]

    for a in area_fracs
        side = max(2, round(Int, sqrt(a) * min(n1, n2)))
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
                       n1::Int, n2::Int, Emin_patch::Int,
                       area_fracs::Vector{Float64}, samples::Int,
                       up::Vector{Int}, down::Vector{Int}, left::Vector{Int}, right::Vector{Int},
                       present_flag::BitVector, present_idx::Vector{Int},
                       seen::Vector{Int}, queue::Vector{Int}, stamp::Base.RefValue{Int})

    S, N = size(A)
    keep_sq = falses(N)

    a_eff = Float64[]
    Sbar  = Float64[]

    for a in area_fracs
        side = max(2, round(Int, sqrt(a) * min(n1, n2)))
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
                # Build present mask inside square: keep_sq & A[i,:]
                empty!(present_idx)
                for k in 1:N
                    if keep_sq[k] & A[i,k]
                        present_flag[k] = true
                        push!(present_idx, k)
                    end
                end
                if length(present_idx) < Emin_patch
                    # fail
                else
                    maxsz = lcc_max_size(present_flag, present_idx, up, down, left, right, seen, queue, stamp)
                    cntS += (maxsz >= Emin_patch)
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
# MVP runner
# ----------------------------
struct MPParams_patch
    n1::Int; n2::Int; S::Int
    basal_frac::Float64; Lmax::Int
    niche_sigma::Float64; niche_cut::Float64
    k_prey::Int; select_sigma::Float64; match_sigma::Float64
    Emin_total::Int
    Emin_patch::Int
    reps::Int
    fgrid::Vector{Float64}
    geoms::Vector{Symbol}
    sar_area_fracs::Vector{Float64}
    sar_samples::Int
end

function run_mvp_vnext(; seed::Int=1234,
    n1::Int=100, n2::Int=100, S::Int=1000,
    basal_frac::Float64=0.25, Lmax::Int=4,
    niche_sigma::Float64=0.8, niche_cut::Float64=0.45,
    k_prey::Int=3, select_sigma::Float64=0.6, match_sigma::Float64=1.0,
    Emin_total::Int=50,
    Emin_patch::Int=50,  # patch threshold (use same as Emin_total by default)
    reps::Int=20,
    fgrid = collect(0.0:0.05:0.90),
    geoms = [:random, :cluster, :front],
    sar_area_fracs = [0.25, 0.35, 0.50, 0.70, 0.85, 1.0],
    sar_samples::Int=20
)
    p = MPParams_patch(n1,n2,S,basal_frac,Lmax,niche_sigma,niche_cut,k_prey,select_sigma,match_sigma,
                  Emin_total, Emin_patch, reps,
                  Vector{Float64}(fgrid), Vector{Symbol}(geoms),
                  Vector{Float64}(sar_area_fracs), sar_samples)

    rng = MersenneTwister(seed)
    N = p.n1*p.n2
    T = length(p.fgrid)

    # Neighbors + BFS buffers
    up, down, left, right = build_neighbors(p.n1, p.n2)
    present_flag = falses(N)
    present_idx  = Int[]
    cons_present = falses(N)
    cons_idx     = Int[]
    seen  = zeros(Int, N)
    queue = Vector{Int}(undef, N)
    stamp = Ref(0)

    # Aggregates per geometry, per f
    EA_total_inc  = Dict{Symbol, Vector{Float64}}()
    EAB_total_inc = Dict{Symbol, Vector{Float64}}()
    EA_patch_inc  = Dict{Symbol, Vector{Float64}}()
    EAB_patch_inc = Dict{Symbol, Vector{Float64}}()

    phi_cons_mean = Dict{Symbol, Vector{Float64}}()
    lcc_cons_frac = Dict{Symbol, Vector{Float64}}()
    fragfail_rate = Dict{Symbol, Vector{Float64}}()  # consumers: survive total-area but fail patch

    for g in p.geoms
        EA_total_inc[g]  = zeros(T)
        EAB_total_inc[g] = zeros(T)
        EA_patch_inc[g]  = zeros(T)
        EAB_patch_inc[g] = zeros(T)
        phi_cons_mean[g] = zeros(T)
        lcc_cons_frac[g] = zeros(T)
        fragfail_rate[g] = zeros(T)
    end

    # SAR fit lists
    c_tot_list = Float64[]; z_tot_list = Float64[]; r2_tot_list = Float64[]
    c_pat_list = Float64[]; z_pat_list = Float64[]; r2_pat_list = Float64[]

    # Baseline richness lists (A-only and AB) under both criteria
    S0A_tot_list  = Float64[]; S0AB_tot_list = Float64[]
    S0A_pat_list  = Float64[]; S0AB_pat_list = Float64[]
    Dprey_list = Float64[]

    # Working arrays
    Acount  = zeros(Int, p.S)
    ABcount = zeros(Int, p.S)
    Amax    = zeros(Int, p.S)
    ABmax   = zeros(Int, p.S)
    P = falses(p.S, N)  # BitMatrix for AB presence

    for rep in 1:p.reps
        rrng = MersenneTwister(rand(rng, UInt))

        env1, env2 = make_environment(rrng; n1=p.n1, n2=p.n2)
        A, mu1, mu2 = make_abiotic_maps(rrng, env1, env2; S=p.S, niche_sigma=p.niche_sigma, niche_cut=p.niche_cut)
        TL, preylist = build_metaweb(rrng; S=p.S, basal_frac=p.basal_frac, Lmax=p.Lmax,
                                     k_prey=p.k_prey, match_sigma=p.match_sigma, select_sigma=p.select_sigma,
                                     mu1=mu1, mu2=mu2)
        push!(Dprey_list, realized_prey_diversity(mu1, mu2, TL, preylist))

        # SAR fits (pre-loss A-only), both criteria
        cT, zT, r2T = sar_fit_total(rrng, A; n1=p.n1, n2=p.n2, Emin_total=p.Emin_total,
                                   area_fracs=p.sar_area_fracs, samples=p.sar_samples)
        isfinite(cT) && isfinite(zT) && (push!(c_tot_list, cT); push!(z_tot_list, zT); push!(r2_tot_list, r2T))

        cP, zP, r2P = sar_fit_patch(rrng, A; n1=p.n1, n2=p.n2, Emin_patch=p.Emin_patch,
                                   area_fracs=p.sar_area_fracs, samples=p.sar_samples,
                                   up=up, down=down, left=left, right=right,
                                   present_flag=present_flag, present_idx=present_idx,
                                   seen=seen, queue=queue, stamp=stamp)
        isfinite(cP) && isfinite(zP) && (push!(c_pat_list, cP); push!(z_pat_list, zP); push!(r2_pat_list, r2P))

        # Precompute nested loss orders per geometry
        ord = Dict{Symbol, Vector{Int}}()
        for g in p.geoms
            ord[g] = make_loss_order(rrng, g, env1)
        end

        # Baselines at f=0 (same keep0 for all geoms)
        keep0 = trues(N)
        compute_P_and_counts!(P, Acount, ABcount, A, keep0, TL, preylist, cons_present, cons_idx, N)
        compute_maxpatches!(Amax, ABmax, A, P, keep0, Acount, ABcount, p.Emin_patch,
                            present_flag, present_idx, up, down, left, right, seen, queue, stamp, N)

        S0A_tot  = count(>=(p.Emin_total), Acount)
        S0AB_tot = count(>=(p.Emin_total), ABcount)
        S0A_pat  = count(>=(p.Emin_patch), Amax)
        S0AB_pat = count(>=(p.Emin_patch), ABmax)

        push!(S0A_tot_list, S0A_tot);   push!(S0AB_tot_list, S0AB_tot)
        push!(S0A_pat_list, S0A_pat);   push!(S0AB_pat_list, S0AB_pat)

        # Sweep f for each geometry
        for g in p.geoms
            for (t, f) in enumerate(p.fgrid)
                keep = (f == 0.0) ? keep0 : keep_from_order(ord[g], f, N)

                compute_P_and_counts!(P, Acount, ABcount, A, keep, TL, preylist, cons_present, cons_idx, N)

                # LCC fraction on consumer-supported overlap footprint
                tot_cons_cells = length(cons_idx)
                if tot_cons_cells == 0
                    lcc = 0.0
                else
                    maxsz = lcc_max_size(cons_present, cons_idx, up, down, left, right, seen, queue, stamp)
                    lcc = maxsz / tot_cons_cells
                end
                lcc_cons_frac[g][t] += lcc

                # Max patch sizes per species (A-only and AB)
                compute_maxpatches!(Amax, ABmax, A, P, keep, Acount, ABcount, p.Emin_patch,
                                    present_flag, present_idx, up, down, left, right, seen, queue, stamp, N)

                SA_tot  = count(>=(p.Emin_total), Acount)
                SAB_tot = count(>=(p.Emin_total), ABcount)
                SA_pat  = count(>=(p.Emin_patch), Amax)
                SAB_pat = count(>=(p.Emin_patch), ABmax)

                EA_total_inc[g][t]  += (S0A_tot  - SA_tot)
                EAB_total_inc[g][t] += (S0AB_tot - SAB_tot)
                EA_patch_inc[g][t]  += (S0A_pat  - SA_pat)
                EAB_patch_inc[g][t] += (S0AB_pat - SAB_pat)

                # mean φ among consumers (TL>1, Acount>0)
                sφ = 0.0; nφ = 0
                @inbounds for i in 1:p.S
                    if TL[i] > 1 && Acount[i] > 0
                        sφ += ABcount[i] / Acount[i]
                        nφ += 1
                    end
                end
                phi_cons_mean[g][t] += (nφ > 0 ? sφ/nφ : NaN)

                # fragmentation-failure rate among consumers:
                # survive by total-area (ABcount>=Emin_total) but fail patch (ABmax<Emin_patch)
                ff = 0; nn = 0
                @inbounds for i in 1:p.S
                    if TL[i] > 1
                        nn += 1
                        if (ABcount[i] >= p.Emin_total) && (ABmax[i] < p.Emin_patch)
                            ff += 1
                        end
                    end
                end
                fragfail_rate[g][t] += (nn > 0 ? ff/nn : NaN)
            end
        end
    end

    # Average over reps
    for g in p.geoms
        EA_total_inc[g]  ./= p.reps
        EAB_total_inc[g] ./= p.reps
        EA_patch_inc[g]  ./= p.reps
        EAB_patch_inc[g] ./= p.reps
        phi_cons_mean[g] ./= p.reps
        lcc_cons_frac[g] ./= p.reps
        fragfail_rate[g] ./= p.reps
    end

    # SAR prediction curves for both criteria (area-only)
    S0A_tot_m = mean(S0A_tot_list)
    S0A_pat_m = mean(S0A_pat_list)

    cT = mean(c_tot_list); zT = mean(z_tot_list); r2T = mean(r2_tot_list)
    cP = mean(c_pat_list); zP = mean(z_pat_list); r2P = mean(r2_pat_list)

    SAR_total_pred = zeros(T)
    SAR_patch_pred = zeros(T)
    @inbounds for (t, f) in enumerate(p.fgrid)
        h = max(1e-9, 1.0 - f)
        SpT = clamp(cT * h^zT, 0.0, S0A_tot_m)
        SpP = clamp(cP * h^zP, 0.0, S0A_pat_m)
        SAR_total_pred[t] = S0A_tot_m - SpT
        SAR_patch_pred[t] = S0A_pat_m - SpP
    end

    # ----------------------------
    # Plot
    # ----------------------------
    fig = Figure(size=(1750, 1250))

    # Row 1: total-area criterion
    for (j, g) in enumerate(p.geoms)
        ax = Axis(fig[1, j],
            xlabel = "habitat loss f",
            ylabel = (j == 1 ? "incremental extinctions" : ""),
            title  = "TOTAL-AREA criterion (Emin=$(p.Emin_total))   geometry=$(String(g))"
        )
        lines!(ax, p.fgrid, EA_total_inc[g];  label="A-only obs", linewidth=3)
        lines!(ax, p.fgrid, EAB_total_inc[g]; label="AB obs (biotic)", linewidth=3)
        lines!(ax, p.fgrid, SAR_total_pred;   label="SAR baseline (area-only)", linestyle=:dash, linewidth=3)
        axislegend(ax; position=:lt)
    end

    # Row 2: patch-viability criterion
    for (j, g) in enumerate(p.geoms)
        ax = Axis(fig[2, j],
            xlabel = "habitat loss f",
            ylabel = (j == 1 ? "incremental extinctions" : ""),
            title  = "PATCH criterion (Epatch=$(p.Emin_patch))   geometry=$(String(g))"
        )
        lines!(ax, p.fgrid, EA_patch_inc[g];  label="A-only obs", linewidth=3)
        lines!(ax, p.fgrid, EAB_patch_inc[g]; label="AB obs (biotic)", linewidth=3)
        lines!(ax, p.fgrid, SAR_patch_pred;   label="SAR baseline (area-only)", linestyle=:dash, linewidth=3)
        axislegend(ax; position=:lt)
    end

    # Row 3: Mechanisms over f
    #   left axis: mean φ(consumers)
    #   right axis: LCC fraction (consumer overlap) and fragmentation-failure rate
    for (j, g) in enumerate(p.geoms)
        axL = Axis(fig[3, j],
            xlabel = "habitat loss f",
            ylabel = (j == 1 ? "mean φ (consumers)" : ""),
            title = "mechanisms (same AB maps; different survival criteria)"
        )
        axR = Axis(fig[3, j],
            yaxisposition = :right,
            ylabel = (j == 3 ? "LCC / frag-fail rate" : ""),
            xgridvisible = false, ygridvisible = false
        )
        linkxaxes!(axL, axR)
        hidespines!(axR)
        hidexdecorations!(axR)

        ylims!(axL, (0.0, 1.0))
        ylims!(axR, (0.0, 1.0))

        lines!(axL, p.fgrid, phi_cons_mean[g]; linewidth=3, label="mean φ")
        lines!(axR, p.fgrid, lcc_cons_frac[g]; linewidth=3, linestyle=:dash, label="LCC fraction")
        lines!(axR, p.fgrid, fragfail_rate[g]; linewidth=3, linestyle=:dot, label="frag-fail rate")

        axislegend(axL; position=:rb)
    end

    # Row 4: Mechanistic summary scatter: SAR error vs LCC (both criteria)
    axS = Axis(fig[4, 1:3],
        xlabel = "LCC fraction of consumer-supported overlap",
        ylabel = "SAR error vs biotic reality (AB_inc − SAR_pred)",
        title  = "mechanistic link: when overlap fragments, SAR error changes (two criteria shown)"
    )
    for g in p.geoms
        # points across f
        x = lcc_cons_frac[g]
        yT = (EAB_total_inc[g] .- SAR_total_pred)
        yP = (EAB_patch_inc[g] .- SAR_patch_pred)

        scatter!(axS, x, yT; markersize=10, label="$(String(g)) total")
        scatter!(axS, x, yP; markersize=10, marker=:utriangle, label="$(String(g)) patch")
        lines!(axS, x, yT; linewidth=2)
        lines!(axS, x, yP; linewidth=2, linestyle=:dash)
    end
    axislegend(axS; position=:rb, nbanks=2)

    display(fig)

    # ----------------------------
    # Print summary
    # ----------------------------
    @printf("\n===== vNext SUMMARY =====\n")
    @printf("grid=%dx%d, S=%d, reps=%d\n", p.n1, p.n2, p.S, p.reps)
    @printf("metaweb: k_prey=%d, select_sigma=%.2f, match_sigma=%.2f, basal_frac=%.2f, Lmax=%d\n",
            p.k_prey, p.select_sigma, p.match_sigma, p.basal_frac, p.Lmax)
    @printf("criteria: Emin_total=%d, Epatch=%d\n", p.Emin_total, p.Emin_patch)
    @printf("mean realized prey diversity Dprey = %.4f\n", mean(Dprey_list))

    @printf("\nBaseline richness at f=0 (means over reps):\n")
    @printf("  TOTAL:  S_A0=%.2f,  S_AB0=%.2f\n", mean(S0A_tot_list), mean(S0AB_tot_list))
    @printf("  PATCH:  S_A0=%.2f,  S_AB0=%.2f\n", mean(S0A_pat_list), mean(S0AB_pat_list))

    @printf("\nSAR fits (pre-loss A-only; means over reps):\n")
    @printf("  TOTAL: z=%.3f, c=%.2f, R2(log)=%.3f\n", zT, cT, r2T)
    @printf("  PATCH: z=%.3f, c=%.2f, R2(log)=%.3f\n", zP, cP, r2P)

    @printf("\nAUC of biotic amplification dE* (AB_inc - A_inc):\n")
    for g in p.geoms
        dE_T = EAB_total_inc[g] .- EA_total_inc[g]
        dE_P = EAB_patch_inc[g] .- EA_patch_inc[g]
        @printf("  %s:  total=%.3f   patch=%.3f\n", String(g),
                auc_trapz(p.fgrid, dE_T), auc_trapz(p.fgrid, dE_P))
    end

    @printf("\nAUC of SAR error (AB_inc - SAR_pred):\n")
    for g in p.geoms
        errT = EAB_total_inc[g] .- SAR_total_pred
        errP = EAB_patch_inc[g] .- SAR_patch_pred
        @printf("  %s:  total=%.3f   patch=%.3f\n", String(g),
                auc_trapz(p.fgrid, errT), auc_trapz(p.fgrid, errP))
    end
    println("")

    return nothing
end

# ----------------------------
# Run
# ----------------------------
run_mvp_vnext(
    seed=1234,
    reps=20,
    k_prey=3,
    select_sigma=0.6,
    Emin_total=50,
    Emin_patch=50
)
