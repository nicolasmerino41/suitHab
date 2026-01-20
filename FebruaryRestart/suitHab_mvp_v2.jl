# suitHab_mvp_v2.jl
#
# Synthetic MVP v2:
#   Compare abiotic-only extinctions (A) vs effective-habitat extinctions requiring trophic support (AB).
#   Improvements:
#     - Extinction threshold Emin (minimum remaining suitable cells)
#     - Baseline-corrected amplification ΔE*(f) to remove f=0 mismatch artifacts
#     - Prey-guild metaweb: select_sigma controls prey-prey similarity ("synchrony")
#     - One loss order per run (nested removals), fast
#     - Plot is displayed (no CSV, no saving)

using Random
using Statistics
using Printf

# Plotting
using CairoMakie

# ----------------------------
# Utilities: simple spatial smoothing (no external deps)
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

@inline idx(i::Int, j::Int, n2::Int) = (i - 1) * n2 + j

# Weighted sampling helpers (no macros, no deps)
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
    if k == 0
        return chosen
    end
    remaining = collect(1:length(cands))
    for _ in 1:k
        ww = w[remaining]
        pick_local = weighted_pick_index(rng, ww)
        pick_idx = remaining[pick_local]
        push!(chosen, cands[pick_idx])
        deleteat!(remaining, pick_local)
        if isempty(remaining)
            break
        end
    end
    return chosen
end

# ----------------------------
# Synthetic landscape + abiotic niches
# ----------------------------

function make_environment(rng::AbstractRNG; n1::Int=80, n2::Int=80, smooth_iters::Int=25)
    env1 = randn(rng, n1, n2)
    env2 = randn(rng, n1, n2)
    smooth_field!(env1; iters=smooth_iters)
    smooth_field!(env2; iters=smooth_iters)
    zscore!(env1)
    zscore!(env2)
    return env1, env2
end

"""
Abiotic suitability maps A[i,cell] as thresholded Gaussian niche in 2D env space.
Returns BitMatrix A (S x Ncells) and species optima mu1,mu2.
"""
function make_abiotic_maps(rng::AbstractRNG, env1::Matrix{Float64}, env2::Matrix{Float64};
                           S::Int=80, niche_sigma::Float64=0.9, niche_cut::Float64=0.35)
    n1, n2 = size(env1)
    N = n1 * n2

    mu1 = randn(rng, S)
    mu2 = randn(rng, S)

    A = BitMatrix(undef, S, N)
    @inbounds for i in 1:S
        for r in 1:n1, c in 1:n2
            d1 = env1[r,c] - mu1[i]
            d2 = env2[r,c] - mu2[i]
            d2sum = d1*d1 + d2*d2
            suit = exp(-d2sum / (2.0 * niche_sigma^2))
            A[i, idx(r,c,n2)] = (suit >= niche_cut)
        end
    end
    return A, mu1, mu2
end

# ----------------------------
# Metaweb generation (acyclic trophic levels) with prey-guild synchrony
# ----------------------------

"""
build_metaweb with prey-guild construction:
- Consumers pick prey from lower trophic levels.
- First pick a "guild seed" prey g close to consumer niche (match_sigma controls feasibility).
- Then pick remaining prey clustered around g (select_sigma controls prey-prey similarity / synchrony).

Interpretation:
- smaller select_sigma => prey are more similar => higher synchrony => redundancy less effective under habitat loss
"""
function build_metaweb(rng::AbstractRNG;
                       S::Int=80,
                       basal_frac::Float64=0.25,
                       Lmax::Int=4,
                       k_prey::Int=3,
                       match_sigma::Float64=0.8,     # consumer-prey feasibility scale
                       select_sigma::Float64=0.6,    # prey-prey clustering scale (synchrony knob)
                       mu1::Vector{Float64},
                       mu2::Vector{Float64})

    nb = max(1, round(Int, basal_frac * S))
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
        if TL[i] == 1
            continue
        end

        # Candidate prey: lower trophic level
        cands = [j for j in 1:S if TL[j] < TL[i]]
        if isempty(cands)
            cands = collect(basals)
        end

        # Weights for consumer-prey niche match (feasibility)
        wmatch = Vector{Float64}(undef, length(cands))
        for (t, j) in enumerate(cands)
            d1 = mu1[i] - mu1[j]
            d2 = mu2[i] - mu2[j]
            d2sum = d1*d1 + d2*d2
            wmatch[t] = exp(-d2sum / (2.0 * match_sigma^2))
        end

        # Pick guild seed g according to match weights
        t_g = weighted_pick_index(rng, wmatch)
        g = cands[t_g]

        if k_prey == 1
            preylist[i] = [g]
            continue
        end

        # For remaining prey, weight by BOTH:
        #   - closeness to guild seed (prey-prey similarity; select_sigma)
        #   - closeness to consumer (feasibility; match_sigma)
        # (this keeps prey feasible while letting select_sigma control synchrony)
        cands2 = [j for j in cands if j != g]
        if isempty(cands2)
            preylist[i] = [g]
            continue
        end

        w = Vector{Float64}(undef, length(cands2))
        for (t, j) in enumerate(cands2)
            d1g = mu1[g] - mu1[j]
            d2g = mu2[g] - mu2[j]
            d2g_sum = d1g*d1g + d2g*d2g
            wguild = exp(-d2g_sum / (2.0 * select_sigma^2))

            d1i = mu1[i] - mu1[j]
            d2i = mu2[i] - mu2[j]
            d2i_sum = d1i*d1i + d2i*d2i
            wfeas = exp(-d2i_sum / (2.0 * match_sigma^2))

            w[t] = wguild * wfeas
        end

        others = sample_weighted_no_replace(rng, cands2, w, k_prey - 1)
        preylist[i] = vcat([g], others)
    end

    return basals, TL, preylist
end

# ----------------------------
# Habitat loss: one removal order per run (nested removals)
# ----------------------------

function make_loss_order(rng::AbstractRNG, geometry::Symbol, env1::Matrix{Float64})
    n1, n2 = size(env1)
    N = n1 * n2
    score = Vector{Float64}(undef, N)

    if geometry == :random
        # nested random removals without replacement
        return randperm(rng, N)

    elseif geometry == :cluster
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

    # highest score removed first
    ord = sortperm(score; rev=true)
    return ord
end

function keep_from_order(ord::Vector{Int}, f::Float64, N::Int)
    f = clamp(f, 0.0, 1.0)
    keep = trues(N)
    nrem = round(Int, f * N)
    nrem = clamp(nrem, 0, N)
    @inbounds for t in 1:nrem
        keep[ord[t]] = false
    end
    return keep
end

# ----------------------------
# Extinction counts
# ----------------------------

# A-only: ex
