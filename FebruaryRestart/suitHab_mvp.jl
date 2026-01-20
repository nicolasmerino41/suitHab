# suitHab_mvp.jl
#
# Synthetic MVP:
#   Compare extinctions predicted by abiotic habitat only (A) vs effective habitat requiring trophic support (AB).
#   Knobs: diet redundancy (k_prey), prey niche synchrony (select_sigma), habitat-loss geometry.

using Random
using Statistics
using Printf
using DelimitedFiles

# ----------------------------
# Utilities
# ----------------------------

# Simple spatial smoothing to create autocorrelated random fields (no external deps)
function smooth_field!(Z::Matrix{Float64}; iters::Int=20)
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

# Normalize to mean 0, std 1
function zscore!(Z::Matrix{Float64})
    μ = mean(Z); σ = std(Z)
    Z .-= μ
    Z ./= (σ > 0 ? σ : 1.0)
    return Z
end

# Convert 2D -> 1D cell index
@inline idx(i, j, n2) = (i-1)*n2 + j

# ----------------------------
# Synthetic landscape + niches
# ----------------------------

"""
make_environment(rng; n1, n2, smooth_iters)

Returns two environmental layers env1, env2 (Float64 matrices), autocorrelated and z-scored.
"""
function make_environment(rng::AbstractRNG; n1::Int=80, n2::Int=80, smooth_iters::Int=25)
    env1 = randn(rng, n1, n2)
    env2 = randn(rng, n1, n2)
    smooth_field!(env1; iters=smooth_iters)
    smooth_field!(env2; iters=smooth_iters)
    zscore!(env1); zscore!(env2)
    return env1, env2
end

"""
make_abiotic_maps(rng, env1, env2; S, niche_sigma, niche_cut)

Each species i has an optimum (mu1[i], mu2[i]) drawn from N(0,1).
Abiotic suitability A_i(x) = exp(-d^2/(2*sigma^2)) >= niche_cut.
Returns A as BitMatrix (S x Ncells) and the optima arrays.
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
# Metaweb generation (acyclic trophic levels)
# ----------------------------

"""
build_metaweb(rng; S, basal_frac, Lmax, k_prey, select_sigma, mu1, mu2)

- Assign trophic level TL[i] in 1..Lmax (basals forced TL=1).
- For each consumer i, choose k_prey prey from lower trophic levels.
- Prey choice biased by niche similarity with width select_sigma:
    weight(j) = exp(-||mu(i)-mu(j)||^2 / (2*select_sigma^2))
Smaller select_sigma => prey niches more similar => higher "prey synchrony".
Returns basals::Vector{Int}, TL::Vector{Int}, preylist::Vector{Vector{Int}}
"""
function build_metaweb(rng::AbstractRNG;
                       S::Int=80,
                       basal_frac::Float64=0.25,
                       Lmax::Int=4,
                       k_prey::Int=3,
                       select_sigma::Float64=0.6,
                       mu1::Vector{Float64},
                       mu2::Vector{Float64})

    nb = max(1, round(Int, basal_frac * S))
    perm = randperm(rng, S)
    basals = perm[1:nb]
    is_basal = falses(S)
    is_basal[basals] .= true

    TL = fill(1, S)
    # Assign non-basals random TL in 2..Lmax
    for i in 1:S
        if !is_basal[i]
            TL[i] = rand(rng, 2:Lmax)
        end
    end

    preylist = [Int[] for _ in 1:S]

    # Helper: sample without replacement from candidates with weights
    function sample_weighted_no_replace(cands::Vector{Int}, w::Vector{Float64}, k::Int)
        k = min(k, length(cands))
        chosen = Int[]
        if k == 0
            return chosen
        end
        # simple sequential sampling: each draw reweights remaining
        remaining = collect(1:length(cands))
        for _ in 1:k
            ww = w[remaining]
            s = sum(ww)
            if s <= 0
                # fallback uniform
                pick_idx = rand(rng, remaining)
            else
                u = rand(rng) * s
                acc = 0.0
                pick_idx = remaining[end]
                for ridx in remaining
                    acc += w[ridx]
                    if acc >= u
                        pick_idx = ridx
                        break
                    end
                end
            end
            push!(chosen, cands[pick_idx])
            deleteat!(remaining, findfirst(==(pick_idx), remaining))
            if isempty(remaining)
                break
            end
        end
        return chosen
    end

    for i in 1:S
        if TL[i] == 1
            continue
        end
        # candidates: species with lower TL
        cands = [j for j in 1:S if TL[j] < TL[i]]
        if isempty(cands)
            # ensure at least one basal exists
            cands = collect(basals)
        end
        # weights by niche similarity
        w = Vector{Float64}(undef, length(cands))
        @inbounds for (t, j) in enumerate(cands)
            d1 = mu1[i] - mu1[j]
            d2 = mu2[i] - mu2[j]
            d2sum = d1*d1 + d2*d2
            w[t] = exp(-d2sum / (2.0 * select_sigma^2))
        end
        prey = sample_weighted_no_replace(cands, w, k_prey)
        preylist[i] = prey
    end

    return basals, TL, preylist
end

# ----------------------------
# Habitat loss operators
# ----------------------------

"""
loss_mask(rng; geometry, f, n1, n2, env1)

Returns keep::BitVector (length Ncells), true = cell remains.
Geometries:
- :random    -> uniform random removal
- :cluster   -> remove top quantile of a smoothed random field (spatially clustered)
- :front     -> remove cells with env1 above a threshold (climate-change-like front)
"""
function loss_mask(rng::AbstractRNG; geometry::Symbol, f::Float64, n1::Int, n2::Int, env1::Matrix{Float64})
    N = n1 * n2
    f = clamp(f, 0.0, 1.0)
    if f == 0.0
        return trues(N)
    elseif f == 1.0
        return falses(N)
    end

    score = Vector{Float64}(undef, N)

    if geometry == :random
        @inbounds for k in 1:N
            score[k] = rand(rng)
        end
    elseif geometry == :cluster
        Z = randn(rng, n1, n2)
        smooth_field!(Z; iters=30)
        zscore!(Z)
        @inbounds for r in 1:n1, c in 1:n2
            score[idx(r,c,n2)] = Z[r,c]
        end
    elseif geometry == :front
        @inbounds for r in 1:n1, c in 1:n2
            score[idx(r,c,n2)] = env1[r,c]  # remove "high env1" first
        end
    else
        error("Unknown geometry: $geometry")
    end

    # Remove the top f fraction (highest scores)
    thr = quantile(score, 1.0 - f)
    keep = BitVector(undef, N)
    @inbounds for k in 1:N
        keep[k] = (score[k] <= thr)
    end
    return keep
end

# Simple quantile without deps (approx via sort)
function quantile(v::Vector{Float64}, p::Float64)
    vv = sort(copy(v))
    n = length(vv)
    p = clamp(p, 0.0, 1.0)
    if p == 0.0
        return vv[1]
    elseif p == 1.0
        return vv[end]
    end
    pos = 1 + (n-1)*p
    lo = floor(Int, pos)
    hi = ceil(Int, pos)
    if lo == hi
        return vv[lo]
    else
        w = pos - lo
        return (1-w)*vv[lo] + w*vv[hi]
    end
end

# ----------------------------
# Biotic assembly (trophic support)
# ----------------------------

"""
assemble_presence(A, keep, TL, preylist; B_mode=:strong, beta=1.0, thetaB=0.2)

Returns P::BitMatrix same size as A (S x Ncells) where:
- Basals: P = A masked by keep
- Consumers:
    - :none   -> P = A masked by keep
    - :strong -> require at least one prey present in cell
    - :soft   -> require (mean prey present)^beta >= thetaB
"""
function assemble_presence(A::BitMatrix, keep::BitVector,
                           TL::Vector{Int}, preylist::Vector{Vector{Int}};
                           B_mode::Symbol=:strong, beta::Float64=1.0, thetaB::Float64=0.2)

    S, N = size(A)

    # Abiotic masked by habitat loss
    A2 = copy(A)
    @inbounds for k in 1:N
        if !keep[k]
            for i in 1:S
                A2[i,k] = false
            end
        end
    end

    if B_mode == :none
        return A2
    end

    # Process species in ascending trophic level
    order = sortperm(TL)
    P = BitMatrix(undef, S, N)
    P .= false

    for i in order
        if TL[i] == 1
            @inbounds for k in 1:N
                P[i,k] = A2[i,k]
            end
        else
            prey = preylist[i]
            if isempty(prey)
                # If no prey defined, treat as basal (failsafe)
                @inbounds for k in 1:N
                    P[i,k] = A2[i,k]
                end
                continue
            end

            if B_mode == :strong
                # support[k] = OR over prey P[prey,k]
                @inbounds for k in 1:N
                    sup = false
                    for pj in prey
                        sup |= P[pj,k]
                        if sup
                            break
                        end
                    end
                    P[i,k] = A2[i,k] & sup
                end

            elseif B_mode == :soft
                m = length(prey)
                @inbounds for k in 1:N
                    cnt = 0
                    for pj in prey
                        cnt += (P[pj,k] ? 1 : 0)
                    end
                    frac = cnt / m
                    sup = (frac^beta >= thetaB)
                    P[i,k] = A2[i,k] & sup
                end
            else
                error("Unknown B_mode: $B_mode")
            end
        end
    end
    return P
end

# Extinction counts
function extinction_count(P::BitMatrix)
    S, N = size(P)
    ext = 0
    @inbounds for i in 1:S
        alive = false
        for k in 1:N
            alive |= P[i,k]
            if alive
                break
            end
        end
        ext += (!alive)
    end
    return ext
end

function extinction_count_by_TL(P::BitMatrix, TL::Vector{Int}, Lmax::Int)
    S, N = size(P)
    extL = zeros(Int, Lmax)
    totL = zeros(Int, Lmax)
    @inbounds for i in 1:S
        l = TL[i]
        totL[l] += 1
        alive = false
        for k in 1:N
            alive |= P[i,k]
            if alive
                break
            end
        end
        extL[l] += (!alive)
    end
    return extL, totL
end

# ----------------------------
# Experiment runner
# ----------------------------
struct RunParams
    n1::Int
    n2::Int
    S::Int
    basal_frac::Float64
    Lmax::Int
    niche_sigma::Float64
    niche_cut::Float64
    k_prey::Int
    select_sigma::Float64
    geometry::Symbol
    B_mode::Symbol
    beta::Float64
    thetaB::Float64
end

function run_one(rng::AbstractRNG, p::RunParams, fgrid::Vector{Float64})
    env1, env2 = make_environment(rng; n1=p.n1, n2=p.n2)
    A, mu1, mu2 = make_abiotic_maps(rng, env1, env2; S=p.S, niche_sigma=p.niche_sigma, niche_cut=p.niche_cut)
    basals, TL, preylist = build_metaweb(rng; S=p.S, basal_frac=p.basal_frac, Lmax=p.Lmax,
                                         k_prey=p.k_prey, select_sigma=p.select_sigma, mu1=mu1, mu2=mu2)

    out = Vector{NamedTuple}(undef, length(fgrid))

    for (t, f) in enumerate(fgrid)
        keep = loss_mask(rng; geometry=p.geometry, f=f, n1=p.n1, n2=p.n2, env1=env1)

        # A-only extinctions: just abiotic masked
        A_only = assemble_presence(A, keep, TL, preylist; B_mode=:none)
        EA = extinction_count(A_only)

        # AB extinctions: abiotic + trophic support
        AB = assemble_presence(A, keep, TL, preylist; B_mode=p.B_mode, beta=p.beta, thetaB=p.thetaB)
        EAB = extinction_count(AB)

        out[t] = (
            f = f,
            EA = EA,
            EAB = EAB,
            dE = EAB - EA
        )
    end
    return out
end

function run_sweep(; seed::Int=1234, reps::Int=20)
    rng = MersenneTwister(seed)

    # Default habitat loss fractions
    fgrid = collect(0.0:0.05:0.90)

    # Sweep dimensions
    kprey_list   = [1, 3, 6]
    selects_list = [0.25, 0.6, 1.5]          # smaller => higher prey synchrony
    geoms        = [:random, :cluster, :front]

    # Fixed "ecology"
    B_mode = :strong

    # Landscape + community defaults
    base = RunParams(
        80, 80,      # grid
        80,          # species
        0.25,        # basal fraction
        4,           # Lmax
        0.9,         # niche_sigma
        0.35,        # niche_cut
        3,           # k_prey (overwritten)
        0.6,         # select_sigma (overwritten)
        :random,     # geometry (overwritten)
        B_mode,
        1.0,         # beta
        0.2          # thetaB
    )

    rows = String[]
    push!(rows, "rep,geometry,k_prey,select_sigma,f,EA,EAB,dE")

    for g in geoms, kp in kprey_list, ss in selects_list
        p = RunParams(base.n1, base.n2, base.S, base.basal_frac, base.Lmax,
                      base.niche_sigma, base.niche_cut,
                      kp, ss, g, base.B_mode, base.beta, base.thetaB)

        EA_sum  = zeros(Float64, length(fgrid))
        EAB_sum = zeros(Float64, length(fgrid))
        dE_sum  = zeros(Float64, length(fgrid))

        for r in 1:reps
            rrng = MersenneTwister(rand(rng, UInt))
            out = run_one(rrng, p, fgrid)

            for (t, nt) in enumerate(out)
                EA_sum[t]  += nt.EA
                EAB_sum[t] += nt.EAB
                dE_sum[t]  += nt.dE

                push!(rows, @sprintf("%d,%s,%d,%.3f,%.2f,%.1f,%.1f,%.1f",
                                     r, String(g), kp, ss, nt.f, nt.EA, nt.EAB, nt.dE))
            end
        end

        @printf("\nDONE geometry=%s k_prey=%d select_sigma=%.3f\n", String(g), kp, ss)
        t50 = findfirst(==(0.50), fgrid)
        @printf("  mean dE at f=0.50: %.3f\n", dE_sum[t50] / reps)
    end

    # ---------- plotting (no @eval, no scope issues) ----------
    fig = Figure(size=(900, 600))
    ax = Axis(fig[1,1],
        xlabel="habitat loss f",
        ylabel="biotic amplification dE = EAB − EA",
        title="Synthetic MVP (geometry = random)"
    )

    for kp in kprey_list, ss in selects_list
        dEmean = zeros(Float64, length(fgrid))
        count  = zeros(Int, length(fgrid))

        for line in rows[2:end]
            parts = split(line, ',')
            geom = Symbol(parts[2])
            kpp  = parse(Int, parts[3])
            sss  = parse(Float64, parts[4])

            if geom == :random && kpp == kp && abs(sss - ss) < 1e-9
                f = parse(Float64, parts[5])
                t = findfirst(==(f), fgrid)
                dEmean[t] += parse(Float64, parts[8])
                count[t]  += 1
            end
        end

        dEmean ./= max.(count, 1)
        lines!(ax, fgrid, dEmean;
                label=@sprintf("k=%d, sel=%.2f", kp, ss))
    end

    axislegend(ax; position=:rb)
    display(fig)

    return nothing
end

# ----------------------------
# Main
# ----------------------------
run_sweep(seed=1234, reps=15)
