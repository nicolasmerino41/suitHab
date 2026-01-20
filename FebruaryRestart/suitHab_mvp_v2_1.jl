# suitHab_mvp_v2_1.jl
#
# Synthetic MVP v2.1:
#   - 3 geometry panels (:random, :cluster, :front)
#   - overlays mean EA (faint dotted), mean EAB (faint dashed), mean dE* (bold)
#   - prints dE*(0.5) and AUC(dE*) for each (geometry, k_prey, select_sigma)
#   - no CSV, no saving; display(fig) only

using Random
using Statistics
using Printf
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

# Weighted sampling helpers (no deps)
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

# Trapezoid AUC
function auc_trapz(x::Vector{Float64}, y::Vector{Float64})
    s = 0.0
    @inbounds for i in 1:length(x)-1
        dx = x[i+1] - x[i]
        s += dx * (y[i] + y[i+1]) / 2.0
    end
    return s
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
                           S::Int=80, niche_sigma::Float64=0.8, niche_cut::Float64=0.45)
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
# Metaweb (acyclic TL) with prey-guild synchrony
# ----------------------------

"""
build_metaweb with prey-guild construction:
- Consumers pick prey from lower trophic levels.
- Pick a guild seed prey g close to consumer niche (match_sigma controls feasibility).
- Pick remaining prey clustered around g (select_sigma controls prey-prey similarity / synchrony).

Interpretation:
- smaller select_sigma => prey more similar => higher synchrony => redundancy less effective
"""
function build_metaweb(rng::AbstractRNG;
                       S::Int=80,
                       basal_frac::Float64=0.25,
                       Lmax::Int=4,
                       k_prey::Int=3,
                       match_sigma::Float64=1.0,
                       select_sigma::Float64=0.6,
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

        cands = [j for j in 1:S if TL[j] < TL[i]]
        if isempty(cands)
            cands = collect(basals)
        end

        # Feasibility weights consumer -> prey
        wmatch = Vector{Float64}(undef, length(cands))
        for (t, j) in enumerate(cands)
            d1 = mu1[i] - mu1[j]
            d2 = mu2[i] - mu2[j]
            d2sum = d1*d1 + d2*d2
            wmatch[t] = exp(-d2sum / (2.0 * match_sigma^2))
        end

        # Guild seed
        t_g = weighted_pick_index(rng, wmatch)
        g = cands[t_g]

        if k_prey == 1
            preylist[i] = [g]
            continue
        end

        cands2 = [j for j in cands if j != g]
        if isempty(cands2)
            preylist[i] = [g]
            continue
        end

        # Combined weights: prey close to guild seed (synchrony) AND feasible for consumer
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
# Habitat loss: nested removals via one order per run
# ----------------------------

function make_loss_order(rng::AbstractRNG, geometry::Symbol, env1::Matrix{Float64})
    n1, n2 = size(env1)
    N = n1 * n2
    score = Vector{Float64}(undef, N)

    if geometry == :random
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

    return sortperm(score; rev=true)  # highest removed first
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
# Extinction counts (with Emin threshold)
# ----------------------------

function extinction_count_Aonly(A::BitMatrix, keep::BitVector; Emin::Int=1)
    S, N = size(A)
    ext = 0
    @inbounds for i in 1:S
        cnt = 0
        for k in 1:N
            if keep[k] & A[i,k]
                cnt += 1
                if cnt >= Emin
                    break
                end
            end
        end
        ext += (cnt < Emin)
    end
    return ext
end

function extinction_count_AB_strong(A::BitMatrix, keep::BitVector,
                                   TL::Vector{Int}, preylist::Vector{Vector{Int}};
                                   Emin::Int=1)
    S, N = size(A)
    order = sortperm(TL)

    P = BitMatrix(undef, S, N)
    P .= false

    @inbounds for i in order
        if TL[i] == 1
            for k in 1:N
                P[i,k] = keep[k] & A[i,k]
            end
        else
            prey = preylist[i]
            if isempty(prey)
                for k in 1:N
                    P[i,k] = keep[k] & A[i,k]
                end
            else
                for k in 1:N
                    sup = false
                    for pj in prey
                        sup |= P[pj,k]
                        if sup
                            break
                        end
                    end
                    P[i,k] = keep[k] & A[i,k] & sup
                end
            end
        end
    end

    ext = 0
    @inbounds for i in 1:S
        cnt = 0
        for k in 1:N
            if P[i,k]
                cnt += 1
                if cnt >= Emin
                    break
                end
            end
        end
        ext += (cnt < Emin)
    end
    return ext
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
    match_sigma::Float64
    geometry::Symbol
    Emin::Int
end

"""
One replicate (fixed env + web) returns mean curves over fgrid:
EA(f), EAB(f), dE*(f) where:
dE*(f) = (EAB(f)-EAB(0)) - (EA(f)-EA(0))
"""
function run_one(rng::AbstractRNG, p::RunParams, fgrid::Vector{Float64})
    env1, env2 = make_environment(rng; n1=p.n1, n2=p.n2)
    A, mu1, mu2 = make_abiotic_maps(rng, env1, env2; S=p.S, niche_sigma=p.niche_sigma, niche_cut=p.niche_cut)

    _, TL, preylist = build_metaweb(rng;
        S=p.S, basal_frac=p.basal_frac, Lmax=p.Lmax,
        k_prey=p.k_prey, match_sigma=p.match_sigma, select_sigma=p.select_sigma,
        mu1=mu1, mu2=mu2
    )

    N = p.n1 * p.n2
    ord = make_loss_order(rng, p.geometry, env1)

    EA  = zeros(Float64, length(fgrid))
    EAB = zeros(Float64, length(fgrid))
    dE  = zeros(Float64, length(fgrid))

    # baseline f=0
    keep0 = trues(N)
    EA0  = extinction_count_Aonly(A, keep0; Emin=p.Emin)
    EAB0 = extinction_count_AB_strong(A, keep0, TL, preylist; Emin=p.Emin)

    for (t, f) in enumerate(fgrid)
        keep = keep_from_order(ord, f, N)
        EA[t]  = extinction_count_Aonly(A, keep; Emin=p.Emin)
        EAB[t] = extinction_count_AB_strong(A, keep, TL, preylist; Emin=p.Emin)
        dE[t]  = (EAB[t] - EAB0) - (EA[t] - EA0)
    end

    return EA, EAB, dE
end

function run_sweep(; seed::Int=1234, reps::Int=20)
    rng = MersenneTwister(seed)

    fgrid = collect(0.0:0.05:0.90)
    t50 = findfirst(==(0.50), fgrid)

    # Sweep
    kprey_list   = [1, 3, 6]
    selects_list = [0.25, 0.6, 1.5]          # smaller => prey more similar (higher synchrony)
    geoms        = [:random, :cluster, :front]

    # Tuned defaults to get earlier action
    base = RunParams(
        80, 80,      # grid
        80,          # species
        0.25,        # basal fraction
        4,           # Lmax
        0.8,         # niche_sigma  (smaller ranges than v2)
        0.45,        # niche_cut    (smaller ranges than v2)
        3,           # k_prey (overwritten)
        0.6,         # select_sigma (overwritten)
        1.0,         # match_sigma  (slightly looser feasibility than v2)
        :random,     # geometry (overwritten)
        50           # Emin (higher than v2)
    )

    # Store mean curves
    EAmean  = Dict{Tuple{Symbol,Int,Float64}, Vector{Float64}}()
    EABmean = Dict{Tuple{Symbol,Int,Float64}, Vector{Float64}}()
    dEmean  = Dict{Tuple{Symbol,Int,Float64}, Vector{Float64}}()

    # For global y-lims (dE*)
    dE_all_min = +Inf
    dE_all_max = -Inf

    println("\n===== SUMMARY (means over reps) =====")
    println("geom,k_prey,select_sigma,dE*(0.50),AUC(dE*)")

    for g in geoms, kp in kprey_list, ss in selects_list
        p = RunParams(base.n1, base.n2, base.S, base.basal_frac, base.Lmax,
                      base.niche_sigma, base.niche_cut,
                      kp, ss, base.match_sigma, g, base.Emin)

        EA_sum  = zeros(Float64, length(fgrid))
        EAB_sum = zeros(Float64, length(fgrid))
        dE_sum  = zeros(Float64, length(fgrid))

        for _ in 1:reps
            rrng = MersenneTwister(rand(rng, UInt))
            EA, EAB, dE = run_one(rrng, p, fgrid)
            EA_sum  .+= EA
            EAB_sum .+= EAB
            dE_sum  .+= dE
        end

        EA_m  = EA_sum  ./ reps
        EAB_m = EAB_sum ./ reps
        dE_m  = dE_sum  ./ reps

        key = (g, kp, ss)
        EAmean[key]  = EA_m
        EABmean[key] = EAB_m
        dEmean[key]  = dE_m

        # stats
        dE50 = dE_m[t50]
        auc  = auc_trapz(fgrid, dE_m)

        @printf("%s,%d,%.2f,%.4f,%.4f\n", String(g), kp, ss, dE50, auc)

        dE_all_min = min(dE_all_min, minimum(dE_m))
        dE_all_max = max(dE_all_max, maximum(dE_m))
    end

    # ---------- plotting ----------
    fig = Figure(size=(1500, 520))

    # color assignment per (kp,ss) consistent across panels
    colors = Makie.wong_colors()
    combs = [(kp, ss) for kp in kprey_list for ss in selects_list]
    color_of = Dict{Tuple{Int,Float64}, Any}()
    for (i, (kp, ss)) in enumerate(combs)
        color_of[(kp, ss)] = colors[1 + (i-1) % length(colors)]
    end

    axes = Dict{Symbol, Axis}()
    for (j, g) in enumerate(geoms)
        ax = Axis(fig[1, j],
            xlabel = "habitat loss f",
            ylabel = (j == 1 ? "EA/EAB (faint) and dE* (bold)" : ""),
            title  = "geometry = $(String(g))"
        )
        axes[g] = ax
    end

    # consistent y-lims for dE* scale visibility (slight padding)
    pad = 0.05 * (dE_all_max - dE_all_min + 1e-9)
    ylow = dE_all_min - pad
    yhigh = dE_all_max + pad

    for g in geoms
        ax = axes[g]

        for kp in kprey_list, ss in selects_list
            key = (g, kp, ss)
            c = color_of[(kp, ss)]

            EA_m  = EAmean[key]
            EAB_m = EABmean[key]
            dE_m  = dEmean[key]

            # faint EA and EAB overlays (same color)
            lines!(ax, fgrid, EA_m;  color=(c, 0.18), linestyle=:dot,  linewidth=1)
            lines!(ax, fgrid, EAB_m; color=(c, 0.18), linestyle=:dash, linewidth=1)

            # bold dE* (same color) + legend label
            lines!(ax, fgrid, dE_m; color=c, linewidth=3,
                   label=@sprintf("k=%d, sel=%.2f", kp, ss))
        end

        # We want to see dE* scale; set y-lims to dE* envelope (EA/EAB will exceed, but faint)
        ylims!(ax, (ylow, yhigh))
    end

    axislegend(axes[:front]; position=:rb)  # one legend, on the rightmost panel
    display(fig)

    return nothing
end

# ----------------------------
# Main
# ----------------------------
run_sweep(seed=1234, reps=20)
