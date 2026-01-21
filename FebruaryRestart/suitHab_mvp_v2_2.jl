# suitHab_mvp_v2_2.jl
#
# Synthetic MVP v2.2:
#   - Twin-axis panels per geometry: EA/EAB (left) + dE* (right)
#   - Realized prey diversity diagnostic (proxy for "anti-synchrony"):
#       Dprey = mean pairwise distance among prey optima within each consumer prey set
#     (smaller Dprey => prey more similar => higher synchrony)
#   - Second row: heatmaps of Dprey across (k_prey, select_sigma) per geometry
#   - Prints summary table: dE*(0.5), AUC(dE*), Dprey_mean
#   - No saving, no CSV, always display(fig)

using Random
using Statistics
using Printf
using CairoMakie

# ----------------------------
# Utilities
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

function auc_trapz(x::Vector{Float64}, y::Vector{Float64})
    s = 0.0
    @inbounds for i in 1:length(x)-1
        dx = x[i+1] - x[i]
        s += dx * (y[i] + y[i+1]) / 2.0
    end
    return s
end

# ----------------------------
# Synthetic environment + abiotic niches
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
# Metaweb + realized prey diversity diagnostic
# ----------------------------
"""
build_metaweb with prey-guild construction:
- pick a guild seed prey g close to consumer (match_sigma)
- pick remaining prey clustered around g (select_sigma controls prey-prey similarity / synchrony)

Returns TL, preylist.
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

        # consumer-prey feasibility weights
        wmatch = Vector{Float64}(undef, length(cands))
        for (t, j) in enumerate(cands)
            d1 = mu1[i] - mu1[j]
            d2 = mu2[i] - mu2[j]
            d2sum = d1*d1 + d2*d2
            wmatch[t] = exp(-d2sum / (2.0 * match_sigma^2))
        end

        # guild seed
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

        # weights: close to guild seed (select_sigma) AND feasible for consumer (match_sigma)
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

    return TL, preylist
end

"""
Realized prey diversity diagnostic:
For each consumer with m>=2 prey, compute mean pairwise Euclidean distance among prey optima.
Return average across consumers.
Smaller => prey more similar => higher synchrony.
"""
function realized_prey_diversity(mu1::Vector{Float64}, mu2::Vector{Float64},
                                 TL::Vector{Int}, preylist::Vector{Vector{Int}})
    S = length(TL)
    vals = Float64[]
    @inbounds for i in 1:S
        if TL[i] == 1
            continue
        end
        prey = preylist[i]
        m = length(prey)
        if m < 2
            continue
        end
        s = 0.0
        np = 0
        for a in 1:m-1
            ja = prey[a]
            for b in a+1:m
                jb = prey[b]
                d1 = mu1[ja] - mu1[jb]
                d2 = mu2[ja] - mu2[jb]
                s += sqrt(d1*d1 + d2*d2)
                np += 1
            end
        end
        if np > 0
            push!(vals, s / np)
        end
    end
    return isempty(vals) ? NaN : mean(vals)
end

# ----------------------------
# Habitat loss: nested removals by a single order
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

    return sortperm(score; rev=true)
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
# Extinction counts with Emin threshold
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
# Runner
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
One replicate returns:
  EA(f), EAB(f), dE*(f), Dprey
where:
  dE*(f) = (EAB(f)-EAB(0)) - (EA(f)-EA(0))
  Dprey  = realized prey diversity (mean pairwise prey distance within consumers)
"""
function run_one(rng::AbstractRNG, p::RunParams, fgrid::Vector{Float64})
    env1, env2 = make_environment(rng; n1=p.n1, n2=p.n2)
    A, mu1, mu2 = make_abiotic_maps(rng, env1, env2; S=p.S, niche_sigma=p.niche_sigma, niche_cut=p.niche_cut)

    TL, preylist = build_metaweb(rng;
        S=p.S, basal_frac=p.basal_frac, Lmax=p.Lmax,
        k_prey=p.k_prey, match_sigma=p.match_sigma, select_sigma=p.select_sigma,
        mu1=mu1, mu2=mu2
    )

    Dprey = realized_prey_diversity(mu1, mu2, TL, preylist)

    N = p.n1 * p.n2
    ord = make_loss_order(rng, p.geometry, env1)

    EA  = zeros(Float64, length(fgrid))
    EAB = zeros(Float64, length(fgrid))
    dE  = zeros(Float64, length(fgrid))

    keep0 = trues(N)
    EA0  = extinction_count_Aonly(A, keep0; Emin=p.Emin)
    EAB0 = extinction_count_AB_strong(A, keep0, TL, preylist; Emin=p.Emin)

    for (t, f) in enumerate(fgrid)
        keep = keep_from_order(ord, f, N)
        EA[t]  = extinction_count_Aonly(A, keep; Emin=p.Emin)
        EAB[t] = extinction_count_AB_strong(A, keep, TL, preylist; Emin=p.Emin)
        dE[t]  = (EAB[t] - EAB0) - (EA[t] - EA0)
    end

    return EA, EAB, dE, Dprey
end

function run_sweep(; seed::Int=1234, reps::Int=40)
    rng = MersenneTwister(seed)

    fgrid = collect(0.0:0.01:0.8)
    t50 = findfirst(==(0.50), fgrid)

    kprey_list   = [1, 3, 6]
    selects_list = [0.25, 0.6, 1.5]
    geoms        = [:random, :cluster, :front]

    # Defaults (kept from v2.1)
    base = RunParams(
        80, 80,
        80,
        0.25,
        4,
        0.8,
        0.45,
        3,      # k overwritten
        0.6,    # sel overwritten
        1.0,    # match_sigma
        :random,
        50      # Emin
    )

    EAmean  = Dict{Tuple{Symbol,Int,Float64}, Vector{Float64}}()
    EABmean = Dict{Tuple{Symbol,Int,Float64}, Vector{Float64}}()
    dEmean  = Dict{Tuple{Symbol,Int,Float64}, Vector{Float64}}()
    Dmean   = Dict{Tuple{Symbol,Int,Float64}, Float64}()

    dE_all_min = +Inf
    dE_all_max = -Inf

    println("\n===== SUMMARY (means over reps) =====")
    println("geom,k_prey,select_sigma,dE*(0.50),AUC(dE*),Dprey_mean  # Dprey small => high synchrony")

    for g in geoms, kp in kprey_list, ss in selects_list
        p = RunParams(base.n1, base.n2, base.S, base.basal_frac, base.Lmax,
                      base.niche_sigma, base.niche_cut,
                      kp, ss, base.match_sigma, g, base.Emin)

        EA_sum  = zeros(Float64, length(fgrid))
        EAB_sum = zeros(Float64, length(fgrid))
        dE_sum  = zeros(Float64, length(fgrid))
        D_sum   = 0.0

        for _ in 1:reps
            rrng = MersenneTwister(rand(rng, UInt))
            EA, EAB, dE, Dprey = run_one(rrng, p, fgrid)
            EA_sum  .+= EA
            EAB_sum .+= EAB
            dE_sum  .+= dE
            D_sum   += Dprey
        end

        EA_m  = EA_sum  ./ reps
        EAB_m = EAB_sum ./ reps
        dE_m  = dE_sum  ./ reps
        D_m   = D_sum / reps

        key = (g, kp, ss)
        EAmean[key]  = EA_m
        EABmean[key] = EAB_m
        dEmean[key]  = dE_m
        Dmean[key]   = D_m

        dE50 = dE_m[t50]
        auc  = auc_trapz(fgrid, dE_m)

        @printf("%s,%d,%.2f,%.4f,%.4f,%.4f\n", String(g), kp, ss, dE50, auc, D_m)

        dE_all_min = min(dE_all_min, minimum(dE_m))
        dE_all_max = max(dE_all_max, maximum(dE_m))
    end

    # ----------------------------
    # Plot: row 1 = twin-axis panels; row 2 = heatmaps of Dprey
    # ----------------------------

    fig = Figure(size=(1650, 900))

    colors = Makie.wong_colors()
    combs = [(kp, ss) for kp in kprey_list for ss in selects_list]
    color_of = Dict{Tuple{Int,Float64}, Any}()
    for (i, (kp, ss)) in enumerate(combs)
        color_of[(kp, ss)] = colors[1 + (i-1) % length(colors)]
    end

    # dE* y-lims shared across geometry (right axis)
    pad = 0.05 * (dE_all_max - dE_all_min + 1e-9)
    yR_lo = dE_all_min - pad
    yR_hi = dE_all_max + pad

    axesL = Dict{Symbol, Axis}()
    axesR = Dict{Symbol, Axis}()
    heat_axes = Dict{Symbol, Axis}()

    for (j, g) in enumerate(geoms)
        # ----- Row 1: twin-axis overlay -----
        axL = Axis(fig[1, j],
            xlabel = "habitat loss f",
            ylabel = (j == 1 ? "EA / EAB" : ""),
            title  = "geometry = $(String(g))"
        )
        axR = Axis(fig[1, j],
            ylabel = (j == 3 ? "dE* (right axis)" : ""),
            yaxisposition = :right,
            xgridvisible = false,
            ygridvisible = false
        )
        linkxaxes!(axL, axR)
        hidespines!(axR)
        hidexdecorations!(axR)

        axesL[g] = axL
        axesR[g] = axR

        # ----- Row 2: heatmap of realized prey diversity -----
        axH = Axis(fig[2, j],
            xlabel = "select_sigma",
            ylabel = (j == 1 ? "k_prey" : ""),
            title  = "realized prey diversity (Dprey)"
        )
        heat_axes[g] = axH
    end

    # Legend only once (rightmost twin-axis panel)
    legend_ax = axesR[:front]

    # Plot curves
    for g in geoms
        axL = axesL[g]
        axR = axesR[g]

        # Right axis y-lims fixed for comparability
        ylims!(axR, (yR_lo, yR_hi))

        for kp in kprey_list, ss in selects_list
            key = (g, kp, ss)
            c = color_of[(kp, ss)]

            EA_m  = EAmean[key]
            EAB_m = EABmean[key]
            dE_m  = dEmean[key]

            # EA/EAB on left axis (faint)
            lines!(axL, fgrid, EA_m;  color=(c, 0.18), linestyle=:dot,  linewidth=1)
            lines!(axL, fgrid, EAB_m; color=(c, 0.18), linestyle=:dash, linewidth=1)

            # dE* on right axis (bold)
            lines!(axR, fgrid, dE_m; color=c, linewidth=3,
                   label=@sprintf("k=%d, sel=%.2f", kp, ss))
        end
    end

    axislegend(legend_ax; position=:lt)

    # Heatmaps: Dprey matrices (k rows, sel cols)
    # We will set consistent color limits across geoms for comparability.
    Dvals = Float64[]
    for g in geoms, kp in kprey_list, ss in selects_list
        v = Dmean[(g, kp, ss)]
        isfinite(v) && push!(Dvals, v)
    end

    @assert !isempty(Dvals) "All Dprey values are NaN — cannot plot heatmap."

    Dmin = minimum(Dvals)
    Dmax = maximum(Dvals)

    # Avoid degenerate color range
    if Dmin == Dmax
        Dmin -= 1e-6
        Dmax += 1e-6
    end

    for g in geoms
        axH = heat_axes[g]

        M = zeros(Float64, length(kprey_list), length(selects_list))
        mask = falses(size(M))

        for (i, kp) in enumerate(kprey_list)
            for (j, ss) in enumerate(selects_list)
                v = Dmean[(g, kp, ss)]
                if isfinite(v)
                    M[i, j] = v
                else
                    M[i, j] = Dmin        # dummy value
                    mask[i, j] = true    # mark as invalid
                end
            end
        end

        heatmap!(
            axH,
            selects_list,
            kprey_list,
            M;
            colorrange = (Dmin, Dmax)
        )

        # visually mark invalid cells
        for I in findall(mask)
            i, j = Tuple(I)
            text!(axH,
                selects_list[j], kprey_list[i],
                text="∅", align=(:center, :center), color=:black
            )
        end

        axH.xticks = (selects_list, string.(selects_list))
        axH.yticks = (kprey_list, string.(kprey_list))
    end

    display(fig)
end

# ----------------------------
# Main
# ----------------------------
run_sweep(seed=1234, reps=40)