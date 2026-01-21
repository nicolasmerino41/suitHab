# suitHab_MVP_static.jl
#
# Static MVP:
#   1) Synthetic landscape (2 smooth env fields) + abiotic niches → A_i(x)
#   2) Habitat loss with geometry g and fraction f → H_{f,g}(x)
#   3) Synthetic metaweb (trophic levels + diets) with k (redundancy) + select_sigma (prey similarity)
#   4) Fixed local prey-support rule (>= 1 prey present in the cell) → AB_i(x)
#   5) Effective habitat via phi_i(f,g)=AB_i/A_i, extinctions via Emin
#   6) SAR baseline fit from pre-loss richness~area, predict extinctions from remaining area only
#   7) Show how geometry affects phi/fragmentation and how SAR error changes (static)

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

# Connected components on a Bool grid (4-neighborhood)
function largest_component_fraction(mask::BitMatrix)
    n1, n2 = size(mask)
    total = count(mask)
    total == 0 && return 0.0

    seen = falses(n1, n2)
    maxsz = 0
    q_i = Vector{Int}(undef, n1*n2)
    q_j = Vector{Int}(undef, n1*n2)

    @inbounds for r in 1:n1, c in 1:n2
        if mask[r,c] && !seen[r,c]
            head = 1
            tail = 1
            q_i[1] = r
            q_j[1] = c
            seen[r,c] = true
            sz = 0

            while head <= tail
                i = q_i[head]
                j = q_j[head]
                head += 1
                sz += 1

                # neighbors
                if i > 1 && mask[i-1,j] && !seen[i-1,j]
                    tail += 1; q_i[tail] = i-1; q_j[tail] = j; seen[i-1,j] = true
                end
                if i < n1 && mask[i+1,j] && !seen[i+1,j]
                    tail += 1; q_i[tail] = i+1; q_j[tail] = j; seen[i+1,j] = true
                end
                if j > 1 && mask[i,j-1] && !seen[i,j-1]
                    tail += 1; q_i[tail] = i; q_j[tail] = j-1; seen[i,j-1] = true
                end
                if j < n2 && mask[i,j+1] && !seen[i,j+1]
                    tail += 1; q_i[tail] = i; q_j[tail] = j+1; seen[i,j+1] = true
                end
            end

            maxsz = max(maxsz, sz)
        end
    end
    return maxsz / total
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
# Metaweb with redundancy + prey similarity control
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
# Habitat loss geometries (nested removal order)
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
# Core: compute A/AB counts, richness, phi, fragmentation
# Fixed prey rule: consumer present in cell if ANY prey present in cell.
# ----------------------------
function compute_counts!(Acount::Vector{Int}, ABcount::Vector{Int},
                         A::BitMatrix, keep::BitVector,
                         TL::Vector{Int}, preylist::Vector{Vector{Int}},
                         n1::Int, n2::Int)

    S, N = size(A)
    fill!(Acount, 0)
    fill!(ABcount, 0)

    # A counts
    @inbounds for i in 1:S
        c = 0
        for k in 1:N
            c += (keep[k] & A[i,k]) ? 1 : 0
        end
        Acount[i] = c
    end

    # Build AB presence map P(i,k) in trophic order
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
                        sup && break
                    end
                    P[i,k] = keep[k] & A[i,k] & sup
                end
            end
        end
    end

    # AB counts
    @inbounds for i in 1:S
        c = 0
        for k in 1:N
            c += P[i,k] ? 1 : 0
        end
        ABcount[i] = c
    end

    # Consumer-supported overlap map for fragmentation metric:
    # any consumer present in cell (TL>1)
    cons_mask = falses(n1, n2)
    @inbounds for k in 1:N
        anyc = false
        for i in 1:S
            if TL[i] > 1 && P[i,k]
                anyc = true
                break
            end
        end
        if anyc
            r = (k - 1) ÷ n2 + 1
            c = (k - 1) % n2 + 1
            cons_mask[r,c] = true
        end
    end

    lcc = largest_component_fraction(BitMatrix(cons_mask))
    return lcc
end

# ----------------------------
# SAR fit (pre-loss): richness vs contiguous square area fractions
# ----------------------------
function sar_fit_from_A(rng::AbstractRNG, A::BitMatrix; n1::Int, n2::Int, Emin::Int,
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

            # richness = count species with >=Emin suitable cells within square
            cnt = 0
            @inbounds for i in 1:S
                cells = 0
                for k in 1:N
                    cells += (keep_sq[k] & A[i,k]) ? 1 : 0
                end
                cnt += (cells >= Emin)
            end
            rich[s] = cnt
        end

        push!(a_eff, frac_eff)
        push!(Sbar, mean(rich))
    end

    # Fit log(S) = log(c) + z log(a)
    xs = Float64[]
    ys = Float64[]
    @inbounds for i in eachindex(a_eff)
        if Sbar[i] > 0 && a_eff[i] > 0
            push!(xs, log(a_eff[i]))
            push!(ys, log(Sbar[i]))
        end
    end
    length(xs) < 2 && return (NaN, NaN, NaN, a_eff, Sbar)

    x̄ = mean(xs); ȳ = mean(ys)
    z = sum((xs .- x̄) .* (ys .- ȳ)) / sum((xs .- x̄).^2)
    b = ȳ - z*x̄
    c = exp(b)

    # R2 on log scale
    yhat = b .+ z .* xs
    ssr = sum((ys .- yhat).^2)
    sst = sum((ys .- ȳ).^2)
    r2 = (sst > 0) ? (1 - ssr/sst) : NaN

    return (c, z, r2, a_eff, Sbar)
end

# ----------------------------
# MVP runner
# ----------------------------
struct MVPParams
    n1::Int; n2::Int; S::Int
    basal_frac::Float64; Lmax::Int
    niche_sigma::Float64; niche_cut::Float64
    k_prey::Int; select_sigma::Float64; match_sigma::Float64
    Emin::Int
    reps::Int
    fgrid::Vector{Float64}
    geoms::Vector{Symbol}
    sar_area_fracs::Vector{Float64}
    sar_samples::Int
end

function run_mvp(; seed::Int=1234,
    n1::Int=80, n2::Int=80, S::Int=80,
    basal_frac::Float64=0.25, Lmax::Int=4,
    niche_sigma::Float64=0.8, niche_cut::Float64=0.45,
    k_prey::Int=3, select_sigma::Float64=0.6, match_sigma::Float64=1.0,
    Emin::Int=50,
    reps::Int=25,
    fgrid = collect(0.0:0.05:0.90),
    geoms = [:random, :cluster, :front],
    sar_area_fracs = [0.25, 0.35, 0.50, 0.70, 0.85, 1.0],
    sar_samples::Int=25
)
    p = MVPParams(n1,n2,S,basal_frac,Lmax,niche_sigma,niche_cut,
                  k_prey,select_sigma,match_sigma,
                  Emin,reps,Vector{Float64}(fgrid),
                  Vector{Symbol}(geoms),
                  Vector{Float64}(sar_area_fracs),sar_samples)

    rng = MersenneTwister(seed)
    N = p.n1*p.n2
    T = length(p.fgrid)

    EA_inc  = Dict{Symbol, Vector{Float64}}()
    EAB_inc = Dict{Symbol, Vector{Float64}}()
    phi_cons = Dict{Symbol, Vector{Float64}}()
    lcc_cons = Dict{Symbol, Vector{Float64}}()

    for g in p.geoms
        EA_inc[g]  = zeros(T)
        EAB_inc[g] = zeros(T)
        phi_cons[g] = zeros(T)
        lcc_cons[g] = zeros(T)
    end

    c_list = Float64[]; z_list = Float64[]; r2_list = Float64[]
    S0_list = Float64[]; SAB0_list = Float64[]; Dprey_list = Float64[]
    Acount = zeros(Int, p.S); ABcount = zeros(Int, p.S)

    for rep in 1:p.reps
        rrng = MersenneTwister(rand(rng, UInt))

        env1, env2 = make_environment(rrng; n1=p.n1, n2=p.n2)
        A, mu1, mu2 = make_abiotic_maps(rrng, env1, env2;
                                       S=p.S, niche_sigma=p.niche_sigma,
                                       niche_cut=p.niche_cut)

        TL, preylist = build_metaweb(rrng;
            S=p.S, basal_frac=p.basal_frac, Lmax=p.Lmax,
            k_prey=p.k_prey, match_sigma=p.match_sigma,
            select_sigma=p.select_sigma,
            mu1=mu1, mu2=mu2
        )

        push!(Dprey_list, realized_prey_diversity(mu1, mu2, TL, preylist))

        c, z, r2, _, _ = sar_fit_from_A(rrng, A;
            n1=p.n1, n2=p.n2, Emin=p.Emin,
            area_fracs=p.sar_area_fracs, samples=p.sar_samples)

        isfinite(c) && isfinite(z) && (push!(c_list, c); push!(z_list, z); push!(r2_list, r2))

        ord = Dict{Symbol, Vector{Int}}()
        for g in p.geoms
            ord[g] = make_loss_order(rrng, g, env1)
        end

        keep0 = trues(N)
        compute_counts!(Acount, ABcount, A, keep0, TL, preylist, p.n1, p.n2)
        S0  = count(>=(p.Emin), Acount)
        SAB0 = count(>=(p.Emin), ABcount)
        push!(S0_list, S0); push!(SAB0_list, SAB0)

        for g in p.geoms
            for (t, f) in enumerate(p.fgrid)
                keep = (f == 0.0) ? keep0 : keep_from_order(ord[g], f, N)
                lcc = compute_counts!(Acount, ABcount, A, keep, TL, preylist, p.n1, p.n2)

                SA  = count(>=(p.Emin), Acount)
                SAB = count(>=(p.Emin), ABcount)

                EA_inc[g][t]  += (S0  - SA)
                EAB_inc[g][t] += (SAB0 - SAB)

                sφ = 0.0; nφ = 0
                for i in 1:p.S
                    if TL[i] > 1 && Acount[i] > 0
                        sφ += ABcount[i] / Acount[i]
                        nφ += 1
                    end
                end
                phi_cons[g][t] += (nφ > 0 ? sφ/nφ : NaN)
                lcc_cons[g][t] += lcc
            end
        end
    end

    for g in p.geoms
        EA_inc[g]  ./= p.reps
        EAB_inc[g] ./= p.reps
        phi_cons[g] ./= p.reps
        lcc_cons[g] ./= p.reps
    end

    S0m = mean(S0_list)
    cm, zm = mean(c_list), mean(z_list)

    SAR_pred = [clamp(cm * max(1e-9, 1-f)^zm, 0, S0m) |> x -> S0m - x
                for f in p.fgrid]

    # ----------------------------
    # Plot
    # ----------------------------
    fig = Figure(size=(1650, 900))

    # 🔑 GLOBAL y-limits for row 1
    yvals = Float64[]
    for g in p.geoms
        append!(yvals, EA_inc[g])
        append!(yvals, EAB_inc[g])
    end
    append!(yvals, SAR_pred)

    ymin = minimum(yvals)
    ymax = maximum(yvals)
    pad = 0.05 * (ymax - ymin + 1e-9)
    ylims_row1 = (ymin - pad, ymax + pad)

    for (j, g) in enumerate(p.geoms)
        ax1 = Axis(fig[1, j],
            xlabel = "habitat loss f",
            ylabel = (j == 1 ? "incremental extinctions" : ""),
            title  = "geometry = $(String(g))"
        )

        ylims!(ax1, ylims_row1)

        lines!(ax1, p.fgrid, EA_inc[g]; linewidth=3, label="A-only obs")
        lines!(ax1, p.fgrid, EAB_inc[g]; linewidth=3, label="AB obs (biotic)")
        lines!(ax1, p.fgrid, SAR_pred; linewidth=3, linestyle=:dash,
               label="SAR baseline")

        axislegend(ax1; position=:lt)

        ax2L = Axis(fig[2, j],
            xlabel = "habitat loss f",
            ylabel = (j == 1 ? "mean φ (consumers)" : "")
        )
        ax2R = Axis(fig[2, j],
            yaxisposition=:right,
            ylabel = (j == 3 ? "LCC fraction" : ""),
            xgridvisible=false, ygridvisible=false
        )
        linkxaxes!(ax2L, ax2R)
        hidespines!(ax2R); hidexdecorations!(ax2R)

        ylims!(ax2L, (0,1)); ylims!(ax2R, (0,1))
        lines!(ax2L, p.fgrid, phi_cons[g]; linewidth=3)
        lines!(ax2R, p.fgrid, lcc_cons[g]; linewidth=3, linestyle=:dash)
    end

    display(fig)
end

# ----------------------------
# Run (edit params here)
# ----------------------------
run_mvp(
    seed=1234,
    reps=25,
    k_prey=3,
    select_sigma=0.6,
    Emin=50
)