###############################################################################
# patchonly_article_ready_independent_v2_metawebER.jl
#
# Single, independent Julia script (no external data):
#   Run it -> generates synthetic environment + niches, builds an ER-like META-WEB
#   (directed, all-species; basal fraction; optional niche assortativity via corr),
#   simulates A-only vs AB across habitat-loss geometries, enforces monotonic sanity,
#   computes raw + relative ΔAUC, and produces 3 figures:
#
#   Fig1 (MAIN): 6 line panels = (low divergence + high divergence) × (3 HL geometries)
#                using Emin_patch metric (connectivity-sensitive).
#                Also produces SI version for viability metric (area-only).
#
#   Fig2 (MAIN): Heatmap sweep (corr × connectance) for multiple niche scenarios
#                showing (i) raw ΔAUC and (ii) relative ΔAUC (ΔAUC/AUC_Aonly).
#                Uses chosen geometry (default: cluster) for main.
#                NOTE: Heatmaps are computed on a coarse grid then bilinear-upsampled
#                for smooth continuous-looking transitions (article-friendly).
#
#   Fig3 (MECH): "support amount vs support connectivity" mechanism panel (3 geometries)
#
# Dependencies: only Random, Statistics, Printf, CairoMakie (standard + Makie).
###############################################################################

using Random, Statistics, Printf
using CairoMakie

# ----------------------------- GLOBAL SETTINGS --------------------------------
const SEED = 12345
Random.seed!(SEED)

# Landscape
const NX = 80
const NY = 80
const NCELLS = NX * NY

# Community
const N_SPECIES = 250
const BASAL_FRAC = 0.33
const N_BASAL = Int(round(BASAL_FRAC * N_SPECIES))
const N_CONS  = N_SPECIES - N_BASAL

# Habitat-loss grid
const F_MAX = 1.0
const NF    = 46
const FGRID = range(0.0, F_MAX, length=NF)

# Parameter sweep (coarse -> upsample for smooth heatmaps)
const CORR_GRID_COARSE = collect(0.0:0.2:1.0)        # 6 values
const PCONN_GRID_COARSE = collect(0.02:0.02:0.20)    # 10 values (target connectance p)

# Heatmap plot resolution (fine grid for display only)
const CORR_GRID_FINE  = collect(0.0:0.05:1.0)
const PCONN_GRID_FINE = collect(0.02:0.01:0.20)

# Geometries
const GEOMS = ["random", "cluster", "front"]
const HEATMAP_GEOM = "cluster"

# Presence thresholds
const A_MIN = 15
const LCC_MIN = 12

# Consumer trophic requirement φ
const PHI_MEAN = 0.80
const PHI_SD   = 0.05
const PHI_MIN  = 0.50
const PHI_MAX  = 0.95

# Niche threshold
const THETA = 0.55

# Link-assortativity length scale in niche space (for corr>0)
const L_NICHE = 0.15

# Monotonic sanity check: richness must not increase with habitat loss
const ENFORCE_MONOTONE = true
const APPLY_MONO_FIX   = true

# Output
const OUT_DIR = joinpath(pwd(), "paper_final_independent_v2_metawebER")
const SAVE_FIGS = true
mkpath(OUT_DIR)

# Style
set_theme!(Theme(
    fontsize = 16,
    Axis = (xlabelsize=15, ylabelsize=15, titlesize=15,
            xticklabelsize=13, yticklabelsize=13),
    Legend = (labelsize=13,)
))

# ----------------------------- HELPERS ----------------------------------------
clamp01(x) = min(max(x, 0.0), 1.0)

function smooth_field!(A::Matrix{Float64}; iters::Int=6)
    ny, nx = size(A)
    B = similar(A)
    for _ in 1:iters
        @inbounds for y in 1:ny, x in 1:nx
            s = 0.0
            c = 0
            for (yy, xx) in ((y,x),
                             (max(y-1,1), x),
                             (min(y+1,ny), x),
                             (y, max(x-1,1)),
                             (y, min(x+1,nx)))
                s += A[yy, xx]; c += 1
            end
            B[y, x] = s / c
        end
        A, B = B, A
    end
    return A
end

@inline linidx(y::Int, x::Int) = (y-1) * NX + x

function build_neighbors()
    neigh = Vector{NTuple{4,Int}}(undef, NCELLS)
    @inbounds for y in 1:NY, x in 1:NX
        i = linidx(y,x)
        up    = (y > 1)  ? linidx(y-1,x) : 0
        down  = (y < NY) ? linidx(y+1,x) : 0
        left  = (x > 1)  ? linidx(y,x-1) : 0
        right = (x < NX) ? linidx(y,x+1) : 0
        neigh[i] = (up, down, left, right)
    end
    return neigh
end
const NEIGH = build_neighbors()

function loss_priority(geom::String; rng=Random.default_rng())
    if geom == "random"
        return rand(rng, Float64, NY, NX)
    elseif geom == "front"
        P = zeros(Float64, NY, NX)
        @inbounds for y in 1:NY, x in 1:NX
            P[y,x] = (x-1) / (NX-1)
        end
        return P
    elseif geom == "cluster"
        P = randn(rng, Float64, NY, NX)
        P = smooth_field!(P; iters=7)
        mn, mx = minimum(P), maximum(P)
        P .= (P .- mn) ./ (mx - mn + 1e-12)
        return P
    else
        error("Unknown geometry: $geom")
    end
end

# Precompute remain masks for all geometries and f to avoid repeated sorting/thresholding.
function precompute_remain_masks(; rng=Random.default_rng())
    masks = Dict{Tuple{String,Int}, BitMatrix}()
    for geom in GEOMS
        P = loss_priority(geom; rng=rng)
        v = vec(copy(P)); sort!(v)
        for (ti, f) in enumerate(FGRID)
            k = clamp(Int(floor(f * length(v))), 0, length(v)-1)
            thr = v[max(k,1)]
            R = trues(NY, NX)
            @inbounds for y in 1:NY, x in 1:NX
                if P[y,x] <= thr && f > 0
                    R[y,x] = false
                end
            end
            masks[(geom, ti)] = R
        end
    end
    return masks
end

function make_environment(; rng=Random.default_rng())
    E = zeros(Float64, NY, NX)
    noise = randn(rng, Float64, NY, NX)
    noise = smooth_field!(noise; iters=5)
    mn, mx = minimum(noise), maximum(noise)
    noise .= (noise .- mn) ./ (mx - mn + 1e-12)
    @inbounds for y in 1:NY, x in 1:NX
        grad = 0.7 * (x-1)/(NX-1) + 0.3 * (y-1)/(NY-1)
        E[y,x] = clamp01(0.65*grad + 0.35*noise[y,x])
    end
    return E
end

@inline function gauss_suit(e::Float64, μ::Float64, σ::Float64)
    z = (e - μ) / (σ + 1e-12)
    return exp(-0.5 * z*z)
end

function potential_mask(E::Matrix{Float64}, μ::Float64, σ::Float64; theta=THETA)
    M = falses(NY, NX)
    @inbounds for y in 1:NY, x in 1:NX
        M[y,x] = gauss_suit(E[y,x], μ, σ) > theta
    end
    return M
end

# Largest connected component for intersection (pot ∧ remain).
# Returns: (area_total, lcc_size, lcc_indices)
function lcc_intersection(pot::BitMatrix, remain::BitMatrix, visited::Vector{UInt8}, stack::Vector{Int})
    fill!(visited, 0x00)
    best_size = 0
    best_indices = Int[]
    area_total = 0

    @inbounds for y in 1:NY, x in 1:NX
        if pot[y,x] & remain[y,x]
            area_total += 1
        end
    end
    area_total == 0 && return (0, 0, best_indices)

    @inbounds for y in 1:NY, x in 1:NX
        if pot[y,x] & remain[y,x]
            i = linidx(y,x)
            if visited[i] == 0x01
                continue
            end
            comp = Int[]
            empty!(stack)
            push!(stack, i)
            visited[i] = 0x01
            while !isempty(stack)
                u = pop!(stack)
                push!(comp, u)
                (up,down,left,right) = NEIGH[u]
                if up != 0
                    yy = (up - 1) ÷ NX + 1; xx = (up - 1) % NX + 1
                    if visited[up] == 0x00 && (pot[yy,xx] & remain[yy,xx])
                        visited[up] = 0x01; push!(stack, up)
                    end
                end
                if down != 0
                    yy = (down - 1) ÷ NX + 1; xx = (down - 1) % NX + 1
                    if visited[down] == 0x00 && (pot[yy,xx] & remain[yy,xx])
                        visited[down] = 0x01; push!(stack, down)
                    end
                end
                if left != 0
                    yy = (left - 1) ÷ NX + 1; xx = (left - 1) % NX + 1
                    if visited[left] == 0x00 && (pot[yy,xx] & remain[yy,xx])
                        visited[left] = 0x01; push!(stack, left)
                    end
                end
                if right != 0
                    yy = (right - 1) ÷ NX + 1; xx = (right - 1) % NX + 1
                    if visited[right] == 0x00 && (pot[yy,xx] & remain[yy,xx])
                        visited[right] = 0x01; push!(stack, right)
                    end
                end
            end
            if length(comp) > best_size
                best_size = length(comp)
                best_indices = comp
            end
        end
    end
    return (area_total, best_size, best_indices)
end

function auc_trapz(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    n < 2 && return 0.0
    s = 0.0
    @inbounds for i in 1:(n-1)
        dx = float(x[i+1] - x[i])
        s += 0.5 * dx * float(y[i] + y[i+1])
    end
    return s
end

function monotone_check_and_fix!(y::Vector{Float64})
    inc = 0.0
    @inbounds for i in 1:(length(y)-1)
        inc = max(inc, y[i+1] - y[i])
    end
    if inc > 1e-12
        if ENFORCE_MONOTONE
            @warn @sprintf("Monotonic violation detected: max increase = %.4g", inc)
        end
        if APPLY_MONO_FIX
            for i in 2:length(y)
                y[i] = min(y[i-1], y[i])
            end
        end
    end
    return inc
end

# ----------------------------- NICHE SCENARIOS --------------------------------
struct Scenario
    id::Int
    name::String
    σ_mode::Symbol
end

const SCENARIOS = Scenario[
    Scenario(1, "Very variable occupancy + very variable sigma", :very_variable),
    Scenario(2, "Same occupancy range, less extreme sigma variability", :less_variable),
    Scenario(3, "Skewed to narrower niches (lower mean occupancy), still variable", :narrow),
    Scenario(4, "Skewed to broader niches (higher mean occupancy), still variable", :broad),
]

function sample_sigma(mode::Symbol, n::Int; rng=Random.default_rng())
    σ = zeros(Float64, n)
    if mode == :very_variable
        for i in 1:n; σ[i] = clamp(0.07 * exp(0.9*randn(rng)), 0.02, 0.30) end
    elseif mode == :less_variable
        for i in 1:n; σ[i] = clamp(0.09 * exp(0.35*randn(rng)), 0.03, 0.28) end
    elseif mode == :narrow
        for i in 1:n; σ[i] = clamp(0.06 * exp(0.45*randn(rng)), 0.02, 0.18) end
    elseif mode == :broad
        for i in 1:n; σ[i] = clamp(0.14 * exp(0.35*randn(rng)), 0.05, 0.35) end
    else
        error("Unknown sigma mode: $mode")
    end
    return σ
end

function sample_phi(n::Int; rng=Random.default_rng())
    φ = PHI_MEAN .+ PHI_SD .* randn(rng, n)
    return clamp.(φ, PHI_MIN, PHI_MAX)
end

# ----------------------------- META-WEB BUILDER --------------------------------
"""
Build a directed metaweb with an ER baseline and optional niche assortativity.

- basal[i] == true  -> species i is basal (no outgoing links)
- For non-basal i, for each j != i:
    w_ij = exp(-|μ_i - μ_j| / L_NICHE), rescaled so mean(w)=1
    p_ij = p * ((1-corr) + corr*w_ij)
  then Bernoulli draw link i -> j with probability clamp(p_ij,0,1).

Guarantees each non-basal has at least 1 prey link (forces one if empty).
Returns:
- prey_lists[i] = Vector{Int} of prey of i
- realised_connectance = edges / (N_nonbasal*(N_SPECIES-1))
"""
function build_metaweb_ER(μ::Vector{Float64}, basal::BitVector, p::Float64, corr::Float64; rng=Random.default_rng())
    prey_lists = [Int[] for _ in 1:length(μ)]
    n = length(μ)

    # Precompute mean(w) over allowed pairs to keep expected p ~ constant.
    # Allowed pairs: i non-basal, j != i.
    wsum = 0.0
    wcnt = 0
    @inbounds for i in 1:n
        basal[i] && continue
        μi = μ[i]
        for j in 1:n
            j == i && continue
            w = exp(-abs(μi - μ[j]) / (L_NICHE + 1e-12))
            wsum += w
            wcnt += 1
        end
    end
    wmean = (wcnt > 0) ? (wsum / wcnt) : 1.0

    edges = 0
    denom = N_CONS * (n - 1)

    @inbounds for i in 1:n
        if basal[i]
            empty!(prey_lists[i])
            continue
        end
        μi = μ[i]
        for j in 1:n
            j == i && continue
            w = exp(-abs(μi - μ[j]) / (L_NICHE + 1e-12))
            w /= (wmean + 1e-12)  # rescale so mean ~ 1
            pij = p * ((1.0 - corr) + corr * w)
            pij = clamp01(pij)
            if rand(rng) < pij
                push!(prey_lists[i], j)
                edges += 1
            end
        end
        if isempty(prey_lists[i])
            # force at least one prey
            j = rand(rng, 1:n)
            j == i && (j = (i % n) + 1)
            push!(prey_lists[i], j)
            edges += 1
        end
    end

    creal = edges / (denom + 1e-12)
    return prey_lists, creal
end

# ----------------------------- CORE SIMULATION --------------------------------
"""
Simulate one (scenario, corr, p_connectance, geometry).

We generate niches (μ, σ), potential habitat masks, an ER/assortative metaweb, then
simulate A-only vs AB across habitat loss:
- A-only: species present if it has enough suitable area (viability) or enough LCC (Emin_patch)
- AB: non-basal species must additionally have enough supported prey fraction >= φ

Mechanism outputs (for Fig3) computed using Emin_patch logic:
- mean φ
- supported_LCC = mean supported prey fraction within consumer LCC
- frag_fail = fraction of consumers viable in A-only (Emin) but failing AB due to prey support
"""
function simulate_one(E::Matrix{Float64}, remain_masks, sc::Scenario, corr::Float64, pconn::Float64, geom::String;
                      rng=Random.default_rng())

    # Basal assignment: first N_BASAL are basal (deterministic ordering for reproducibility)
    basal = falses(N_SPECIES)
    for i in 1:N_BASAL
        basal[i] = true
    end

    # Niche parameters
    μ = rand(rng, N_SPECIES)
    σ = sample_sigma(sc.σ_mode, N_SPECIES; rng=rng)

    # Potential masks
    pot = Vector{BitMatrix}(undef, N_SPECIES)
    for i in 1:N_SPECIES
        pot[i] = potential_mask(E, μ[i], σ[i])
    end

    # Metaweb + φ for consumers
    prey_lists, creal = build_metaweb_ER(μ, basal, pconn, corr; rng=rng)
    φ_req = sample_phi(N_SPECIES; rng=rng)  # only used for non-basal; basal ignored
    meanφ = mean(φ_req[.!basal])

    visited = Vector{UInt8}(undef, NCELLS)
    stack = Int[]

    # outputs
    rA_E = zeros(Float64, NF)
    rB_E = zeros(Float64, NF)
    rA_V = zeros(Float64, NF)
    rB_V = zeros(Float64, NF)

    mean_phi_req = fill(meanφ, NF)
    supported_lcc = zeros(Float64, NF)
    frag_fail = zeros(Float64, NF)

    # temp presence vectors
    presentV = falses(N_SPECIES)
    presentE = falses(N_SPECIES)

    for ti in 1:NF
        remain = remain_masks[(geom, ti)]

        # Compute A-only presence (V and E) + keep each consumer LCC indices for Emin
        # (we recompute per f; simplest and consistent)
        lcc_idx_cache = Vector{Vector{Int}}(undef, N_SPECIES)
        area_cache = zeros(Int, N_SPECIES)
        lcc_cache  = zeros(Int, N_SPECIES)

        for i in 1:N_SPECIES
            area, lcc, idx = lcc_intersection(pot[i], remain, visited, stack)
            area_cache[i] = area
            lcc_cache[i]  = lcc
            lcc_idx_cache[i] = idx
            presentV[i] = (area >= A_MIN)
            presentE[i] = (area >= A_MIN) && (lcc >= LCC_MIN)
        end

        # AB presence initialised as A-only for basal (they don't require prey)
        presentV_AB = copy(presentV)
        presentE_AB = copy(presentE)

        # Mechanism tallies (Emin)
        denom_cons = 0
        supp_sum = 0.0
        fail_frag = 0

        # Apply trophic constraint to non-basal species
        for i in 1:N_SPECIES
            basal[i] && continue

            # Emin-based AB (within consumer LCC)
            eA = presentE[i]
            if eA && !isempty(lcc_idx_cache[i])
                denom_cons += 1
                prey = prey_lists[i]
                deg = length(prey)
                supported = 0
                for j in prey
                    presentE[j] || continue
                    has = false
                    @inbounds for u in lcc_idx_cache[i]
                        yy = (u - 1) ÷ NX + 1
                        xx = (u - 1) % NX + 1
                        if pot[j][yy,xx] & remain[yy,xx]
                            has = true
                            break
                        end
                    end
                    supported += has ? 1 : 0
                end
                supp = (deg > 0) ? (supported / deg) : 0.0
                supp_sum += supp

                eAB = (supp >= φ_req[i])
                presentE_AB[i] = eAB
                if eA && !eAB
                    fail_frag += 1
                end
            else
                presentE_AB[i] = false
            end

            # Viability-based AB (area-only): check prey overlap anywhere in consumer remaining suitable area.
            # To avoid scanning full grid for each prey, we reuse consumer LCC indices as a proxy for "accessible patch".
            # (This keeps viability as area-only via A_MIN, but uses patch co-location proxy for prey encounter.)
            vA = presentV[i]
            if vA && !isempty(lcc_idx_cache[i])
                prey = prey_lists[i]
                deg = length(prey)
                supported = 0
                for j in prey
                    presentV[j] || continue
                    has = false
                    @inbounds for u in lcc_idx_cache[i]
                        yy = (u - 1) ÷ NX + 1
                        xx = (u - 1) % NX + 1
                        if pot[j][yy,xx] & remain[yy,xx]
                            has = true
                            break
                        end
                    end
                    supported += has ? 1 : 0
                end
                supp = (deg > 0) ? (supported / deg) : 0.0
                presentV_AB[i] = (supp >= φ_req[i])
            else
                presentV_AB[i] = false
            end
        end

        # richness
        rA_E[ti] = sum(presentE)
        rB_E[ti] = sum(presentE_AB)
        rA_V[ti] = sum(presentV)
        rB_V[ti] = sum(presentV_AB)

        # mechanism summaries
        if denom_cons > 0
            supported_lcc[ti] = supp_sum / denom_cons
            frag_fail[ti] = fail_frag / denom_cons
        else
            supported_lcc[ti] = 0.0
            frag_fail[ti] = 0.0
        end
    end

    if ENFORCE_MONOTONE
        monotone_check_and_fix!(rA_E); monotone_check_and_fix!(rB_E)
        monotone_check_and_fix!(rA_V); monotone_check_and_fix!(rB_V)
    end

    return (rA_E=rA_E, rB_E=rB_E, rA_V=rA_V, rB_V=rB_V,
            mean_phi_req=mean_phi_req, supported_lcc=supported_lcc, frag_fail=frag_fail,
            creal=creal)
end

# ----------------------------- AUC METRICS ------------------------------------
function auc_metrics(rA::Vector{Float64}, rB::Vector{Float64})
    aucA = auc_trapz(collect(FGRID), rA)
    aucB = auc_trapz(collect(FGRID), rB)
    d = rA .- rB
    aucD = auc_trapz(collect(FGRID), d)
    rel = aucA > 1e-9 ? (aucD / aucA) : 0.0
    return (aucA=aucA, aucB=aucB, dAUC_raw=aucD, dAUC_rel=rel)
end

# ----------------------------- HEATMAP UPSAMPLING -----------------------------
# Bilinear interpolation of Z defined on (xcoarse, ycoarse) to (xfine, yfine)
function bilinear_upsample(xc::Vector{Float64}, yc::Vector{Float64}, Zc::Matrix{Float64},
                           xf::Vector{Float64}, yf::Vector{Float64})
    Zf = fill(NaN, length(yf), length(xf))

    # helper to find bracketing index
    function bracket(v::Vector{Float64}, x::Float64)
        if x <= v[1]; return 1, 1, 0.0 end
        if x >= v[end]; return length(v), length(v), 0.0 end
        hi = searchsortedfirst(v, x)
        lo = hi - 1
        t = (x - v[lo]) / (v[hi] - v[lo] + 1e-12)
        return lo, hi, t
    end

    for (ix, x) in enumerate(xf)
        x0, x1, tx = bracket(xc, x)
        for (iy, y) in enumerate(yf)
            y0, y1, ty = bracket(yc, y)
            z00 = Zc[y0, x0]
            z10 = Zc[y0, x1]
            z01 = Zc[y1, x0]
            z11 = Zc[y1, x1]
            # if any NaN, propagate conservatively
            if any(isnan, (z00,z10,z01,z11))
                Zf[iy, ix] = NaN
            else
                z0 = (1-tx)*z00 + tx*z10
                z1 = (1-tx)*z01 + tx*z11
                Zf[iy, ix] = (1-ty)*z0 + ty*z1
            end
        end
    end
    return Zf
end

# ----------------------------- FIGURE BUILDERS --------------------------------
function fig1_examples!(results_dict, ex_low, ex_high; metric=:Emin)
    fig = Figure(size=(1450, 720))
    geoms = GEOMS

    function plot_case!(row::Int, label::String, ex)
        sc, corr, pconn = ex.sc, ex.corr, ex.p
        for (j, g) in enumerate(geoms)
            ax = Axis(fig[row, j],
                title = "$(label) — $(g)",
                xlabel = (row == 2 ? "habitat loss f" : ""),
                ylabel = (j == 1 ? "patch richness" : "")
            )
            res = results_dict[(sc.id, corr, pconn, g)]
            if metric == :Emin
                rA, rB = res.rA_E, res.rB_E
            else
                rA, rB = res.rA_V, res.rB_V
            end
            lines!(ax, collect(FGRID), rA, linewidth=3, label="A-only")
            lines!(ax, collect(FGRID), rB, linewidth=3, label="AB")
            xlims!(ax, 0, F_MAX)
            ylims!(ax, 0, max(maximum(rA), maximum(rB)) * 1.05)
            if row == 1 && j == 1
                axislegend(ax; position=:lb, framevisible=true)
            end
        end
    end

    plot_case!(1, "Low divergence", ex_low)
    plot_case!(2, "High divergence", ex_high)

    ttl = metric == :Emin ? "Fig1 — Example behaviours (Emin_patch; connectivity-sensitive)" :
                            "FigS1 — Example behaviours (viability; area-only)"
    Label(fig[0, :], ttl, fontsize=18)
    return fig
end

function finite_colorrange(Z)
    vals = Z[isfinite.(Z)]
    isempty(vals) && return (0.0, 1.0)   # safe fallback
    lo = minimum(vals)
    hi = maximum(vals)
    lo == hi && (lo -= 1e-6; hi += 1e-6)
    return (lo, hi)
end

function fig2_heatmaps(auc_table; which=:Emin)
    scs = SCENARIOS
    nshow = min(4, length(scs))

    fig_raw = Figure(size=(1500, 720))
    fig_rel = Figure(size=(1500, 720))

    for (col, sc) in enumerate(scs[1:nshow])
        xc = CORR_GRID_COARSE
        yc = PCONN_GRID_COARSE

        Zraw_c = fill(NaN, length(yc), length(xc))
        Zrel_c = fill(NaN, length(yc), length(xc))

        for (xi, corr) in enumerate(xc), (yi, p) in enumerate(yc)
            key = (sc.id, corr, p, HEATMAP_GEOM, which)
            if haskey(auc_table, key)
                Zraw_c[yi, xi] = auc_table[key].dAUC_raw
                Zrel_c[yi, xi] = auc_table[key].dAUC_rel
            end
        end

        # upsample to fine grid for smooth visual transitions
        Zraw = bilinear_upsample(xc, yc, Zraw_c, CORR_GRID_FINE, PCONN_GRID_FINE)
        Zrel = bilinear_upsample(xc, yc, Zrel_c, CORR_GRID_FINE, PCONN_GRID_FINE)

        ax1 = Axis(fig_raw[1, col],
            title = "$(sc.name) — $(HEATMAP_GEOM)",
            xlabel = "corr (link assortativity)",
            ylabel = (col == 1 ? "target connectance p" : ""),
            titlesize = 9
        )
        cr = finite_colorrange(Zraw)
        hm1 = image!(ax1, CORR_GRID_FINE, PCONN_GRID_FINE, Zraw; colorrange = cr, colormap = :viridis)
        xlims!(ax1, minimum(CORR_GRID_FINE), maximum(CORR_GRID_FINE))
        ylims!(ax1, minimum(PCONN_GRID_FINE), maximum(PCONN_GRID_FINE))

        ax2 = Axis(fig_rel[1, col],
            title = "$(sc.name) — $(HEATMAP_GEOM)",
            xlabel = "corr (link assortativity)",
            ylabel = (col == 1 ? "target connectance p" : ""),
            titlesize = 9
        )
        cr = finite_colorrange(Zrel)
        hm2 = image!(ax2, CORR_GRID_FINE, PCONN_GRID_FINE, Zrel; colorrange = cr, colormap = :viridis)
        xlims!(ax2, minimum(CORR_GRID_FINE), maximum(CORR_GRID_FINE))
        ylims!(ax2, minimum(PCONN_GRID_FINE), maximum(PCONN_GRID_FINE))

        if col == nshow
            Colorbar(fig_raw[1, col+1], hm1, label="raw ΔAUC = ∫(A-only − AB) df")
            Colorbar(fig_rel[1, col+1], hm2, label="relative ΔAUC = raw ΔAUC / AUC(A-only)")
        end
    end

    ttl_raw = which == :Emin ? "Fig2a — Raw ΔAUC (Emin_patch)" : "FigS2a — Raw ΔAUC (viability)"
    ttl_rel = which == :Emin ? "Fig2b — Relative ΔAUC (Emin_patch)" : "FigS2b — Relative ΔAUC (viability)"
    Label(fig_raw[0, :], ttl_raw, fontsize=18)
    Label(fig_rel[0, :], ttl_rel, fontsize=18)
    return fig_raw, fig_rel
end

function fig3_mechanism(ex_high, results_dict)
    fig = Figure(size=(1450, 360))
    sc, corr, p = ex_high.sc, ex_high.corr, ex_high.p

    for (j, g) in enumerate(GEOMS)
        ax = Axis(fig[1, j],
            title="Mechanism: support amount vs support connectivity — $(g)",
            xlabel="habitat loss f",
            ylabel="mean φ (consumers)"
        )
        axr = Axis(fig[1, j], yaxisposition=:right, ylabel="supported-LCC / frag-fail")
        hidespines!(axr); hidexdecorations!(axr)

        res = results_dict[(sc.id, corr, p, g)]
        lines!(ax, collect(FGRID), res.mean_phi_req, linewidth=3, label="mean φ")
        lines!(axr, collect(FGRID), res.supported_lcc, linewidth=3, linestyle=:dash, label="supported-LCC")
        lines!(axr, collect(FGRID), res.frag_fail, linewidth=3, linestyle=:dot, label="frag-fail")

        xlims!(ax, 0, F_MAX)
        ylims!(ax, 0, 1)
        ylims!(axr, -0.1, 1.1)

        if j == 1
            axislegend(ax; position=:lb, framevisible=true)
        end
    end
    Label(fig[0, :], "Fig3 (test) — Mechanism panel", fontsize=18)
    return fig
end

# ----------------------------- RUN EVERYTHING ---------------------------------
@info "Generating synthetic environment + precomputing remain masks..."
E = make_environment()
remain_masks = precompute_remain_masks(; rng=MersenneTwister(SEED + 999))

# Example points:
# Low divergence: broad niches + high corr + higher connectance
# High divergence: narrow niches + low corr + low connectance
const EX_LOW  = (sc=SCENARIOS[4], corr=1.0, p=0.16)
const EX_HIGH = (sc=SCENARIOS[3], corr=0.0, p=0.04)

@info "Running sweep for Fig2 (coarse grid) and storing example curves for Fig1/Fig3..."

# AUC table:
# auc_table[(scenario_id, corr, pconn, geom, :Emin/:Viab)] -> NamedTuple metrics
auc_table = Dict{Tuple{Int,Float64,Float64,String,Symbol}, NamedTuple}()

# Store full curves only for EX_LOW and EX_HIGH across all geometries
results_dict = Dict{Tuple{Int,Float64,Float64,String}, NamedTuple}()

# Deterministic RNG per (scenario,corr,p,geom) for reproducibility
function local_rng(scid, corr, p, geom)
    h = hash((SEED, scid, round(Int, corr*1000), round(Int, p*10000), geom))
    return MersenneTwister(UInt32(mod(h, typemax(UInt32))))
end

# -------------------------------------------------
# 1) Build job list (pure data, no computation)
# -------------------------------------------------
jobs = [
    (sc, corr, p, geom)
    for sc in SCENARIOS
    for corr in CORR_GRID_COARSE
    for p in PCONN_GRID_COARSE
    for geom in (HEATMAP_GEOM,)
]

# -------------------------------------------------
# 2) Thread-local result containers
# -------------------------------------------------
# each thread writes to its own Dict
thread_tables = [Dict{Tuple, Any}() for _ in 1:Threads.nthreads()]

# -------------------------------------------------
# 3) Parallel loop
# -------------------------------------------------
Threads.@threads for idx in eachindex(jobs)
    tid = Threads.threadid()
    local_table = thread_tables[tid]

    sc, corr, p, geom = jobs[idx]

    rng = local_rng(sc.id, corr, p, geom)

    res = simulate_one(
        E, remain_masks,
        sc, corr, p, geom;
        rng = rng
    )

    local_table[(sc.id, corr, p, geom, :Emin)] =
        auc_metrics(res.rA_E, res.rB_E)

    local_table[(sc.id, corr, p, geom, :Viab)] =
        auc_metrics(res.rA_V, res.rB_V)

end

# -------------------------------------------------
# 4) Merge results (serial, safe)
# -------------------------------------------------
auc_table = Dict{Tuple, Any}()

for tbl in thread_tables
    merge!(auc_table, tbl)
end

# Run and store example curves for Fig1 + Fig3 across all geometries
for ex in (EX_LOW, EX_HIGH)
    for geom in GEOMS
        rng = local_rng(ex.sc.id, ex.corr, ex.p, geom)
        res = simulate_one(E, remain_masks, ex.sc, ex.corr, ex.p, geom; rng=rng)
        results_dict[(ex.sc.id, ex.corr, ex.p, geom)] = res

        # Also store AUC summaries for console reporting
        auc_table[(ex.sc.id, ex.corr, ex.p, geom, :Emin)] = auc_metrics(res.rA_E, res.rB_E)
        auc_table[(ex.sc.id, ex.corr, ex.p, geom, :Viab)] = auc_metrics(res.rA_V, res.rB_V)
    end
end

@info "Building figures..."

# Fig1 (Emin main + viability SI)
fig1 = fig1_examples!(results_dict, EX_LOW, EX_HIGH; metric=:Emin)
display(fig1)
if SAVE_FIGS
    save(joinpath(OUT_DIR, "Fig1_examples_Emin.png"), fig1, px_per_unit=2)
    save(joinpath(OUT_DIR, "Fig1_examples_Emin.pdf"), fig1)
end

figS1 = fig1_examples!(results_dict, EX_LOW, EX_HIGH; metric=:Viab)
display(figS1)
if SAVE_FIGS
    save(joinpath(OUT_DIR, "FigS1_examples_viability.png"), figS1, px_per_unit=2)
    save(joinpath(OUT_DIR, "FigS1_examples_viability.pdf"), figS1)
end

# Fig2 (heatmaps) Emin main + viability SI (same computed coarse grid)
fig2_raw, fig2_rel = fig2_heatmaps(auc_table; which=:Emin)
display(fig2_raw); display(fig2_rel)
if SAVE_FIGS
    save(joinpath(OUT_DIR, "Fig2_raw_dAUC_Emin.png"), fig2_raw, px_per_unit=2)
    save(joinpath(OUT_DIR, "Fig2_raw_dAUC_Emin.pdf"), fig2_raw)
    save(joinpath(OUT_DIR, "Fig2_rel_dAUC_Emin.png"), fig2_rel, px_per_unit=2)
    save(joinpath(OUT_DIR, "Fig2_rel_dAUC_Emin.pdf"), fig2_rel)
end

figS2_raw, figS2_rel = fig2_heatmaps(auc_table; which=:Viab)
display(figS2_raw); display(figS2_rel)
if SAVE_FIGS
    save(joinpath(OUT_DIR, "FigS2_raw_dAUC_viability.png"), figS2_raw, px_per_unit=2)
    save(joinpath(OUT_DIR, "FigS2_raw_dAUC_viability.pdf"), figS2_raw)
    save(joinpath(OUT_DIR, "FigS2_rel_dAUC_viability.png"), figS2_rel, px_per_unit=2)
    save(joinpath(OUT_DIR, "FigS2_rel_dAUC_viability.pdf"), figS2_rel)
end

# Fig3 (mechanism) using EX_HIGH
fig3 = fig3_mechanism(EX_HIGH, results_dict)
display(fig3)
if SAVE_FIGS
    save(joinpath(OUT_DIR, "Fig3_mechanism.png"), fig3, px_per_unit=2)
    save(joinpath(OUT_DIR, "Fig3_mechanism.pdf"), fig3)
end

# Console report for key examples
function report_summary(sc::Scenario, geom::String, corr::Float64, p::Float64)
    mE = auc_table[(sc.id, corr, p, geom, :Emin)]
    mV = auc_table[(sc.id, corr, p, geom, :Viab)]
    @info @sprintf("[%s | %s | corr=%.2f | p=%.3f] Emin: rawΔAUC=%.2f relΔAUC=%.3f | Viab: rawΔAUC=%.2f relΔAUC=%.3f",
        sc.name, geom, corr, p, mE.dAUC_raw, mE.dAUC_rel, mV.dAUC_raw, mV.dAUC_rel)
end

@info "Key example summaries (EX_LOW):"
for g in GEOMS
    report_summary(EX_LOW.sc, g, EX_LOW.corr, EX_LOW.p)
end
@info "Key example summaries (EX_HIGH):"
for g in GEOMS
    report_summary(EX_HIGH.sc, g, EX_HIGH.corr, EX_HIGH.p)
end

@info "DONE. Outputs in: $OUT_DIR"
###############################################################################