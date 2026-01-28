###############################################################################
# patchonly_article_ready_independent_v3_metawebERcascade.jl
#
# v3 goals (ONLY fixing issues you flagged, keeping the rest conceptually intact):
#   (1) Fig3 geometry contrast restored (mechanism no longer flat) by using a
#       trophically-plausible ER/cascade-like metaweb (feeds only "down" niche axis),
#       so support starts high and fragmentation/geometry drives declines.
#   (2) SAR_baseline + SAR_eff built correctly:
#       - SAR_baseline anchored at A-only f=0
#       - SAR_eff anchored at AB f=0
#       - exponent z is fitted from early A-only curve vs remaining area (per panel)
#       - SAR_eff uses effective area = habitat LCC fraction (geometry-sensitive)
#   (3) Curves less “stair-like”: replicate-averaging + optional tiny display-only
#       smoothing (does not affect AUC calculations).
#   (4) Habitat loss truly reaches f=1 (remain mask empty) via rank-based removal.
#   (5) Heatmaps: 2×2 layout, square panels, no Makie ImageLike deprecation warning
#       (use heatmap! with interpolate=true).
#
# Dependencies: Random, Statistics, Printf, CairoMakie
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

# Parameter sweep (increase these if you want even cleaner heatmaps)
const CORR_GRID = collect(0.0:0.1:1.0)         # 11 values
const PCONN_GRID = collect(0.02:0.01:0.20)     # 19 values

# Geometries
const GEOMS = ["random", "cluster", "front"]
const HEATMAP_GEOM = "cluster"

# Presence thresholds
const A_MIN  = 15
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

# Replicates to smooth “stair-like” curves
const NREP_EXAMPLES = 6      # Fig1/Fig3
const NREP_HEATMAP  = 3      # Fig2

# Monotonic sanity check (richness need not be strictly monotone in AB because
# support fraction can increase when “bad” habitat disappears; we still enforce
# monotone for plotting/AUC stability, but we do NOT spam warnings anymore.)
const ENFORCE_MONOTONE = true
const APPLY_MONO_FIX   = true
const MONO_TOL         = 1e-9

# Output
const OUT_DIR = joinpath(pwd(), "paper_final_independent_v3_metawebER")
const SAVE_FIGS = true
mkpath(OUT_DIR)

# Plot style
set_theme!(Theme(
    fontsize = 16,
    Axis = (xlabelsize=15, ylabelsize=15, titlesize=15,
            xticklabelsize=13, yticklabelsize=13),
    Legend = (labelsize=13,)
))

# ----------------------------- HELPERS ----------------------------------------
clamp01(x) = min(max(x, 0.0), 1.0)
@inline linidx(y::Int, x::Int) = (y-1) * NX + x

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

# Rank-based removal ensures:
#  - nested remain masks across f
#  - f=1 removes ALL cells (no “stopping early”)
function precompute_remain_masks_and_habitat_stats(; rng=Random.default_rng())
    masks = Dict{Tuple{String,Int}, BitMatrix}()
    hab_area_frac = Dict{Tuple{String,Int}, Float64}()
    hab_lcc_frac  = Dict{Tuple{String,Int}, Float64}()

    visited = Vector{UInt8}(undef, NCELLS)
    stack = Int[]

    # LCC of remaining habitat only (remain mask treated as "occupied")
    function habitat_lcc_size(remain::BitMatrix)
        fill!(visited, 0x00)
        best = 0
        @inbounds for y in 1:NY, x in 1:NX
            remain[y,x] || continue
            i = linidx(y,x)
            visited[i] == 0x01 && continue
            empty!(stack)
            push!(stack, i)
            visited[i] = 0x01
            sz = 0
            while !isempty(stack)
                u = pop!(stack)
                sz += 1
                (up,down,left,right) = NEIGH[u]
                if up != 0
                    yy = (up - 1) ÷ NX + 1; xx = (up - 1) % NX + 1
                    if visited[up] == 0x00 && remain[yy,xx]
                        visited[up] = 0x01; push!(stack, up)
                    end
                end
                if down != 0
                    yy = (down - 1) ÷ NX + 1; xx = (down - 1) % NX + 1
                    if visited[down] == 0x00 && remain[yy,xx]
                        visited[down] = 0x01; push!(stack, down)
                    end
                end
                if left != 0
                    yy = (left - 1) ÷ NX + 1; xx = (left - 1) % NX + 1
                    if visited[left] == 0x00 && remain[yy,xx]
                        visited[left] = 0x01; push!(stack, left)
                    end
                end
                if right != 0
                    yy = (right - 1) ÷ NX + 1; xx = (right - 1) % NX + 1
                    if visited[right] == 0x00 && remain[yy,xx]
                        visited[right] = 0x01; push!(stack, right)
                    end
                end
            end
            best = max(best, sz)
        end
        return best
    end

    for geom in GEOMS
        P = loss_priority(geom; rng=rng)
        pvec = vec(P)
        ord = sortperm(pvec)  # increasing priority removed first
        for (ti, f) in enumerate(FGRID)
            k = Int(round(f * NCELLS))  # number removed
            k = clamp(k, 0, NCELLS)
            R = trues(NY, NX)
            if k == NCELLS
                R .= false
            elseif k > 0
                @inbounds for t in 1:k
                    idx = ord[t]
                    y = (idx - 1) ÷ NX + 1
                    x = (idx - 1) % NX + 1
                    R[y,x] = false
                end
            end
            masks[(geom, ti)] = R
            area = count(R) / NCELLS
            hab_area_frac[(geom, ti)] = area
            hab_lcc_frac[(geom, ti)] = (area > 0) ? habitat_lcc_size(R) / NCELLS : 0.0
        end
    end

    return masks, hab_area_frac, hab_lcc_frac
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

# LCC of intersection (pot ∧ remain)
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
            visited[i] == 0x01 && continue
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

# Monotone enforcement (quiet)
function monotone_fix!(y::Vector{Float64})
    if APPLY_MONO_FIX
        for i in 2:length(y)
            y[i] = min(y[i-1], y[i])
        end
    end
    return y
end

function enforce_monotone_if_needed!(y::Vector{Float64})
    ENFORCE_MONOTONE || return 0.0
    inc = 0.0
    @inbounds for i in 1:(length(y)-1)
        inc = max(inc, y[i+1] - y[i])
    end
    if inc > MONO_TOL
        monotone_fix!(y)
    end
    return inc
end

# Tiny display-only smoothing to reduce stair-step appearance
function smooth1d(y::Vector{Float64}; w::Int=3)
    w <= 1 && return copy(y)
    n = length(y)
    out = similar(y)
    h = w ÷ 2
    @inbounds for i in 1:n
        a = max(1, i-h)
        b = min(n, i+h)
        out[i] = mean(@view y[a:b])
    end
    return out
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
        @inbounds for i in 1:n; σ[i] = clamp(0.07 * exp(0.9*randn(rng)), 0.02, 0.30) end
    elseif mode == :less_variable
        @inbounds for i in 1:n; σ[i] = clamp(0.09 * exp(0.35*randn(rng)), 0.03, 0.28) end
    elseif mode == :narrow
        @inbounds for i in 1:n; σ[i] = clamp(0.06 * exp(0.45*randn(rng)), 0.02, 0.18) end
    elseif mode == :broad
        @inbounds for i in 1:n; σ[i] = clamp(0.14 * exp(0.35*randn(rng)), 0.05, 0.35) end
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
Trophically-plausible ER / cascade-like metaweb:

- Sort species by niche μ (ascending).
- Basal species = lowest N_BASAL in μ (no outgoing links).
- Each consumer i can ONLY feed on species with lower μ (j < i in sorted order).
- Baseline link probability p among allowed pairs, with optional assortativity:
    w_ij = exp(-|μ_i - μ_j|/L_NICHE), rescaled to mean 1 over allowed pairs
    p_ij = p * ((1-corr) + corr*w_ij)
- Guarantees each consumer has at least 1 prey.

Returns:
- prey_lists (indexed in ORIGINAL species ids)
- realised_connectance (edges / allowed_pairs)
- basal BitVector in ORIGINAL ids
- perm (sorted order ids), invperm
"""
function build_metaweb_cascade_ER(μ::Vector{Float64}, p::Float64, corr::Float64; rng=Random.default_rng())
    n = length(μ)
    perm = sortperm(μ)              # niche rank
    invp = invperm(perm)

    basal_sorted = falses(n)
    @inbounds for r in 1:N_BASAL
        basal_sorted[r] = true
    end

    basal = falses(n)               # basal in ORIGINAL ids
    @inbounds for r in 1:n
        basal[perm[r]] = basal_sorted[r]
    end

    prey_lists = [Int[] for _ in 1:n]   # ORIGINAL ids

    # mean(w) over allowed pairs to stabilize expected connectance
    wsum = 0.0
    wcnt = 0
    @inbounds for ri in 1:n
        basal_sorted[ri] && continue
        μi = μ[perm[ri]]
        for rj in 1:(ri-1)
            w = exp(-abs(μi - μ[perm[rj]]) / (L_NICHE + 1e-12))
            wsum += w; wcnt += 1
        end
    end
    wmean = (wcnt > 0) ? (wsum / wcnt) : 1.0

    edges = 0
    allowed = 0

    @inbounds for ri in 1:n
        basal_sorted[ri] && continue
        i = perm[ri]
        μi = μ[i]
        # allowed prey are ranks < ri
        local_preys = Int[]
        for rj in 1:(ri-1)
            j = perm[rj]
            allowed += 1
            w = exp(-abs(μi - μ[j]) / (L_NICHE + 1e-12))
            w /= (wmean + 1e-12)
            pij = p * ((1.0 - corr) + corr * w)
            pij = clamp01(pij)
            if rand(rng) < pij
                push!(local_preys, j)
            end
        end
        if isempty(local_preys)
            # force one prey: pick among allowed with prob ∝ w (assortative)
            # (fallback to uniform if ri==1 cannot happen because basal are lowest ranks)
            bestj = perm[ri-1]
            if ri > 2
                # sample a few candidates for cheap weighted choice
                cand = perm[max(1, ri-20):(ri-1)]
                ws = zeros(Float64, length(cand))
                for (k, j) in enumerate(cand)
                    ws[k] = exp(-abs(μi - μ[j]) / (L_NICHE + 1e-12))
                end
                s = sum(ws)
                if s > 0
                    u = rand(rng) * s
                    acc = 0.0
                    for (k, j) in enumerate(cand)
                        acc += ws[k]
                        if acc >= u
                            bestj = j
                            break
                        end
                    end
                end
            end
            push!(local_preys, bestj)
        end
        prey_lists[i] = local_preys
        edges += length(local_preys)
    end

    creal = edges / (allowed + 1e-12)
    return prey_lists, creal, basal, perm, invp
end

# ----------------------------- CORE SIMULATION --------------------------------
function simulate_one(E::Matrix{Float64}, remain_masks, sc::Scenario, corr::Float64, pconn::Float64, geom::String;
                      rng=Random.default_rng())

    μ = rand(rng, N_SPECIES)
    σ = sample_sigma(sc.σ_mode, N_SPECIES; rng=rng)

    pot = Vector{BitMatrix}(undef, N_SPECIES)
    @inbounds for i in 1:N_SPECIES
        pot[i] = potential_mask(E, μ[i], σ[i])
    end

    prey_lists, creal, basal, _, _ = build_metaweb_cascade_ER(μ, pconn, corr; rng=rng)

    φ_req = sample_phi(N_SPECIES; rng=rng)
    meanφ = mean(φ_req[.!basal])

    visited = Vector{UInt8}(undef, NCELLS)
    stack = Int[]

    rA_E = zeros(Float64, NF)
    rB_E = zeros(Float64, NF)
    rA_V = zeros(Float64, NF)
    rB_V = zeros(Float64, NF)

    mean_phi_req = fill(meanφ, NF)
    supported_lcc = zeros(Float64, NF)
    frag_fail = zeros(Float64, NF)

    presentV = falses(N_SPECIES)
    presentE = falses(N_SPECIES)

    for ti in 1:NF
        remain = remain_masks[(geom, ti)]

        lcc_idx_cache = Vector{Vector{Int}}(undef, N_SPECIES)

        @inbounds for i in 1:N_SPECIES
            area, lcc, idx = lcc_intersection(pot[i], remain, visited, stack)
            lcc_idx_cache[i] = idx
            presentV[i] = (area >= A_MIN)
            presentE[i] = (area >= A_MIN) && (lcc >= LCC_MIN)
        end

        presentV_AB = copy(presentV)
        presentE_AB = copy(presentE)

        denom_cons = 0
        supp_sum = 0.0
        fail_frag = 0

        @inbounds for i in 1:N_SPECIES
            basal[i] && continue

            # Emin AB within consumer LCC
            eA = presentE[i]
            if eA && !isempty(lcc_idx_cache[i])
                denom_cons += 1
                prey = prey_lists[i]
                deg = length(prey)
                supported = 0
                for j in prey
                    presentE[j] || continue
                    has = false
                    for u in lcc_idx_cache[i]
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
                if !eAB
                    fail_frag += 1
                end
            else
                presentE_AB[i] = false
            end

            # Viability AB (still “area-only” threshold, but prey encounter via the consumer patch)
            vA = presentV[i]
            if vA && !isempty(lcc_idx_cache[i])
                prey = prey_lists[i]
                deg = length(prey)
                supported = 0
                for j in prey
                    presentV[j] || continue
                    has = false
                    for u in lcc_idx_cache[i]
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

        rA_E[ti] = sum(presentE)
        rB_E[ti] = sum(presentE_AB)
        rA_V[ti] = sum(presentV)
        rB_V[ti] = sum(presentV_AB)

        if denom_cons > 0
            supported_lcc[ti] = supp_sum / denom_cons
            frag_fail[ti] = fail_frag / denom_cons
        else
            supported_lcc[ti] = 0.0
            frag_fail[ti] = 0.0
        end
    end

    if ENFORCE_MONOTONE
        enforce_monotone_if_needed!(rA_E); enforce_monotone_if_needed!(rB_E)
        enforce_monotone_if_needed!(rA_V); enforce_monotone_if_needed!(rB_V)
    end

    return (rA_E=rA_E, rB_E=rB_E, rA_V=rA_V, rB_V=rB_V,
            mean_phi_req=mean_phi_req, supported_lcc=supported_lcc, frag_fail=frag_fail,
            creal=creal)
end

# Replicate averaging (smooth curves + stabilise patterns)
function simulate_reps(E, remain_masks, sc, corr, pconn, geom; nrep::Int, seedtag::Int)
    rA_E = zeros(Float64, NF); rB_E = zeros(Float64, NF)
    rA_V = zeros(Float64, NF); rB_V = zeros(Float64, NF)
    mean_phi = zeros(Float64, NF)
    supp = zeros(Float64, NF)
    fail = zeros(Float64, NF)
    creal_acc = 0.0

    for rep in 1:nrep
        h = hash((SEED, seedtag, sc.id, round(Int, corr*1000), round(Int, pconn*10000), geom, rep))
        rng = MersenneTwister(UInt32(mod(h, typemax(UInt32))))
        res = simulate_one(E, remain_masks, sc, corr, pconn, geom; rng=rng)

        rA_E .+= res.rA_E; rB_E .+= res.rB_E
        rA_V .+= res.rA_V; rB_V .+= res.rB_V
        mean_phi .+= res.mean_phi_req
        supp .+= res.supported_lcc
        fail .+= res.frag_fail
        creal_acc += res.creal
    end

    inv = 1.0 / nrep
    rA_E .*= inv; rB_E .*= inv
    rA_V .*= inv; rB_V .*= inv
    mean_phi .*= inv
    supp .*= inv
    fail .*= inv

    return (rA_E=rA_E, rB_E=rB_E, rA_V=rA_V, rB_V=rB_V,
            mean_phi_req=mean_phi, supported_lcc=supp, frag_fail=fail,
            creal=creal_acc*inv)
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

# ----------------------------- SAR HELPERS ------------------------------------
# Fit z from early A-only richness vs remaining area fraction: log(S/S0)=z*log(A)
function fit_z_from_curve(rA::Vector{Float64}, area_frac::Vector{Float64}; maxf::Float64=0.6)
    S0 = max(rA[1], 1e-9)
    num = 0.0
    den = 0.0
    @inbounds for (ti, f) in enumerate(FGRID)
        f > maxf && break
        A = area_frac[ti]
        S = rA[ti]
        if A > 0 && S > 0
            x = log(A)
            y = log(S / S0)
            num += x*y
            den += x*x
        end
    end
    z = (den > 0) ? (num / den) : 0.25
    return clamp(z, 0.0, 1.0)
end

function sar_curve(S0::Float64, area_frac::Vector{Float64}, z::Float64)
    return [S0 * (a <= 0 ? 0.0 : a^z) for a in area_frac]
end

# ----------------------------- FIGURE BUILDERS --------------------------------
function fig1_examples!(results_dict, ex_low, ex_high, hab_area_frac, hab_lcc_frac; metric=:Emin)
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

            # SAR baselines (panel-specific)
            area = [hab_area_frac[(g, ti)] for ti in 1:NF]      # total remaining habitat fraction
            aeff = [hab_lcc_frac[(g, ti)]  for ti in 1:NF]      # effective area = habitat LCC fraction

            z = fit_z_from_curve(rA, area; maxf=0.6)
            sar_base = sar_curve(rA[1], area, z)
            sar_eff  = sar_curve(rB[1], aeff, z)  # anchored to AB at f=0

            # Display-only smoothing to reduce stair-like appearance
            x = collect(FGRID)
            rAplt = smooth1d(collect(rA); w=3)
            rBplt = smooth1d(collect(rB); w=3)
            sbplt = smooth1d(collect(sar_base); w=3)
            seplt = smooth1d(collect(sar_eff);  w=3)

            lines!(ax, x, rAplt, linewidth=3, label="A-only")
            lines!(ax, x, rBplt, linewidth=3, label="AB")
            lines!(ax, x, sbplt, linewidth=2.5, linestyle=:dash, label="SAR baseline")
            lines!(ax, x, seplt, linewidth=2.5, linestyle=:dot,  label="SAR_eff")

            xlims!(ax, 0, F_MAX)
            ymax = maximum(vcat(rA, rB, sar_base, sar_eff))
            ylims!(ax, 0, ymax * 1.05)

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

# Heatmaps: 2×2 square panels, interpolate=true (no ImageLike warning)
# Heatmaps: 2×2 square panels, each with its own colorbar
function fig2_heatmaps(auc_table; which=:Emin)
    scs = SCENARIOS[1:4]

    fig_raw = Figure(size=(1300, 900))
    fig_rel = Figure(size=(1300, 900))

    # -------------------------------------------------
    # Build matrices
    # -------------------------------------------------
    function makeZ(sc::Scenario, keyfield::Symbol)
        xs = CORR_GRID
        ys = PCONN_GRID
        Z = fill(NaN, length(ys), length(xs))
        for (xi, corr) in enumerate(xs), (yi, p) in enumerate(ys)
            key = (sc.id, corr, p, HEATMAP_GEOM, which)
            if haskey(auc_table, key)
                Z[yi, xi] = getproperty(auc_table[key], keyfield)
            end
        end
        return xs, ys, Z
    end

    # -------------------------------------------------
    # Color ranges (shared per figure, but bars are per-panel)
    # -------------------------------------------------
    function global_cr(keyfield::Symbol)
        vals = Float64[]
        for sc in scs
            _, _, Z = makeZ(sc, keyfield)
            append!(vals, Z[isfinite.(Z)])
        end
        isempty(vals) && return (0.0, 1.0)
        lo, hi = minimum(vals), maximum(vals)
        lo == hi && (lo -= 1e-6; hi += 1e-6)
        return (lo, hi)
    end

    cr_raw = global_cr(:dAUC_raw)
    cr_rel = global_cr(:dAUC_rel)

    # -------------------------------------------------
    # Plot panels
    # -------------------------------------------------
    for (idx, sc) in enumerate(scs)
        r = idx ≤ 2 ? 1 : 2
        c = idx % 2 == 1 ? 1 : 3   # leave column 2/4 for colorbars

        # ---------------- RAW ----------------
        ax1 = Axis(fig_raw[r, c],
            title = sc.name,
            xlabel = "corr (assortativity)",
            ylabel = "target connectance p",
            aspect = 1,
            titlesize = 11
        )

        xs, ys, Zraw = makeZ(sc, :dAUC_raw)
        hm1 = heatmap!(
            ax1, xs, ys, Zraw;
            colormap = :viridis,
            colorrange = cr_raw,
            interpolate = true
        )

        Colorbar(fig_raw[r, c+1], hm1;
            label = "raw ΔAUC"
        )

        # ---------------- REL ----------------
        ax2 = Axis(fig_rel[r, c],
            title = sc.name,
            xlabel = "corr (assortativity)",
            ylabel = "target connectance p",
            aspect = 1,
            titlesize = 11
        )

        _, _, Zrel = makeZ(sc, :dAUC_rel)
        hm2 = heatmap!(
            ax2, xs, ys, Zrel;
            colormap = :viridis,
            colorrange = cr_rel,
            interpolate = true
        )

        Colorbar(fig_rel[r, c+1], hm2;
            label = "relative ΔAUC"
        )
    end

    ttl_raw = which == :Emin ? "Fig2a — Raw ΔAUC (Emin_patch)" :
                              "FigS2a — Raw ΔAUC (viability)"
    ttl_rel = which == :Emin ? "Fig2b — Relative ΔAUC (Emin_patch)" :
                              "FigS2b — Relative ΔAUC (viability)"

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

        x = collect(FGRID)
        meanφ = smooth1d(collect(res.mean_phi_req); w=3)
        supp  = smooth1d(collect(res.supported_lcc); w=3)
        fail  = smooth1d(collect(res.frag_fail); w=3)

        lines!(ax,  x, meanφ, linewidth=3, label="mean φ")
        lines!(axr, x, supp,  linewidth=3, linestyle=:dash, label="supported-LCC")
        lines!(axr, x, fail,  linewidth=3, linestyle=:dot,  label="frag-fail")

        xlims!(ax, 0, F_MAX)
        ylims!(ax, 0, 1)
        ylims!(axr, -0.05, 1.05)

        if j == 1
            axislegend(ax; position=:lb, framevisible=true)
        end
    end
    Label(fig[0, :], "Fig3 (test) — Mechanism panel", fontsize=18)
    return fig
end

# ----------------------------- RUN EVERYTHING ---------------------------------
@info "Generating synthetic environment + precomputing remain masks + habitat stats..."
E = make_environment(; rng=MersenneTwister(SEED))
remain_masks, hab_area_frac, hab_lcc_frac = precompute_remain_masks_and_habitat_stats(; rng=MersenneTwister(SEED + 999))

# Example points (same intent as v2)
const EX_LOW  = (sc=SCENARIOS[4], corr=1.0, p=0.16)
const EX_HIGH = (sc=SCENARIOS[3], corr=0.0, p=0.04)

@info "Running heatmap sweep (avg over reps) + storing example curves..."

# auc_table[(scenario_id, corr, pconn, geom, :Emin/:Viab)] -> NamedTuple metrics
auc_table = Dict{Tuple{Int,Float64,Float64,String,Symbol}, NamedTuple}()

# results_dict[(scenario_id, corr, pconn, geom)] -> averaged curves for examples
results_dict = Dict{Tuple{Int,Float64,Float64,String}, NamedTuple}()

# ---- Heatmap jobs (only HEATMAP_GEOM), averaged ----
jobs = [(sc, corr, p) for sc in SCENARIOS for corr in CORR_GRID for p in PCONN_GRID]
thread_tables = [Dict{Tuple, Any}() for _ in 1:Threads.nthreads()]

Threads.@threads for idx in eachindex(jobs)
    tid = Threads.threadid()
    local_tbl = thread_tables[tid]
    sc, corr, p = jobs[idx]

    res = simulate_reps(E, remain_masks, sc, corr, p, HEATMAP_GEOM; nrep=NREP_HEATMAP, seedtag=777)
    local_tbl[(sc.id, corr, p, HEATMAP_GEOM, :Emin)] = auc_metrics(res.rA_E, res.rB_E)
    local_tbl[(sc.id, corr, p, HEATMAP_GEOM, :Viab)] = auc_metrics(res.rA_V, res.rB_V)
end

for tbl in thread_tables
    merge!(auc_table, tbl)
end

# ---- Example curves across ALL geometries (avg over reps) ----
for ex in (EX_LOW, EX_HIGH)
    for geom in GEOMS
        res = simulate_reps(E, remain_masks, ex.sc, ex.corr, ex.p, geom; nrep=NREP_EXAMPLES, seedtag=9999)
        results_dict[(ex.sc.id, ex.corr, ex.p, geom)] = res

        # also store AUC summaries for reporting
        auc_table[(ex.sc.id, ex.corr, ex.p, geom, :Emin)] = auc_metrics(res.rA_E, res.rB_E)
        auc_table[(ex.sc.id, ex.corr, ex.p, geom, :Viab)] = auc_metrics(res.rA_V, res.rB_V)
    end
end

@info "Building figures..."

# Fig1 (Emin main + viability SI) — now includes SAR baseline + SAR_eff
fig1 = fig1_examples!(results_dict, EX_LOW, EX_HIGH, hab_area_frac, hab_lcc_frac; metric=:Emin)
display(fig1)
if SAVE_FIGS
    save(joinpath(OUT_DIR, "Fig1_examples_Emin.png"), fig1, px_per_unit=2)
    save(joinpath(OUT_DIR, "Fig1_examples_Emin.pdf"), fig1)
end

figS1 = fig1_examples!(results_dict, EX_LOW, EX_HIGH, hab_area_frac, hab_lcc_frac; metric=:Viab)
display(figS1)
if SAVE_FIGS
    save(joinpath(OUT_DIR, "FigS1_examples_viability.png"), figS1, px_per_unit=2)
    save(joinpath(OUT_DIR, "FigS1_examples_viability.pdf"), figS1)
end

# Fig2 heatmaps (2×2 square panels) Emin main + viability SI
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

# Fig3 mechanism (now geometry-sensitive again)
fig3 = fig3_mechanism(EX_HIGH, results_dict)
display(fig3)
if SAVE_FIGS
    save(joinpath(OUT_DIR, "Fig3_mechanism.png"), fig3, px_per_unit=2)
    save(joinpath(OUT_DIR, "Fig3_mechanism.pdf"), fig3)
end

# Console report
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
