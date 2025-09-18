module Metawebs

using Random, Statistics, StatsBase, Distributions

export SpeciesPool, build_metaweb_archetype, metaweb_diagnostics, rewire_placebo
export build_metaweb_niche, build_metaweb_modular, build_metaweb_powerlaw

struct SpeciesPool
    S::Int
    masses::Vector{Float64}
    basal::BitVector          # length S
    prey::Vector{Vector{Int}} # prey lists (indices), basal: empty
end

"Compute quick diagnostics used for acceptance: connectance, median diet, basal share."
function metaweb_diagnostics(pool::SpeciesPool)
    S = pool.S
    E = sum(length.(pool.prey))
    C = E / max(S*(S-1)/2,1)             # directed on lower trophic side
    diets = [length(pool.prey[s]) for s in 1:S if !pool.basal[s]]
    meddiet = isempty(diets) ? 0.0 : median(diets)
    basal_share = count(pool.basal)/S
    (; C, meddiet, basal_share)
end

"Draw prey from a Gaussian kernel on log mass ratios; enforce cap."
function _draw_prey!(rng, pool::SpeciesPool; R0_mean=100.0, sigma=0.6, pmax=0.9,
                     density=0.12, max_prey::Int=9999)
    S = pool.S
    logR0 = log(R0_mean)
    for s in 1:S
        pool.basal[s] && (pool.prey[s] = Int[]; continue)
        plist = Int[]
        for q in 1:S
            (q==s || pool.masses[q] >= pool.masses[s]) && continue
            z = (log(pool.masses[s]/pool.masses[q]) - logR0)/sigma
            p = pmax * exp(-0.5*z^2) * density
            if rand(rng) < p
                push!(plist, q)
            end
        end
        # ensure at least one prey if possible
        if isempty(plist)
            if s > 1
                cand = 1:(s-1)
                target = log(pool.masses[s]) - logR0
                q = cand[argmin(abs.(log.(pool.masses[cand]) .- target))]
                push!(plist, q)
            else
                # no smaller species exist; leave empty
                plist = Int[]
            end
        end
        if length(plist) > max_prey
            shuffle!(rng, plist)
            plist = plist[1:max_prey]
        end
        pool.prey[s] = plist
    end
    pool
end

"Generate masses + basal; then build prey sets to match redundancy archetype."
function build_metaweb_archetype(rng::AbstractRNG; S::Int=150, basal_frac::Float64=0.45,
        archetype::Symbol=:mid)
    # masses (log-uniform)
    logm = collect(range(log(1e-2), log(10.0); length=S)); shuffle!(rng, logm)
    masses = exp.(logm)
    # basal set
    order  = sortperm(masses)
    nB     = clamp(round(Int, basal_frac*S), 0, S)
    basal  = falses(S); basal[order[1:nB]] .= true
    pool = SpeciesPool(S, masses, basal, [Int[] for _ in 1:S])

    # redundancy settings
    if archetype === :high
        density, max_prey = 0.20, 8
    elseif archetype === :mid
        density, max_prey = 0.12, 5
    elseif archetype === :low
        density, max_prey = 0.07, 2
    else
        error("unknown archetype $archetype")
    end

    # draw prey sets
    _draw_prey!(rng, pool; density=density, max_prey=max_prey)
    pool
end

"Degree-preserving placebo rewiring among consumers (simple edge swaps)."
function rewire_placebo(rng::AbstractRNG, pool::SpeciesPool; nswap::Int=5000)
    S = pool.S
    prey = deepcopy(pool.prey)
    consumers = [s for s in 1:S if !pool.basal[s] && length(prey[s])>=1]
    for _ in 1:nswap
        s1, s2 = rand(rng, consumers, 2)
        isempty(prey[s1]) && continue
        isempty(prey[s2]) && continue
        q1 = rand(rng, prey[s1])
        q2 = rand(rng, prey[s2])
        if (q2 ∉ prey[s1]) && (q1 ∉ prey[s2]) && (q1 != s2) && (q2 != s1)
            # swap
            i1 = findfirst(==(q1), prey[s1]); prey[s1][i1] = q2
            i2 = findfirst(==(q2), prey[s2]); prey[s2][i2] = q1
        end
    end
    SpeciesPool(pool.S, copy(pool.masses), copy(pool.basal), prey)
end

# =============== helpers (local to this file) =================================
# build a pool from masses + prey lists; basal are those with 0 prey
function _pool_from_prey(masses, prey; basal_frac=0.35)
    S = length(masses)
    order = sortperm(masses)
    nB = clamp(round(Int, basal_frac*S), 1, S)  # at least one basal
    basal = falses(S)
    basal[order[1:nB]] .= true
    return SpeciesPool(S, masses, basal, prey)
end

# always respect trophic ordering: prey indices must be < predator index
# (we sort species by mass ascending and work in that order)
function _sorted_masses(rng::AbstractRNG, S::Int)
    logm = collect(range(log(1e-2), log(10.0); length=S))
    shuffle!(rng, logm)
    masses = exp.(logm)
    sort(masses)  # ascending
end

# =============== 1) Classic niche-model-ish generator =========================
"""
    build_metaweb_niche(rng; S=175, Ctarget=0.10, beta=1.5)

Classic niche-like interval model on [0,1].
- Each species draws a niche value `n ~ U(0,1)` (we set it proportional to rank/S).
- Feeding range `r ~ Beta(1,beta) * n` (controls connectance).
- Center `c ~ U(r/2, n)`; prey are those with niche `nq` in `[c-r/2, c+r/2]` and lower rank.
`Ctarget` is only a rough guide via `beta`.
"""
function build_metaweb_niche(rng::AbstractRNG; S::Int=175, Ctarget::Float64=0.10, beta::Float64=1.5)
    masses = _sorted_masses(rng, S)
    # niche coordinate proportional to rank (intervality)
    n = collect(range(0, 1; length=S+2))[2:end-1]  # avoid exact 0/1
    prey = [Int[] for _ in 1:S]
    # simple mapping: larger beta -> narrower ranges -> lower connectance
    β = max(beta, 1e-3)
    for s in 2:S
        ns = n[s]
        r  = rand(rng, Beta(1, β)) * ns
        c  = rand(rng) * (ns - r/2) + r/2
        lo, hi = c - r/2, c + r/2
        for q in 1:s-1
            if (n[q] ≥ lo) & (n[q] ≤ hi)
                push!(prey[s], q)
            end
        end
        isempty(prey[s]) && push!(prey[s], rand(rng, 1:s-1))   # ensure ≥1 prey
    end
    _pool_from_prey(masses, prey)
end

# =============== 2) Simple modular (blocky) web ===============================
"""
    build_metaweb_modular(rng; S=175, K=3, p_in=0.25, p_out=0.03)

Split the mass-ordered species into `K` contiguous modules.
Consumers connect to lower-rank prey within-module with prob `p_in`,
and to other modules with prob `p_out` (both capped to lower ranks).
Gives tunable modularity with a bottom-heavy trophic order.
"""
function build_metaweb_modular(rng::AbstractRNG; S::Int=175, K::Int=3, p_in::Float64=0.25, p_out::Float64=0.03)
    masses = _sorted_masses(rng, S)
    prey   = [Int[] for _ in 1:S]
    # module boundaries (approximately equal sizes)
    edges = round.(Int, range(1, S; length=K+1))
    modid = zeros(Int, S)
    for k in 1:K
        for i in edges[k]:min(edges[k+1], S)
            modid[i] = k
        end
    end
    for s in 2:S
        for q in 1:s-1
            pin = (modid[s] == modid[q]) ? p_in : p_out
            if rand(rng) < pin
                push!(prey[s], q)
            end
        end
        isempty(prey[s]) && push!(prey[s], rand(rng, 1:s-1))
    end
    _pool_from_prey(masses, prey)
end

# =============== 3) Power-law diet heterogeneity ==============================
"""
    build_metaweb_powerlaw(rng; S=175, kmin=1, kmax=8, alpha=2.4)

Assign each consumer a diet size `k ~ Zipf(α)` truncated to [kmin,kmax],
then choose that many lower-rank prey uniformly without replacement.
Captures heavy-tailed diet heterogeneity without changing trophic order.
"""
function build_metaweb_powerlaw(rng::AbstractRNG; S::Int=175, kmin::Int=1, kmax::Int=8, alpha::Float64=2.4)
    masses = _sorted_masses(rng, S)
    prey   = [Int[] for _ in 1:S]
    for s in 2:S
        support = max(1, min(kmax, s-1))
        k = clamp(Int(floor((rand(rng)^(-1/(alpha-1))))), kmin, support)
        k = max(kmin, min(k, support))
        cand = randperm(rng, s-1)
        prey[s] = cand[1:k]
        isempty(prey[s]) && push!(prey[s], s-1)
    end
    _pool_from_prey(masses, prey)
end



end # module
