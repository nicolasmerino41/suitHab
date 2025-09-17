module Metawebs

using Random, Statistics, StatsBase, Distributions

export SpeciesPool, build_metaweb_archetype, metaweb_diagnostics, rewire_placebo

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

end # module
