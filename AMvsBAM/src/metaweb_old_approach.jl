module MetaWeb

using Random, Statistics, StatsBase, Distributions

export Metaweb, build_metaweb

"""
    Metaweb(S, basal_frac, A, trophic_role)

Container for a directed feeding network on S species.
- `A` is S×S adjacency (A[p,q]==1 if predator p eats prey q).
- `trophic_role` is `:basal` or `:consumer` per species.
"""
struct Metaweb
    S::Int
    basal_frac::Float64
    A::BitMatrix
    trophic_role::Vector{Symbol}
end

"""
    build_metaweb(rng; S, basal_frac, connectance, R95, motif_mix)

Build a metaweb with approximate **connectance** (C),
approximate diet breadth (via predator out-degree) with target 95th percentile `R95`,
and a coarse **motif mix**:
- `:chains` (fewer shared prey),
- `:mixed`,
- `:omnivory` (more shared prey / cross-level links).
"""
function build_metaweb(rng::AbstractRNG; S::Int=150, basal_frac::Float64=0.3,
                       connectance::Float64=0.10, R95::Int=5, motif_mix::Symbol=:mixed)
    nb = max(1, round(Int, basal_frac*S))
    nc = S - nb
    trophic_role = vcat(fill(:basal, nb), fill(:consumer, nc))

    A = falses(S,S)

    # candidates per predator (kept as in your code)
    function prey_candidates(j)
        if motif_mix == :chains
            return 1:nb
        elseif motif_mix == :mixed
            return 1:(j-1)
        else # :omnivory
            return 1:(j-1)
        end
    end

    # --- target edges are specified over S^2 (like your figures) ---
    total_allowed = sum(j->length(prey_candidates(j)), nb+1:S)
    target_edges  = clamp(round(Int, connectance * S^2), 0, total_allowed)

    # degree draws (same idea as yours)
    λ = (max(target_edges,1) / max(nc,1)) / 3 + 1
    degs = Int[]
    for j in nb+1:S
        kj = clamp(rand(rng, Poisson(λ)), 0, length(prey_candidates(j)))
        push!(degs, kj)
    end
    # stretch toward target
    scale = target_edges / max(1, sum(degs))
    degs = [clamp(round(Int, d*scale), 0, length(prey_candidates(nb+i))) for (i,d) in enumerate(degs)]

    # place links (unchanged logic)
    for (i,j) in enumerate(nb+1:S)
        cands = collect(prey_candidates(j))
        if motif_mix == :chains
            kk = min(length(cands), degs[i])
            sel = kk==0 ? Int[] : sample(rng, cands, kk; replace=false)
        elseif motif_mix == :mixed
            kk = min(length(cands), degs[i])
            sel = kk==0 ? Int[] : sample(rng, cands, kk; replace=false)
        else
            kk = min(length(cands), degs[i])
            if kk>0
                w = [1 + sum(A[:,q]) for q in cands]
                sel = sample(rng, cands, Weights(w), kk; replace=false)
            else
                sel = Int[]
            end
        end
        for q in sel
            A[j,q] = true
        end
    end

    # --- coarse R95 nudge (kept) ---
    outdeg = [sum(A[j,:]) for j in 1:S]
    current_R95 = quantile(outdeg, 0.95)
    tries = 0
    while abs(current_R95 - R95) > 1 && tries < 2000
        j = rand(rng, nb+1:S)
        cands = collect(1:(j-1))
        if current_R95 < R95
            addable = filter(q->!A[j,q], cands)
            if !isempty(addable)
                A[j, rand(rng, addable)] = true
            end
        else
            remable = filter(q->A[j,q], cands)
            if !isempty(remable)
                A[j, rand(rng, remable)] = false
            end
        end
        outdeg = [sum(A[j,:]) for j in 1:S]
        current_R95 = quantile(outdeg, 0.95)
        tries += 1
    end

    # --- FINAL EXACT CONNECTANCE MATCHER (new) ---
    # respects motif-specific sampling; adjusts only by single-edge flips
    function add_one!()
        # pick a predator with available new prey
        js = [j for j in nb+1:S if any(!A[j,q] for q in prey_candidates(j))]
        isempty(js) && return false
        j = rand(rng, js)
        cands = [q for q in prey_candidates(j) if !A[j,q]]
        if isempty(cands); return false; end
        if motif_mix == :omnivory
            w = [1 + sum(A[:,q]) for q in cands]
            q = sample(rng, cands, Weights(w))
        else
            q = rand(rng, cands)
        end
        A[j,q] = true
        return true
    end
    function remove_one!()
        # pick an existing edge uniformly
        edges = [(j,q) for j in nb+1:S for q in prey_candidates(j) if A[j,q]]
        isempty(edges) && return false
        j,q = rand(rng, edges)
        A[j,q] = false
        return true
    end

    cur_edges = count(A)
    goal = target_edges
    guard = 0
    while cur_edges != goal && guard < 100000
        if cur_edges < goal
            add_one!() || break
            cur_edges += 1
        else
            remove_one!() || break
            cur_edges -= 1
        end
        guard += 1
    end
    # ---------------------------------------------------------------

    return Metaweb(S, basal_frac, BitMatrix(A), trophic_role)
end

end # module

# ########## I you run this instead, you'll get the old version where connectance is not enforced
# function build_metaweb(rng::AbstractRNG; S::Int=150, basal_frac::Float64=0.3,
#                        connectance::Float64=0.10, R95::Int=5, motif_mix::Symbol=:mixed)
#     nb = max(1, round(Int, basal_frac*S))
#     nc = S - nb
#     trophic_role = vcat(fill(:basal, nb), fill(:consumer, nc))

#     # target edges ~ C * S*(S-1) but only predator->prey with prey index < predator index (acyclic-ish)
#     A = falses(S,S)

#     # prey pool per consumer depends on motif_mix
#     function prey_candidates(j)
#         # consumers can eat basals and possibly consumers (omnivory)
#         if motif_mix == :chains
#             return 1:nb
#         elseif motif_mix == :mixed
#             return 1:(j-1)  # DAG-ish
#         else # :omnivory
#             return 1:(j-1)
#         end
#     end

#     # stochastic degrees to match connectance and R95 roughly
#     total_possible = sum(j->length(prey_candidates(j)), nb+1:S)
#     target_edges = max(1, round(Int, connectance*total_possible))

#     # degree distribution: truncated Poisson tuned by λ, then clamp by candidates
#     λ = (target_edges / nc) / 3 + 1
#     degs = Int[]
#     for j in nb+1:S
#         kj = clamp(rand(rng, Poisson(λ)), 0, length(prey_candidates(j)))
#         push!(degs, kj)
#     end
#     # stretch to hit target edges approximately
#     scale = target_edges / max(1,sum(degs))
#     degs = [clamp(round(Int, d*scale), 0, length(prey_candidates(nb+i))) for (i,d) in enumerate(degs)]

#     # place links
#     for (i,j) in enumerate(nb+1:S)
#         cands = collect(prey_candidates(j))
#         if motif_mix == :chains
#             # sample without replacement, little sharing
#             kk = min(length(cands), degs[i])
#             sel = kk==0 ? Int[] : sample(rng, cands, kk; replace=false)
#         elseif motif_mix == :mixed
#             kk = min(length(cands), degs[i])
#             sel = kk==0 ? Int[] : sample(rng, cands, kk; replace=false)
#         else
#             # omnivory: bias toward already-popular prey to create sharing
#             kk = min(length(cands), degs[i])
#             if kk>0
#                 w = [1 + sum(A[:,q]) for q in cands]
#                 sel = sample(rng, cands, Weights(w), kk; replace=false)
#             else
#                 sel = Int[]
#             end
#         end
#         for q in sel
#             A[j,q] = true
#         end
#     end

#     # coarse adjust of R95 using simple add/remove on predators with too-low/high out-degree
#     outdeg = [sum(A[j,:]) for j in 1:S]
#     current_R95 = quantile(outdeg, 0.95)
#     tries = 0
#     while abs(current_R95 - R95) > 1 && tries < 2000
#         j = rand(rng, nb+1:S)
#         cands = collect(1:(j-1))
#         if current_R95 < R95
#             # add a prey if possible
#             addable = filter(q->!A[j,q], cands)
#             if !isempty(addable)
#                 A[j, rand(rng, addable)] = true
#             end
#         else
#             # remove a prey if possible
#             remable = filter(q->A[j,q], cands)
#             if !isempty(remable)
#                 A[j, rand(rng, remable)] = false
#             end
#         end
#         outdeg = [sum(A[j,:]) for j in 1:S]
#         current_R95 = quantile(outdeg, 0.95)
#         tries += 1
#     end

#     return Metaweb(S, basal_frac, BitMatrix(A), trophic_role)
# end
# end # module
