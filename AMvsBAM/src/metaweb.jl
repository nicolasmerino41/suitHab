module MetaWeb

using Random, Statistics

# ------------------------
# Types
# ------------------------
struct Metaweb
    S::Int
    basal_frac::Float64
    A::BitMatrix
    trophic_role::Vector{Symbol}
end

# ------------------------
# helpers
# ------------------------
_sorted_masses(rng::AbstractRNG, S::Int) = sort(rand(rng, S))

function _trophic_roles(S::Int, basal_frac::Float64)
    nb = max(1, round(Int, basal_frac*S))
    v  = fill(:consumer, S)
    v[1:nb] .= :basal
    return v
end

_allowed_count(S::Int) = (S*(S-1)) ÷ 2    # strict lower triangle

# --- classic shape matrices (only lower triangle used) ---
function _shape_chains(S::Int; band_frac::Float64=0.08)
    P = zeros(Float64, S, S)
    for p in 2:S
        band = max(1, round(Int, band_frac*p))
        for q in 1:(p-1)
            P[p,q] = (q >= p-band) ? 1.0 : 0.0
        end
    end
    P
end

function _shape_omnivory(S::Int; decay::Float64=0.02)
    P = zeros(Float64, S, S)
    for p in 2:S, q in 1:(p-1)
        d = p - q
        P[p,q] = 1.0 / (1.0 + decay*d)
    end
    P
end

function _shape_mixed(S::Int; band_frac::Float64=0.10, decay::Float64=0.03, w_chain::Float64=0.6)
    Pc = _shape_chains(S; band_frac=band_frac)
    Po = _shape_omnivory(S; decay=decay)
    w  = clamp(w_chain, 0.0, 1.0)
    w*Pc .+ (1-w)*Po
end

# --- simple Bernoulli sampler (previous behavior) ---
function _sample_to_connectance(rng::AbstractRNG, P::Matrix{Float64}, Ctarget::Float64; nb::Int=0)
    S = size(P,1)
    mask = falses(S,S)
    for p in max(nb+1,2):S, q in 1:(p-1)
        mask[p,q] = true
    end
    target_within = Ctarget # interpret C over S^2 consistently with figures
    lo, hi = 0.0, 1e6
    @inbounds for _ in 1:40
        mid = 0.5*(lo+hi)
        m   = mean(@view min.(1.0, mid .* P)[mask])
        (m < target_within) ? (lo = mid) : (hi = mid)
    end
    s = 0.5*(lo+hi)
    Prob = min.(1.0, s .* P)

    A = falses(S,S)
    @inbounds for p in max(nb+1,2):S, q in 1:(p-1)
        (rand(rng) < Prob[p,q]) && (A[p,q] = true)
    end
    BitMatrix(A)
end

# ------------------------
# Heavy-tailed degree tools
# ------------------------
@inline function _zipf_trunc(rng::AbstractRNG, α::Float64, kmin::Int, kmax::Int)
    kmin >= kmax && return kmin
    xm = max(kmin, 1)
    u  = rand(rng)
    x  = xm / (1 - u)^(1/α)
    k  = clamp(Int(floor(x)), kmin, kmax)
    k
end

function _sample_without_replacement!(rng::AbstractRNG, chosen::Vector{Int}, w::Vector{Float64}, k::Int)
    k = min(k, count(>(0.0), w))
    k <= 0 && return
    sumw = sum(w)
    for _ in 1:k
        sumw <= 0 && break
        r = rand(rng) * sumw
        acc = 0.0
        sel = 0
        @inbounds for i in eachindex(w)
            wi = w[i]
            wi <= 0 && continue
            acc += wi
            if acc + eps(acc) >= r
                sel = i
                break
            end
        end
        if sel == 0
            sel_idx = findlast(>(0.0), w)
            sel_idx === nothing && break
            sel = sel_idx::Int
        end
        push!(chosen, sel)
        sumw -= w[sel]
        w[sel] = 0.0
    end
end

# --- optional refinement toward heavier tails (keeps total edges) ---
function _refine_powerlaw!(rng::AbstractRNG, A::BitMatrix, P::Matrix{Float64},
                           nb::Int, attract::Vector{Float64}, kout_target::Vector{Int};
                           sweeps::Int=3*size(A,1))
    S = size(A,1); first_cons = max(nb+1,2)
    for _ in 1:sweeps
        p = rand(rng, first_cons:S)
        kt = kout_target[p]
        koutp = 0
        @inbounds for q in 1:(p-1); koutp += A[p,q]; end
        if koutp < kt
            w = zeros(Float64, p-1)
            @inbounds for q in 1:(p-1)
                w[q] = (!A[p,q]) ? attract[q]*P[p,q] : 0.0
            end
            if sum(w) > 0
                chosen = Int[]; _sample_without_replacement!(rng, chosen, w, 1)
                !isempty(chosen) && (A[p,chosen[1]] = true)
            end
        elseif koutp > kt
            w = zeros(Float64, p-1)
            @inbounds for q in 1:(p-1)
                w[q] = A[p,q] ? 1.0/(attract[q] + 1e-12) : 0.0
            end
            if sum(w) > 0
                chosen = Int[]; _sample_without_replacement!(rng, chosen, w, 1)
                !isempty(chosen) && (A[p,chosen[1]] = false)
            end
        end
    end
    return A
end

# --- NEW: exact-connectance adjuster (adds/removes once) ---
function _match_connectance_exact!(rng::AbstractRNG, A::BitMatrix, P::Matrix{Float64}, nb::Int, target_edges::Int)
    S = size(A,1); first_cons = max(nb+1,2)
    cur = count(A)
    guard = 0
    while cur != target_edges && guard < 200000
        if cur < target_edges
            # add one in allowed empty cell; bias by P
            p = rand(rng, first_cons:S)
            cands = Int[]
            w     = Float64[]
            @inbounds for q in 1:(p-1)
                if !A[p,q] && P[p,q] > 0.0
                    push!(cands, q); push!(w, P[p,q])
                end
            end
            if !isempty(cands)
                # simple weighted choice
                s = sum(w); r = rand(rng)*s; acc=0.0; idx=1
                for i in eachindex(cands)
                    acc += w[i]
                    if acc + eps(acc) >= r
                        idx = i; break
                    end
                end
                A[p, cands[idx]] = true
                cur += 1
            end
        else
            # remove one existing edge; light bias against high P
            p = rand(rng, first_cons:S)
            cands = Int[]
            w     = Float64[]
            @inbounds for q in 1:(p-1)
                if A[p,q]
                    push!(cands, q); push!(w, 1.0/(P[p,q] + 1e-12))
                end
            end
            if !isempty(cands)
                s = sum(w); r = rand(rng)*s; acc=0.0; idx=1
                for i in eachindex(cands)
                    acc += w[i]
                    if acc + eps(acc) >= r
                        idx = i; break
                    end
                end
                A[p, cands[idx]] = false
                cur -= 1
            end
        end
        guard += 1
    end
    return A
end

# --- NEW: R95 tuner via edge reassignments (keeps total edges fixed) ---
function _redistribute_edges_for_R95!(rng::AbstractRNG, A::BitMatrix, P::Matrix{Float64}, nb::Int, R95_target::Int; sweeps::Int=5000)
    S = size(A,1); first_cons = max(nb+1,2)
    outdeg = sum(A; dims=2)[:]
    function r95()
        quantile(outdeg[first_cons:end], 0.95)
    end
    for _ in 1:sweeps
        cur = r95()
        # choose donor & receiver depending on direction needed
        if cur < R95_target
            # increase tail: move one edge from a lower-degree predator to a higher-degree one
            order = sortperm(outdeg[first_cons:end])
            donor_rel = order[1]             # lowest
            recv_rel  = order[end]           # highest
        else
            # decrease tail: move from highest to lowest
            order = sortperm(outdeg[first_cons:end])
            donor_rel = order[end]
            recv_rel  = order[1]
        end
        donor = donor_rel + first_cons - 1
        recv  = recv_rel  + first_cons - 1
        donor == recv && continue

        # pick an existing edge in donor
        don_qs = Int[]
        @inbounds for q in 1:(donor-1)
            A[donor,q] && push!(don_qs, q)
        end
        isempty(don_qs) && continue
        q_rm = rand(rng, don_qs)

        # pick a candidate new prey for receiver (allowed & not already linked)
        recv_cands = Int[]
        @inbounds for q in 1:(recv-1)
            (!A[recv,q] && P[recv,q] > 0.0) && push!(recv_cands, q)
        end
        isempty(recv_cands) && continue
        q_add = rand(rng, recv_cands)

        # try the reassignment
        A[donor,q_rm] = false
        A[recv,q_add] = true
        outdeg[donor] -= 1
        outdeg[recv]  += 1

        newr = r95()
        improve = abs(newr - R95_target) < abs(cur - R95_target)
        if !improve && rand(rng) > 0.05  # small exploration chance
            # revert
            A[recv,q_add]  = false
            A[donor,q_rm]  = true
            outdeg[donor] += 1
            outdeg[recv]  -= 1
        end
    end
    return A
end

# --- main heavy-tail builder (returns BitMatrix) ---
function _build_heavytail(rng::AbstractRNG;
        S::Int, nb::Int, Ctarget::Float64, P::Matrix{Float64},
        α_out::Float64=2.2, α_in::Float64=1.6, kmin_out::Int=1, R95::Int=5)

    # target links over S×S but cap by allowed cells
    Ltarget = round(Int, Ctarget * (S^2))
    allowed = sum(p->(p-1), max(nb+1,2):S)
    Ltarget = min(Ltarget, allowed)

    # 1) draw heavy-tailed out-degree targets
    kout_raw = zeros(Int, S)
    for p in max(nb+1,2):S
        kout_raw[p] = _zipf_trunc(rng, α_out, min(kmin_out, p-1), p-1)
    end
    function total_links(ρ)
        s = 0
        @inbounds for p in max(nb+1,2):S
            s += min(p-1, max(0, round(Int, ρ*kout_raw[p])))
        end
        s
    end
    if sum(kout_raw) == 0
        return falses(S,S)
    end
    lo, hi = 0.0, 10.0
    for _ in 1:40
        mid = 0.5*(lo+hi)
        (total_links(mid) < Ltarget) ? (lo = mid) : (hi = mid)
    end
    ρ = 0.5*(lo+hi)
    kout = zeros(Int, S)
    @inbounds for p in max(nb+1,2):S
        kout[p] = min(p-1, max(0, round(Int, ρ*kout_raw[p])))
    end
    @inbounds for p in max(nb+1,2):S
        kout[p] = clamp(kout[p] + rand(rng, -1:1), 0, p-1)
    end
    diff = Ltarget - sum(kout)
    while diff != 0
        p = rand(rng, max(nb+1,2):S)
        if diff > 0 && kout[p] < (p-1)
            kout[p] += 1; diff -= 1
        elseif diff < 0 && kout[p] > 0
            kout[p] -= 1; diff += 1
        end
    end

    # 2) heavy-tailed prey "attractiveness" (for in-degree skew)
    attract = ones(Float64, S)
    for q in 1:S
        u = rand(rng)
        attract[q] = u^(-1/α_in)
    end

    # 3) initial placement honoring P and kout
    A = falses(S,S)
    @inbounds for p in max(nb+1,2):S
        k = kout[p]; k == 0 && continue
        w = similar(attract, p-1)
        for q in 1:(p-1)
            w[q] = attract[q] * P[p,q]
        end
        k = min(k, count(>(0.0), w)); k == 0 && continue
        chosen = Int[]; _sample_without_replacement!(rng, chosen, w, k)
        for q in chosen; A[p,q] = true; end
    end

    # 4) OPTIONAL: nudge within rows toward target kout while favoring popular prey
    _refine_powerlaw!(rng, A, P, nb, attract, kout; sweeps=2*S)

    # 5) EXACT CONNECTANCE (already in your code)
    _match_connectance_exact!(rng, A, P, nb, Ltarget)

    # 6) R95 tuning by reassignments (already in your code)
    _redistribute_edges_for_R95!(rng, A, P, nb, R95; sweeps=5*S)

    # 7) NEW: sculpt predator diet overlap without changing kout/kin/C
    _shape_overlap_via_swaps!(rng, A, P, nb; sweeps=25*S, temp0=0.08)

    return BitMatrix(A)
end

# --- chains variant that widens band to hit target ---
function _chains_bandwidth_for_target(S::Int, Ctarget::Float64)
    f_allowed = (_allowed_count(S)) / (S^2)
    target_within = min(1.0, Ctarget / max(f_allowed, 1e-12))

    L_of = function (κ)
        κi = floor(Int, κ)
        t = min(S-1, κi+1)
        L_small = (t <= 1) ? 0 : sum(p->(p-1), 2:t)
        L_large = (t >= S) ? 0 : (S - t) * κi
        L_small + L_large
    end

    lo, hi = 1.0, float(S-1)
    allowed = _allowed_count(S)
    for _ in 1:40
        mid = 0.5*(lo+hi)
        m_within = L_of(mid)/allowed
        (m_within < target_within) ? (lo = mid) : (hi = mid)
    end
    κ = 0.5*(lo+hi)
    max(1, floor(Int, κ))
end

function _build_metaweb_chains(rng::AbstractRNG; S::Int, basal_frac::Float64, Ctarget::Float64)
    nb = max(1, round(Int, basal_frac*S))
    κ = _chains_bandwidth_for_target(S, Ctarget)
    A = falses(S,S)
    for p in (nb+1):S
        k = min(κ, p-1)
        if k > 1
            δ = rand(rng, (-1):1)
            k = clamp(k + δ, 1, p-1)
        end
        qlo = p - k
        @inbounds for q in qlo:(p-1)
            A[p,q] = true
        end
    end
    return A
end

# ---- fast overlap score on a small predator sample (keeps it cheap)
function _overlap_objective(A::BitMatrix; first_cons::Int, sample_n::Int=30)
    S = size(A,1)
    preds = first_cons:S
    if length(preds) <= 1; return 0.0; end
    # sample a subset of predators for speed
    idx = (length(preds) <= sample_n) ? collect(preds) :
                                        rand(preds, sample_n)
    # collect their prey sets as bitvectors
    vs = [view(A, i, :) for i in idx]
    # variance of pairwise Jaccard overlaps (higher var => mixed guilds/specialists)
    acc = Float64[]; push!(acc, )  # dummy to keep type stable
    empty!(acc)
    for ii in 1:length(vs)-1, jj in ii+1:length(vs)
        a = vs[ii]; b = vs[jj]
        inter = 0; uni = 0
        @inbounds for q in eachindex(a)
            ai = a[q]; bi = b[q]
            inter += (ai & bi)
            uni   += (ai | bi)
        end
        j = (uni == 0) ? 0.0 : inter / uni
        push!(acc, j)
    end
    isempty(acc) && return 0.0
    return var(acc)      # use variance as our objective
end

# ---- masked double-edge swap that preserves kout/kin and C
# Swaps (p,a),(r,b) -> (p,b),(r,a) if the two new edges are allowed by mask P
function _edge_swap_preserving_degrees!(A::BitMatrix, P::Matrix{Float64},
                                        p::Int, r::Int, a::Int, b::Int)::Bool
    @inbounds begin
        if p == r || a == b; return false; end
        # existing edges must exist
        if !(A[p,a] & A[r,b]); return false; end
        # proposed edges must not already exist and must be allowed by mask
        if A[p,b] || A[r,a]; return false; end
        if P[p,b] <= 0.0 || P[r,a] <= 0.0; return false; end
        # perform swap
        A[p,a] = false; A[r,b] = false
        A[p,b] = true;  A[r,a] = true
    end
    return true
end

# ---- anneal diet-overlap while preserving degrees/connectance
# ---- anneal diet-overlap while preserving degrees/connectance
function _shape_overlap_via_swaps!(rng::AbstractRNG, A::BitMatrix, P::Matrix{Float64},
                                   nb::Int; sweeps::Int=20_000, temp0::Float64=0.05)
    S = size(A,1); first_cons = max(nb+1,2)

    # current objective
    score = _overlap_objective(A; first_cons=first_cons)

    for t in 1:sweeps
        T = temp0 / (1 + 0.0005*t)

        # pick predators that have at least one existing link
        preds = [p for p in first_cons:S if any(@view A[p,1:p-1])]
        length(preds) < 2 && break

        ok = false
        p = 0; r = 0; a = 0; b = 0

        # try a few proposals to find a feasible double-swap
        @inbounds for _try in 1:20
            p = rand(rng, preds)
            r = rand(rng, preds)
            p == r && continue

            pl = findall(q->A[p,q]==1, 1:(p-1))
            rl = findall(q->A[r,q]==1, 1:(r-1))
            if isempty(pl) || isempty(rl)
                continue
            end

            a = rand(rng, pl)
            b = rand(rng, rl)
            if a != b && !A[p,b] && !A[r,a] && P[p,b] > 0.0 && P[r,a] > 0.0
                ok = true
                break
            end
        end
        ok || continue

        # perform swap
        A[p,a] = false;  A[r,b] = false
        A[p,b] = true;   A[r,a] = true

        newscore = _overlap_objective(A; first_cons=first_cons)
        Δ = newscore - score
        if Δ >= 0 || rand(rng) < exp(Δ / max(T, 1e-6))
            score = newscore   # accept
        else
            # revert
            A[p,b] = false;  A[r,a] = false
            A[p,a] = true;   A[r,b] = true
        end
    end
    return A
end

# ------------------------
# Public API
# ------------------------
function build_metaweb(
    rng; S::Int=175, basal_frac::Float64=0.30,
    connectance::Float64=0.10, R95::Int=5, motif_mix::Symbol=:mixed,
    degree_tail::Symbol=:powerlaw, α_out::Float64=2.2, α_in::Float64=1.6,
    kmin_out::Int=1, align::Float64=0.5
)

    masses = _sorted_masses(rng, S)
    roles  = _trophic_roles(S, basal_frac)
    nb     = max(1, round(Int, basal_frac*S))

    if motif_mix == :chains
        A = _build_metaweb_chains(rng; S=S, basal_frac=basal_frac, Ctarget=connectance)
        return Metaweb(S, basal_frac, BitMatrix(A), roles)
    end

    # map R95 to mild shape tweaks (original knobs)
    r     = clamp(R95, 1, 10)
    band0 = 0.06 + 0.006*(r-5)   # baseline band for chains component
    dec0  = 0.03 - 0.002*(r-5)   # baseline decay for omnivory component

    # --- NEW: alignment-dependent motif shaping (ecological compensation) ---
    # phi controls strength of compensation; 0.0 = off, 1.0 = strong
    phi = 0.8
    cent = align - 0.5              # [-0.5, +0.5]
    # when align is LOW (cent<0):
    #   - widen band (more vertical reach),
    #   - lower decay (more long-range/omnivory),
    #   - mix more toward omnivory
    band   = band0 * (1.0 + phi * (-cent))              # wider if align low
    decay  = max(1e-3, dec0  * (1.0 - 0.9*phi * (-cent))) # shallower if align low
    w_chain_eff = clamp(0.6 - 0.5*phi * (-cent), 0.1, 0.9) # shift weight toward omnivory at low align

    P = if motif_mix == :omnivory
        _shape_omnivory(S; decay=decay)
    elseif motif_mix == :mixed
        _shape_mixed(S; band_frac=band, decay=decay, w_chain=w_chain_eff)
    else
        # :chains – handled later by _build_metaweb_chains
        _shape_chains(S; band_frac=band)  # only used for degree_tail=:none branch
    end

    # zero-out basal rows
    P[1:nb, :] .= 0.0
    
    if degree_tail == :none
        A = _sample_to_connectance(rng, P, connectance; nb=nb)
        return Metaweb(S, basal_frac, A, roles)
    end
    # inside _build_heavytail call, compute α_out_eff from align:
    # α_out_eff = clamp( (α_out + 1.2*(align - 0.5)), 1.2, 3.5)
    α_out_eff = clamp(α_out + 2.0*(align - 0.5), 1.2, 4.0)  # ±1 shift across 0..1
    # low align (~0) → α_out_eff ≈ α_out - 0.6 (more generalists)
    # high align (~1) → α_out_eff ≈ α_out + 0.6 (more specialists)

    A = _build_heavytail(rng; S=S, nb=nb, Ctarget=connectance, P=P,
                     α_out=α_out_eff, α_in=α_in, kmin_out=kmin_out, R95=R95)

    return Metaweb(S, basal_frac, A, roles)
end

# convenience
global_connectance(A::AbstractMatrix{<:Integer}) = count(==(1), A) / (size(A,1)^2)

end # module
