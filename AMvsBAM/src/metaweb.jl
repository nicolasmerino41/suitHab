module MetaWeb

using Random, Statistics

# ---------------- Types ----------------
struct Metaweb
    S::Int
    basal_frac::Float64
    A::BitMatrix
    trophic_role::Vector{Symbol}
end

# --------------- helpers ---------------
_sorted_masses(rng::AbstractRNG, S::Int) = sort(rand(rng, S))

function _trophic_roles(S::Int, basal_frac::Float64)
    nb = max(1, round(Int, basal_frac*S))
    v  = fill(:consumer, S)
    v[1:nb] .= :basal
    v
end

_allowed_count(S::Int) = (S*(S-1)) ÷ 2

# --- motif masks (lower triangle only used) ---
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

# --- Bernoulli sampler used when degree_tail=:none ---
function _sample_to_connectance(rng::AbstractRNG, P::Matrix{Float64}, Ctarget::Float64; nb::Int=0)
    S = size(P,1)
    mask = falses(S,S)
    for p in max(nb+1,2):S, q in 1:(p-1)
        mask[p,q] = true
    end
    target_within = Ctarget # C defined over S^2
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

# ------------- heavy-tail toys -------------
@inline function _zipf_trunc(rng::AbstractRNG, α::Float64, kmin::Int, kmax::Int)
    kmin >= kmax && return kmin
    xm = max(kmin, 1)
    u  = rand(rng)
    x  = xm / (1 - u)^(1/α)
    clamp(Int(floor(x)), kmin, kmax)
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
            wi = w[i]; wi <= 0 && continue
            acc += wi
            if acc + eps(acc) >= r
                sel = i; break
            end
        end
        if sel == 0
            sel_idx = findlast(>(0.0), w); sel_idx === nothing && break
            sel = sel_idx::Int
        end
        push!(chosen, sel)
        sumw -= w[sel]; w[sel] = 0.0
    end
end

# --- exact-C matcher (adds/removes until hit target) ---
function _match_connectance_exact!(rng::AbstractRNG, A::BitMatrix, P::Matrix{Float64}, nb::Int, target_edges::Int)
    S = size(A,1); first_cons = max(nb+1,2)
    cur = count(A); guard = 0
    while cur != target_edges && guard < 200000
        if cur < target_edges
            p = rand(rng, first_cons:S)
            cands = Int[]; w = Float64[]
            @inbounds for q in 1:(p-1)
                if !A[p,q] && P[p,q] > 0.0
                    push!(cands,q); push!(w,P[p,q])
                end
            end
            if !isempty(cands)
                s = sum(w); r = rand(rng)*s; acc=0.0; idx=1
                for i in eachindex(cands)
                    acc += w[i]; if acc + eps(acc) >= r; idx=i; break; end
                end
                A[p,cands[idx]] = true; cur += 1
            end
        else
            p = rand(rng, first_cons:S)
            cands = Int[]; w = Float64[]
            @inbounds for q in 1:(p-1)
                if A[p,q]
                    push!(cands,q); push!(w, 1.0/(P[p,q] + 1e-12))
                end
            end
            if !isempty(cands)
                s = sum(w); r = rand(rng)*s; acc=0.0; idx=1
                for i in eachindex(cands)
                    acc += w[i]; if acc + eps(acc) >= r; idx=i; break; end
                end
                A[p,cands[idx]] = false; cur -= 1
            end
        end
        guard += 1
    end
    A
end

# --- degree-preserving overlap shaper (optional) ---
function _edge_swap_preserving_degrees!(A::BitMatrix, P::Matrix{Float64}, p::Int, r::Int, a::Int, b::Int)::Bool
    @inbounds begin
        (p == r || a == b) && return false
        !(A[p,a] & A[r,b]) && return false
        (A[p,b] || A[r,a]) && return false
        (P[p,b] <= 0.0 || P[r,a] <= 0.0) && return false
        A[p,a] = false; A[r,b] = false
        A[p,b] = true;  A[r,a] = true
    end
    true
end

function _shape_overlap_via_swaps!(rng::AbstractRNG, A::BitMatrix, P::Matrix{Float64}, nb::Int; sweeps::Int=10_000)
    S = size(A,1); first_cons = max(nb+1,2)
    for _ in 1:sweeps
        preds = [p for p in first_cons:S if any(@view A[p,1:p-1])]
        length(preds) < 2 && break
        ok=false; p=0; r=0; a=0; b=0
        @inbounds for _try in 1:20
            p = rand(rng, preds); r = rand(rng, preds); p==r && continue
            pl = findall(q->A[p,q]==1, 1:(p-1)); isempty(pl) && continue
            rl = findall(q->A[r,q]==1, 1:(r-1)); isempty(rl) && continue
            a = rand(rng, pl); b = rand(rng, rl)
            if a != b && !A[p,b] && !A[r,a] && P[p,b] > 0.0 && P[r,a] > 0.0
                ok=true; break
            end
        end
        ok || continue
        _edge_swap_preserving_degrees!(A,P,p,r,a,b) || continue
    end
    A
end

# --------- CHAIN convenience ----------
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
    for _ in 1:40
        mid = 0.5*(lo+hi)
        m_within = L_of(mid) / _allowed_count(S)
        (m_within < target_within) ? (lo=mid) : (hi=mid)
    end
    max(1, floor(Int, 0.5*(lo+hi)))
end

function _chains_bandwidth_for_R95(S::Int, targetR95::Int)
    # crude but monotone: wider κ -> larger outdegrees -> larger R95
    lo, hi = 1, S-1
    for _ in 1:30
        mid = (lo+hi) ÷ 2
        # build synthetic outdegrees: row p has min(mid, p-1)
        ks = [min(mid, p-1) for p in 2:S]
        r95 = quantile!(ks, 0.95; sorted=false)
        (r95 < targetR95) ? (lo = mid+1) : (hi = mid-1)
    end
    max(1, lo)
end

function _build_metaweb_chains(rng::AbstractRNG; S::Int, basal_frac::Float64,
                                targetC::Union{Nothing,Float64}, targetR95::Union{Nothing,Int})
    nb = max(1, round(Int, basal_frac*S))
    κ = isnothing(targetC) ? _chains_bandwidth_for_R95(S, targetR95) :
                             _chains_bandwidth_for_target(S, targetC)
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
    A
end

# ---------- main heavy-tail builders ----------
# mode=:C  -> enforce exact connectance, let R95 float
# mode=:R95-> enforce target R95 on out-degree, let C float
function _build_heavytail(rng::AbstractRNG;
        S::Int, nb::Int, P::Matrix{Float64},
        mode::Symbol, targetC::Union{Nothing,Float64}, targetR95::Union{Nothing,Int},
        α_out::Float64=2.2, α_in::Float64=1.6, kmin_out::Int=1)

    first_cons = max(nb+1,2)

    # prey attractiveness (for in-degree skew)
    attract = ones(Float64, S)
    for q in 1:S
        u = rand(rng)
        attract[q] = u^(-1/α_in)
    end

    # 1) pick row targets
    kout_raw = zeros(Int, S)
    for p in first_cons:S
        kout_raw[p] = _zipf_trunc(rng, α_out, min(kmin_out, p-1), p-1)
    end

    kout = zeros(Int, S)
    if mode == :C
        # scale to hit total links = C*S^2 within allowed
        Ltarget = round(Int, targetC * (S^2))
        allowed = sum(p->(p-1), first_cons:S)
        Ltarget = min(Ltarget, allowed)
        function total_links(ρ)
            s = 0
            @inbounds for p in first_cons:S
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
            (total_links(mid) < Ltarget) ? (lo=mid) : (hi=mid)
        end
        ρ = 0.5*(lo+hi)
        @inbounds for p in first_cons:S
            kout[p] = min(p-1, max(0, round(Int, ρ*kout_raw[p])))
        end
    else
        # scale to hit R95 (monotone in ρ)
        function r95_of(ρ)
            ks = Int[]
            @inbounds for p in first_cons:S
                push!(ks, min(p-1, max(0, round(Int, ρ*kout_raw[p]))))
            end
            quantile!(ks, 0.95; sorted=false)
        end
        lo, hi = 0.0, 20.0
        for _ in 1:40
            mid = 0.5*(lo+hi)
            (r95_of(mid) < targetR95) ? (lo=mid) : (hi=mid)
        end
        ρ = 0.5*(lo+hi)
        @inbounds for p in first_cons:S
            kout[p] = min(p-1, max(0, round(Int, ρ*kout_raw[p])))
        end
    end

    # slight jitter
    @inbounds for p in first_cons:S
        kout[p] = clamp(kout[p] + rand(rng, -1:1), 0, p-1)
    end

    # 2) build adjacency honoring P and kout
    A = falses(S,S)
    @inbounds for p in first_cons:S
        k = kout[p]; k == 0 && continue
        w = similar(attract, p-1)
        for q in 1:(p-1)
            w[q] = attract[q] * P[p,q]
        end
        k = min(k, count(>(0.0), w)); k == 0 && continue
        chosen = Int[]; _sample_without_replacement!(rng, chosen, w, k)
        for q in chosen; A[p,q] = true; end
    end

    # 3) optional overlap shaping (degree-preserving)
    _shape_overlap_via_swaps!(rng, A, P, nb; sweeps=5*S)

    # 4) final constraint step depending on mode
    if mode == :C
        allowed = sum(p->(p-1), first_cons:S)
        Ltarget = min(round(Int, targetC*(S^2)), allowed)
        _match_connectance_exact!(rng, A, P, nb, Ltarget)
    else
        # ensure the realized R95 is close (cheap nudges by reassignment)
        outdeg = sum(A; dims=2)[:]
        function r95() quantile(outdeg[first_cons:end], 0.95) end
        sweeps = 5*S
        for _ in 1:sweeps
            cur = r95()
            if cur == targetR95; break; end
            # move one edge from some predator to another to adjust tail
            order = sortperm(outdeg[first_cons:end])
            donor_rel = (cur > targetR95) ? order[end] : order[1]
            recv_rel  = (cur > targetR95) ? order[1]  : order[end]
            donor = donor_rel + first_cons - 1
            recv  = recv_rel  + first_cons - 1
            # donor must have an edge; receiver must have a free allowed slot
            don_q = [q for q in 1:(donor-1) if A[donor,q]]
            isempty(don_q) && continue
            add_q = [q for q in 1:(recv-1) if !A[recv,q] && P[recv,q] > 0.0]
            isempty(add_q) && continue
            q_rm  = rand(rng, don_q)
            q_add = rand(rng, add_q)
            A[donor,q_rm] = false; outdeg[donor]-=1
            A[recv,q_add] = true;  outdeg[recv] +=1
        end
    end

    BitMatrix(A)
end

# --------------- Public API ---------------
"""
build_metaweb(...; control=:auto, connectance=0.10, R95=5, motif_mix=:mixed, degree_tail=:powerlaw, α_out=2.2, α_in=1.6, kmin_out=1, align=0.5)

- `control=:C`   → enforce exact `connectance`; R95 is free (reported ex-post).
- `control=:R95` → enforce target `R95` (95th out-degree); C is free.
- `control=:auto` picks `:R95` if `connectance===nothing` and `R95!==nothing`,
  otherwise `:C`.
"""
function build_metaweb(rng::AbstractRNG; S::Int=175, basal_frac::Float64=0.30,
    control::Symbol=:auto,
    connectance::Union{Nothing,Float64}=0.10,
    R95::Union{Nothing,Int}=5,
    motif_mix::Symbol=:mixed, degree_tail::Symbol=:powerlaw,
    α_out::Float64=2.2, α_in::Float64=1.6, kmin_out::Int=1, align::Float64=0.5)

    roles  = _trophic_roles(S, basal_frac)
    nb     = max(1, round(Int, basal_frac*S))

    # resolve control mode
    mode = control
    if control == :auto
        mode = (connectance === nothing && R95 !== nothing) ? :R95 : :C
    end

    # motif mask (you can keep your alignment-dependent tweaks if you like;
    # here we keep it simple/stable)
    P = if motif_mix == :chains
        _shape_chains(S)
    elseif motif_mix == :omnivory
        _shape_omnivory(S)
    else
        _shape_mixed(S)
    end
    P[1:nb, :] .= 0.0

    if degree_tail == :none
        # Bernoulli case only meaningful for control=:C
        Ctar = (mode == :C) ? (connectance::Float64) : 0.05
        A = _sample_to_connectance(rng, P, Ctar; nb=nb)
        return Metaweb(S, basal_frac, A, roles)
    end

    # heavy-tailed build
    A = if motif_mix == :chains
        _build_metaweb_chains(rng; S=S, basal_frac=basal_frac,
                              targetC = (mode==:C ? connectance::Float64 : nothing),
                              targetR95 = (mode==:R95 ? R95::Int : nothing))
    else
        _build_heavytail(rng; S=S, nb=nb, P=P, mode=mode,
                         targetC= (mode==:C ? connectance::Float64 : nothing),
                         targetR95=(mode==:R95 ? R95::Int : nothing),
                         α_out=α_out, α_in=α_in, kmin_out=kmin_out)
    end

    Metaweb(S, basal_frac, A, roles)
end

# convenience
global_connectance(A::AbstractMatrix{<:Integer}) = count(==(1), A) / (size(A,1)^2)

end # module
