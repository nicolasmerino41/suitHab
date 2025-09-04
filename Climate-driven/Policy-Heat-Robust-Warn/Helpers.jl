nan_to_zero!(v) = (for i in eachindex(v); if isnan(v[i]); v[i]=0.0; end; end; v)
# ΔBSH decomposition (fraction basis; mean over consumers)
# Returns NamedTuple: dB, dA_only, dInt_only, synergy
function decomp_at_mask_fraction(pool, grid; τ::Float64, keepmask::BitVector)
    Z0 = climate_pass(pool, grid; τ=τ);   P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask);        P1 = assemble(Z1, pool)

    cons = .!pool.basal

    A0 = climate_area_per_species(Z0)[cons]
    A1 = climate_area_per_species(Z1)[cons]
    p0 = bsh1_cond_per_species(P0, Z0, pool)[cons]
    p1 = bsh1_cond_per_species(P1, Z1, pool)[cons]
    nan_to_zero!(p0); nan_to_zero!(p1)

    B0 = A0 .* p0
    B1 = A1 .* p1

    dA   = A1 .- A0
    dp   = p1 .- p0
    Abar = 0.5 .* (A0 .+ A1)
    pbar = 0.5 .* (p0 .+ p1)

    dB        = B1 .- B0
    dA_only   = pbar .* dA
    dInt_only = Abar .* dp
    synergy   = dB .- dA_only .- dInt_only  # ~0 (floating error)

    mean(dB), mean(dA_only), mean(dInt_only), mean(synergy)
    return (dB=mean(dB), dA_only=mean(dA_only),
            dInt_only=mean(dInt_only), synergy=mean(synergy))
end

# ---- cell-level scores from a baseline (before any loss) ----
# A_cell  = climate suitability of consumers in each cell (mean over consumers)
# prey_cell = mean, across consumers, of "fraction of that consumer's prey present" in the cell
function cell_scores(Z0::BitMatrix, P0::BitMatrix, pool::SpeciesPool)
    cons = .!pool.basal
    C    = size(Z0,2)
    A_cell = vec(mean(view(Z0, cons, :); dims=1))   # climate score
    # prey support per cell: mean fraction of diet present across consumers
    preyfrac = zeros(sum(cons), C); row = 1
    for s in findall(cons)
        prey = pool.E[s]; diet = max(1, length(prey))
        if diet == 0; row += 1; continue; end
        @inbounds for c in 1:C
            k = 0; @inbounds for q in prey; k += P0[q,c] ? 1 : 0; end
            preyfrac[row,c] = k / diet
        end
        row += 1
    end
    prey_cell = vec(mean(preyfrac; dims=1))
    return A_cell, prey_cell
end

# keep the top fraction of cells by a score (ties broken deterministically)
topk_keepmask(score::AbstractVector{<:Real}, keep_frac::Real) = begin
    C  = length(score)
    k  = clamp(round(Int, keep_frac*C), 1, C)
    ord = sortperm(score; rev=true, by=identity)
    keep = falses(C); keep[ord[1:k]] .= true; keep
end

# climate-decile–matched prey hotspot (prevents “cheating” by climate)
function a_matched_hotspot_mask(A_cell::Vector{Float64}, prey_cell::Vector{Float64};
                                keep_frac::Float64, nbins::Int=10, seed::Int=1)
    rng = MersenneTwister(seed)
    C = length(A_cell); keep = falses(C)
    qs = [quantile(A_cell, i/nbins) for i in 0:nbins]
    for b in 1:nbins
        idx = findall(c -> (A_cell[c] ≥ qs[b]) && (A_cell[c] ≤ qs[b+1]), 1:C)
        isempty(idx) && continue
        kkeep = round(Int, keep_frac * length(idx))
        ord   = sortperm(prey_cell[idx]; rev=true)
        keep[idx[ord[1:kkeep]]] .= true
    end
    return keep
end

# handy collector → always works whether your decomp uses (:dA_only…) or (:clim…)
_get(nt::NamedTuple, alts::Tuple) = (for k in alts; hasproperty(nt,k) && return getfield(nt,k); end; error("Missing $(alts)"))
function collect_decomp_series(vec::Vector{<:NamedTuple})
    clim   = Float64[_get(nt, (:clim,   :dA_only))   for nt in vec]
    inter  = Float64[_get(nt, (:inter,  :dInt_only)) for nt in vec]
    syn    = Float64[_get(nt, (:syn,    :synergy))   for nt in vec]
    tot    = Float64[_get(nt, (:tot,    :dB))        for nt in vec]
    (; clim, inter, syn, tot)
end
function excess_from_decomps(R,S)
    ED_clim = S.clim .- R.clim; ED_int = S.inter .- R.inter; ED_syn = S.syn .- R.syn
    ED_tot  = ED_clim .+ ED_int .+ ED_syn
    RED_tot = ED_tot ./ (abs.(R.tot) .+ 1e-12)
    share   = [abs(ED_tot[i])<1e-12 ? NaN : ED_int[i]/ED_tot[i] for i in eachindex(ED_tot)]
    (; ED_clim, ED_int, ED_syn, ED_tot, RED_tot, share)
end

# Simple degree-preserving rewiring: each consumer keeps diet size; prey replaced uniformly at random
function rewire_metaweb_fraction(pool::SpeciesPool; frac::Float64, seed::Int=42)
    rng = MersenneTwister(seed)
    nb = count(pool.basal); S = pool.S
    newE = [Int[] for _ in 1:S]
    for s in 1:nb; newE[s] = Int[]; end
    for s in nb+1:S
        k   = length(pool.E[s])
        rew = rand(rng) < frac
        if !rew
            newE[s] = copy(pool.E[s])
        else
            newE[s] = rand(rng, 1:nb, k)
        end
    end
    SpeciesPool(pool.S, copy(pool.masses), copy(pool.basal),
                copy(pool.mu), copy(pool.b), newE)
end

# Collect a vector of NamedTuples (one per loss) into arrays of damage components
function collect_damage_series(vec_nt::Vector{<:NamedTuple})
    Dtot = Float64[]; Dcl = Float64[]; Dint = Float64[]; Dsyn = Float64[]
    for nt in vec_nt
        d = damage_from_decomp(nt)
        push!(Dtot, d.D_tot); push!(Dcl, d.D_clim); push!(Dint, d.D_int); push!(Dsyn, d.D_syn)
    end
    return (; D_tot=Dtot, D_clim=Dcl, D_int=Dint, D_syn=Dsyn)
end

# Excess damage vs random + robust interaction share
function excess_from_damage(R, S; eps=1e-8)
    ED_clim = @. S.D_clim - R.D_clim
    ED_int  = @. S.D_int  - R.D_int
    ED_syn  = @. S.D_syn  - R.D_syn
    ED_tot  = @. S.D_tot  - R.D_tot
    denom   = @. abs(ED_clim) + abs(ED_int) + abs(ED_syn)
    share   = [denom[i] < eps ? NaN : abs(ED_int[i]) / denom[i] for i in eachindex(ED_tot)]
    return (; ED_tot, ED_clim, ED_int, ED_syn, share, ED_int_abs=abs.(ED_int))
end
# Convert decomposition → DAMAGE components (up = worse)
damage_from_decomp(nt) = (
    D_tot = -nt.dB,
    D_clim = -nt.dA_only,
    D_int  = -nt.dInt_only,
    D_syn  = -nt.synergy
)

