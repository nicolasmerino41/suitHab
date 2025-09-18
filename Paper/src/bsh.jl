module BSH

using Random, Statistics, Distributions
using ..Grids
using ..Metawebs
using ..HL

export BAMParams, abiotic_matrix, movement_gate, assemble_AM, assemble_BAM,
       relative_loss_curves, per_species_relative_loss, pfail_curve,
       placebo_curves

struct BAMParams
    τA::Float64      # abiotic cell suitability threshold
    τB::Float64      # biotic sufficiency threshold (mean prey presence)
    τocc::Float64    # occupancy threshold (here we use τB directly; keep for future)
    γ::Float64       # smoothness of prey sufficiency (we use mean; γ kept for extension)
    movement::Symbol # :off or :component
    T::Int           # component size for movement gate
end
BAMParams(; τA=0.5, τB=0.35, τocc=0.2, γ=3.0, movement::Symbol=:off, T::Int=8) =
    BAMParams(τA,τB,τocc,γ,movement,T)

struct MovementParams
    mode::Symbol
    T::Int
end
MovementParams(; mode::Symbol=:none, T::Int=8) = MovementParams(mode,T)

# Abiotic suitability A[s,c] (Gaussian niche around species "optimum" equal to climate value)
function abiotic_matrix(pool::SpeciesPool, grid::Grid; niche_width::Float64=0.08, seed::Int=1)
    rng = MersenneTwister(seed)
    # map each species to a niche center uniformly in [0,1]
    mu = rand(rng, pool.S)
    S, C = pool.S, grid.C
    A = Matrix{Float64}(undef, S, C)
    @inbounds for s in 1:S, i in 1:C
        A[s,i] = exp(-((grid.climate[i] - mu[s])^2) / (2*niche_width^2 + 1e-12))
    end
    A
end

"Movement gate per species: among kept & A≥τA cells, allow only components of size≥T."
function movement_gate(grid::Grid, A::Matrix{Float64}, keep::BitVector; τA::Float64, T::Int)
    S, C = size(A)
    nx, ny = grid.nx, grid.ny
    M = falses(S, C)
    @inbounds for s in 1:S
        ok = BitVector(undef, C)
        for i in 1:C; ok[i] = keep[i] & (A[s,i] ≥ τA); end
        visited = falses(C)
        for i in 1:C
            (ok[i] & !visited[i]) || continue
            comp = Int[]; q=[i]; visited[i]=true
            while !isempty(q)
                v = popfirst!(q); push!(comp,v)
                for nb in HL.neighbors4(v, nx, ny)
                    if ok[nb] & !visited[nb]
                        visited[nb]=true; push!(q, nb)
                    end
                end
            end
            if length(comp) ≥ T
                for v in comp; M[s,v]=true; end
            end
        end
    end
    M
end

"AM occupancy: A≥τA & kept; if movement=:component, also require component size ≥ T."
function assemble_AM(pool::SpeciesPool, grid::Grid, A::Matrix{Float64}, keep::BitVector; pars::BAMParams)
    S, C = pool.S, grid.C
    # movement mask per species (or all true if OFF)
    M = pars.movement === :off ? trues(S, C) :
        movement_gate(grid, A, keep; τA=pars.τA, T=pars.T)

    P = falses(S, C)
    @inbounds for s in 1:S, i in 1:C
        P[s,i] = keep[i] & (A[s,i] ≥ pars.τA) & M[s,i]
    end
    bsh = [sum(@view P[s,:]) / C for s in 1:S]
    P, bsh
end

"BAM occupancy: requires A≥τA, kept, movement gate, and prey sufficiency ≥ τB."
# Biotic occupancy with optional prey sufficiency aggregator.
function assemble_BAM(
    pool::Metawebs.SpeciesPool, grid::Grids.Grid, A::Matrix{Float64}, keep::BitVector;
    pars::BAMParams, agg::Symbol=:mean, kreq::Int=1
)
    S, C = pool.S, grid.C
    M = pars.movement === :off ? trues(S,C) : movement_gate(grid, A, keep; τA=pars.τA, T=pars.T)

    P = falses(S, C)
    preyfrac = zeros(Float64, S, C)

    order = sortperm(pool.masses)
    for s in order
        if pool.basal[s]
            @inbounds for i in 1:C
                P[s,i] = keep[i] & (A[s,i] ≥ pars.τA) & M[s,i]
            end
        else
            pr = pool.prey[s]
            if isempty(pr)
                @inbounds for i in 1:C; P[s,i] = false; end
            else
                @inbounds for i in 1:C
                    if keep[i] & (A[s,i] ≥ pars.τA) & M[s,i]
                        if agg === :mean
                            m = mean(@view P[pr, i]); preyfrac[s,i] = m
                            P[s,i] = (m ≥ pars.τB)
                        elseif agg === :min
                            # most stringent: all prey present
                            ok = all(@view P[pr, i]); preyfrac[s,i] = ok ? 1.0 : 0.0
                            P[s,i] = ok
                        elseif agg === :kofn
                            # require at least kreq prey present
                            n_ok = count(@view P[pr, i]); n_tot = length(pr)
                            preyfrac[s,i] = n_tot == 0 ? 0.0 : n_ok / n_tot
                            P[s,i] = (n_ok ≥ kreq)
                        else
                            error("Unknown agg = $agg (use :mean, :min, or :kofn)")
                        end
                    else
                        P[s,i] = false
                    end
                end
            end
        end
    end
    bsh = [sum(@view P[s,:])/C for s in 1:S]
    (; P, bsh, preyfrac)
end


"Relative loss curves for AM and BAM across geometries and loss fractions."
function relative_loss_curves(
    rng::AbstractRNG, pool::SpeciesPool, grid::Grid, pars::BAMParams;
    loss_fracs=0.2:0.2:0.8, seed_A::Int=1, A_fn=abiotic_matrix, 
    agg::Symbol=:mean, kreq::Int=1
)
    # A = abiotic_matrix(pool, grid; seed=seed_A)
    A = A_fn(pool, grid; seed=seed_A)
    geometries = (:random, :clustered, :front)
    nx, ny = grid.nx, grid.ny

    # baseline (no HL)
    keep0 = trues(grid.C)
    _, b0_AM  = assemble_AM(pool, grid, A, keep0; pars=pars)
    # bam = assemble_BAM(pool, grid, A, keep0; pars=pars)
    # b0_BAM = bam.bsh
    bam = assemble_BAM(pool, grid, A, keep0; pars=pars, agg=agg, kreq=kreq)
    b0_BAM = bam.bsh

    res = Dict{Symbol,Any}()
    for g in geometries
        xs = Float64[]; yAM = Float64[]; yBAM = Float64[]
        for f in loss_fracs
            keepfrac = 1 - f
            keep = g===:random    ? HL.random_mask(rng, grid.C, keepfrac) :
                   g===:clustered ? HL.clustered_mask(rng, nx, ny, keepfrac; nseeds=8) :
                                    HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)
            _, bAM = assemble_AM(pool, grid, A, keep; pars=pars)
            # bBAM = assemble_BAM(pool, grid, A, keep; pars=pars).bsh
            bBAM = assemble_BAM(pool, grid, A, keep; pars=pars, agg=agg, kreq=kreq).bsh

            ΔAM  = mean(bAM)  - mean(b0_AM)
            ΔBAM = mean(bBAM) - mean(b0_BAM)
            relAM  = ΔAM / (mean(b0_AM)  + eps())
            relBAM = ΔBAM / (mean(b0_BAM) + eps())

            push!(xs, f); push!(yAM, relAM); push!(yBAM, relBAM)
        end
        res[g] = (; x=xs, relAM=yAM, relBAM=yBAM)
    end
    res
end

"Per-species relative losses at a chosen f* for CDFs."
function per_species_relative_loss(
    rng::AbstractRNG, pool::SpeciesPool, grid::Grid, pars::BAMParams;
    fstar::Float64=0.6, geometry::Symbol=:random, seed_A::Int=1,
    A_fn=abiotic_matrix, agg::Symbol=:mean, kreq::Int=1
)
    # A = abiotic_matrix(pool, grid; seed=seed_A)
    A = A_fn(pool, grid; seed=seed_A)

    keep0 = trues(grid.C)
    _, b0_AM  = assemble_AM(pool, grid, A, keep0; pars=pars)
    # b0_BAM    = assemble_BAM(pool, grid, A, keep0; pars=pars).bsh
    b0_BAM = assemble_BAM(pool, grid, A, keep0; pars=pars, agg=agg, kreq=kreq).bsh

    keepfrac = 1 - fstar
    keep = geometry===:random    ? HL.random_mask(rng, grid.C, keepfrac) :
           geometry===:clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                   HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)
    _, bAM = assemble_AM(pool, grid, A, keep; pars=pars)
    # bBAM   = assemble_BAM(pool, grid, A, keep; pars=pars).bsh
    bBAM   = assemble_BAM(pool, grid, A, keep;  pars=pars, agg=agg, kreq=kreq).bsh

    relAM  = (bAM .- b0_AM) ./ (b0_AM .+ eps())
    relBAM = (bBAM .- b0_BAM) ./ (b0_BAM .+ eps())
    relAM, relBAM
end

"P_fail curve among A-suitable kept cells (BAM only)."
function pfail_curve(;
    rng::AbstractRNG, pool::SpeciesPool, grid::Grid, pars::BAMParams,
    loss_fracs=0.2:0.2:0.8, seed_A::Int=1, A_fn=abiotic_matrix, agg::Symbol=:mean, kreq::Int=1
)
    # A = abiotic_matrix(pool, grid; seed=seed_A)
    A = A_fn(pool, grid; seed=seed_A)
    
    nx, ny = grid.nx, grid.ny
    geometries = (:random, :clustered, :front)
    out = Dict{Symbol,Any}()
    for g in geometries
        xs = Float64[]; ys = Float64[]
        for f in loss_fracs
            keepfrac = 1 - f
            keep = g===:random    ? HL.random_mask(rng, grid.C, keepfrac) :
                   g===:clustered ? HL.clustered_mask(rng, nx, ny, keepfrac; nseeds=8) :
                                    HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)
            
            # bam = assemble_BAM(pool, grid, A, keep; pars=pars)
            bam = assemble_BAM(pool, grid, A, keep; pars=pars, agg=agg, kreq=kreq)
            
            # among A-suitable kept cells (and movement-passing if M is on), compute prey-gate failure
            cons = [s for s in 1:pool.S if !pool.basal[s]]

            # movement mask used for the denominator if M is on
            M = pars.movement === :off ? trues(pool.S, grid.C) :
                movement_gate(grid, A, keep; τA=pars.τA, T=pars.T)

            maskA = falses(pool.S, grid.C)
            for s in cons, i in 1:grid.C
                maskA[s,i] = keep[i] & (A[s,i] ≥ pars.τA) & M[s,i]
            end

            fail = 0; denom=0
            for s in cons, i in 1:grid.C
                if maskA[s,i]
                    denom+=1
                    if bam.preyfrac[s,i] < pars.τB; fail+=1; end
                end
            end
            push!(xs, f); push!(ys, denom==0 ? 0.0 : fail/denom)
        end
        out[g] = (; x=xs, y=ys)
    end
    out
end

"Placebo: rewire metaweb and recompute relative loss curves (BAM only)."
function placebo_curves(
    ; rng::AbstractRNG, pool::SpeciesPool, grid::Grid, pars::BAMParams,
    loss_fracs=0.2:0.2:0.8, seed_A::Int=1, A_fn=abiotic_matrix, agg::Symbol=:mean, kreq::Int=1
)
    placebo = Metawebs.rewire_placebo(rng, pool; nswap=5000)
    
    # A = abiotic_matrix(placebo, grid; seed=seed_A)
    A = A_fn(placebo, grid; seed=seed_A)

    keep0 = trues(grid.C)

    # b0 = assemble_BAM(placebo, grid, A, keep0; pars=pars).bsh
    b0 = assemble_BAM(placebo, grid, A, keep0; pars=pars, agg=agg, kreq=kreq).bsh

    geometries = (:random, :clustered, :front)
    nx, ny = grid.nx, grid.ny
    res = Dict{Symbol,Any}()
    for g in geometries
        xs = Float64[]; y = Float64[]
        for f in loss_fracs
            keepfrac = 1 - f
            keep = g===:random    ? HL.random_mask(rng, grid.C, keepfrac) :
                   g===:clustered ? HL.clustered_mask(rng, nx, ny, keepfrac; nseeds=8) :
                                    HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)
            # b = assemble_BAM(placebo, grid, A, keep; pars=pars).bsh
            b = assemble_BAM(placebo, grid, A, keep;  pars=pars, agg=agg, kreq=kreq).bsh
            Δ = mean(b) - mean(b0)
            push!(xs, f); push!(y, Δ/(mean(b0)+eps()))
        end
        res[g] = (; x=xs, rel=y)
    end
    res
end

# --- NEW: abiotic_matrix_aligned ----------------------------------------------------------
export abiotic_matrix_aligned

"""
Abiotic suitability with optional (i) basal climate bias and (ii) prey–consumer alignment.

- `bias_basal` in [0,1]: 0.0 = uniform; 0.8 puts most basal optima on the low end of climate.
- `align` in [0,1]: 0 = independent consumer optima; 1 = consumer optimum = mean(prey optima).
- Separate niche widths for basal vs consumers allowed.
"""
function abiotic_matrix_aligned(
    pool::Metawebs.SpeciesPool, grid::Grids.Grid;
    niche_basal::Float64=0.10, niche_cons::Float64=0.12,
    bias_basal::Float64=0.0, align::Float64=0.0,
    seed::Int=1
)
    rng = MersenneTwister(seed)
    S, C = pool.S, grid.C

    μ = rand(rng, S)  # start uniform [0,1]

    # bias basal to the low-climate side using a Beta distribution
    if bias_basal > 0
        # mean ≈ bias_basal; moderately concentrated (α+β≈10)
        α = 1.0 + 9.0*bias_basal
        β = 1.0 + 9.0*(1.0 - bias_basal)
        for s in 1:S
            pool.basal[s] && (μ[s] = rand(rng, Beta(α, β)))
        end
    end

    # align consumers' optima toward the mean of their prey optima
    if align > 0
        for s in 1:S
            pool.basal[s] && continue
            pr = pool.prey[s]
            if !isempty(pr)
                μp = mean(@view μ[pr])
                μ[s] = clamp((1-align)*μ[s] + align*μp, 0.0, 1.0)
            end
        end
    end

    A = Matrix{Float64}(undef, S, C)
    @inbounds for s in 1:S, i in 1:C
        σ = pool.basal[s] ? niche_basal : niche_cons
        A[s,i] = exp(-((grid.climate[i] - μ[s])^2) / (2σ^2 + 1e-12))
    end
    A
end

end # module

