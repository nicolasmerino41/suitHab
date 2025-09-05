# --- BAM PARAMS --------------------------------------------------------------
Base.@kwdef struct BAMParams
    α::Float64 = 1.0       # abiotic weight (importance of A)
    β::Float64 = 1.0       # biotic weight  (importance of B)
    μ::Float64 = 0.0       # movement weight (importance of M)
    γ::Float64 = 2.0       # steepness of saturating prey support h(x)=1-exp(-γ x)
    τA::Float64 = 0.5      # climate suit. threshold to define "suitable" for sources & stats
    τocc::Float64 = 0.35   # occupancy threshold for binary P after scoring
end

Base.@kwdef struct MovementParams
    mode::Symbol = :none    # :none | :access | :component
    λ::Float64 = 0.10       # distance kernel scale (for :access)
    T::Int = 4              # minimum connected-suitable component size (for :component)
    metric::Symbol = :euclid  # :euclid | :manhattan (neighbors4 uses manhattan)
end

# Continuous suitability 0..1 for every species × cell
function abiotic_matrix(pool::SpeciesPool, grid::Grid)
    S, C = pool.S, grid.C
    A = Matrix{Float64}(undef, S, C)
    @inbounds for s in 1:S, c in 1:C
        A[s,c] = exp(-((grid.climate[c] - pool.mu[s])^2) / (2 * pool.b[s]^2 + 1e-12))
    end
    A
end

# Sources: climate-suitable kept cells for species s
function _source_cells(A::Matrix{Float64}, keep::BitVector, τA::Float64, s::Int)
    findall(i -> keep[i] && A[s,i] ≥ τA, 1:size(A,2))
end

# Euclidean distance (in [0,1] coords)
@inline function _dist2(xy::Matrix{Float64}, i::Int, j::Int)
    dx = xy[1,i]-xy[1,j]; dy = xy[2,i]-xy[2,j]; dx*dx + dy*dy
end

# Accessibility mask: M = 1 - ∏ (1 - exp(-d/λ))
function movement_accessibility(pool::SpeciesPool, grid::Grid, A::Matrix{Float64},
                                keep::BitVector, mp::MovementParams, τA::Float64)
    S, C = pool.S, grid.C
    λ = max(mp.λ, 1e-6)
    M = fill(0.0, S, C)
    @inbounds for s in 1:S
        src = _source_cells(A, keep, τA, s)
        if isempty(src)
            # no sources; species can't reach anything
            continue
        end
        for i in 1:C
            keep[i] || continue
            acc = 1.0
            for j in src
                # small optimization: ignore very far cells
                d2 = _dist2(grid.xy, i, j)
                w  = exp(-sqrt(d2)/λ)
                acc *= (1 - w)
                if acc < 1e-6; break; end
            end
            M[s,i] = 1 - acc
        end
    end
    M
end

# Component requirement: M=1 inside connected components (on kept+climate-suitable) with size ≥T
function movement_component(pool::SpeciesPool, grid::Grid, A::Matrix{Float64},
                            keep::BitVector, mp::MovementParams, τA::Float64)
    S, C = pool.S, grid.C
    M = zeros(Float64, S, C)
    nx, ny = grid.nx, grid.ny
    @inbounds for s in 1:S
        ok = BitVector(undef, C)
        @inbounds for i in 1:C
            ok[i] = keep[i] && (A[s,i] ≥ τA)
        end
        visited = falses(C)
        for i in 1:C
            (ok[i] && !visited[i]) || continue
            # BFS to get one component
            comp = Int[]
            q = [i]; visited[i] = true
            while !isempty(q)
                v = popfirst!(q)
                push!(comp, v)
                for nb in neighbors4(v, nx, ny)
                    if ok[nb] && !visited[nb]
                        visited[nb] = true
                        push!(q, nb)
                    end
                end
            end
            if length(comp) ≥ mp.T
                for v in comp
                    M[s,v] = 1.0
                end
            end
        end
    end
    M
end

# Dispatcher
function movement_matrix(pool::SpeciesPool, grid::Grid, A::Matrix{Float64},
                         keep::BitVector; mp::MovementParams=MovementParams(), τA::Float64=0.5)
    if mp.mode === :none
        return ones(Float64, pool.S, grid.C)
    elseif mp.mode === :access
        return movement_accessibility(pool, grid, A, keep, mp, τA)
    elseif mp.mode === :component
        return movement_component(pool, grid, A, keep, mp, τA)
    else
        error("Unknown movement mode $(mp.mode)")
    end
end
# prey weights: equal or normalized by #prey
function prey_weights(pool::SpeciesPool)
    S = pool.S
    W = Vector{Vector{Float64}}(undef, S)
    for s in 1:S
        prey = pool.E[s]
        if isempty(prey)
            W[s] = Float64[]
        else
            W[s] = fill(1.0/length(prey), length(prey))
        end
    end
    W
end

# Assemble BAM: returns (P::BitMatrix, Bsupport::Matrix{Float64}, Score::Matrix{Float64})
function assemble_BAM(pool::SpeciesPool, grid::Grid, A::Matrix{Float64}, keep::BitVector;
                      bam::BAMParams=BAMParams(), mp::MovementParams=MovementParams())
    S, C = pool.S, grid.C
    # Movement depends on kept cells; compute once
    M = movement_matrix(pool, grid, A, keep; mp=mp, τA=bam.τA)

    # Score placeholders
    P = falses(S, C)                         # binary occupancy
    Bsup = ones(Float64, S, C)               # biotic support 0..1 (1 for basal)
    Score = zeros(Float64, S, C)

    # Topological order = by body mass (your metaweb construction guarantees prey lighter)
    order = sortperm(pool.masses)  # increasing mass: basal first
    W = prey_weights(pool)

    @inbounds for s in order
        prey = pool.E[s]
        for i in 1:C
            keep[i] || continue
            if pool.basal[s] || isempty(prey)
                b = 1.0
            else
                # availability: sum of occupied prey in cell (weighted)
                x = 0.0
                for (k, q) in enumerate(prey)
                    x += W[s][k] * (P[q,i] ? 1.0 : 0.0)
                end
                b = 1.0 - exp(-bam.γ * x)        # saturating, in [0,1)
            end
            Bsup[s,i] = b
            Score[s,i] = (A[s,i]^bam.α) * (b^bam.β) * (M[s,i]^bam.μ)
            P[s,i] = Score[s,i] ≥ bam.τocc
        end
    end

    return P, Bsup, Score, M
end

# Keep the top keep_frac along an axis with noisy threshold => front-like removal
function frontlike_mask(grid::Grid, keep_frac::Float64; axis::Symbol=:x, noise::Float64=0.0, seed::Int=0)
    Random.seed!(seed)
    C = grid.C
    coord = axis === :x ? view(grid.xy, 1, :) : view(grid.xy, 2, :)
    # threshold so that ~keep_frac are kept
    q = quantile(coord, 1 - keep_frac)
    # add small noise to avoid a razor straight front
    thr = q .+ noise .* (rand(C) .- 0.5)
    keep = BitVector(undef, C)
    @inbounds for i in 1:C
        keep[i] = coord[i] ≥ thr[i]
    end
    keep
end

Base.@kwdef struct BAMStats
    A::Vector{Float64}
    B::Vector{Float64}
    M::Vector{Float64}
    F::Vector{Float64}       # A * B^β * M^μ
    BSH::Vector{Float64}     # realized fraction of occupied cells over FULL area
end

function species_stats(pool::SpeciesPool, grid::Grid, A::Matrix{Float64}, keep::BitVector,
                       P::BitMatrix, Bsup::Matrix{Float64}, M::Matrix{Float64}, bam::BAMParams)
    S, C = pool.S, grid.C
    As = zeros(Float64, S); Bs = zeros(Float64, S); Ms = zeros(Float64, S)
    F  = zeros(Float64, S); BSH = zeros(Float64, S)

    @inbounds for s in 1:S
        idx = Int[]
        for i in 1:C
            if keep[i] && A[s,i] ≥ bam.τA
                push!(idx, i)
            end
        end
        if isempty(idx)
            As[s]  = 0.0; Bs[s] = 0.0; Ms[s] = 0.0
            F[s]   = 0.0
            BSH[s] = 0.0
        else
            As[s]  = length(idx) / C
            Bs[s]  = mean(@view Bsup[s, idx])
            Ms[s]  = mean(@view M[s, idx])
            F[s]   = As[s] * (Bs[s]^bam.β) * (Ms[s]^bam.μ)
            BSH[s] = sum(@view P[s, :]) / C
        end
    end
    BAMStats(As, Bs, Ms, F, BSH)
end

# --- 3-factor Shapley main effects (A,B,M); residual = synergy ----------------
Base.@kwdef struct Shapley3{T}
    dF::T
    dA::T
    dB::T
    dM::T
    synergy::T
end

# Evaluate F given scalar A,B,M and exponents
@inline _F(A,B,M,β,μ) = A * (B^β) * (M^μ)

function shapley3(A0,B0,M0, A1,B1,M1; β::Float64=1.0, μ::Float64=0.0)
    # weights for permutations (see derivation in analysis)
    f = _F
    ΔA = (1/3) * ( f(A1,B0,M0,β,μ) - f(A0,B0,M0,β,μ) ) +
         (1/6) * ( f(A1,B1,M0,β,μ) - f(A0,B1,M0,β,μ) ) +
         (1/6) * ( f(A1,B0,M1,β,μ) - f(A0,B0,M1,β,μ) ) +
         (1/3) * ( f(A1,B1,M1,β,μ) - f(A0,B1,M1,β,μ) )

    ΔB = (1/3) * ( f(A0,B1,M0,β,μ) - f(A0,B0,M0,β,μ) ) +
         (1/6) * ( f(A1,B1,M0,β,μ) - f(A1,B0,M0,β,μ) ) +
         (1/6) * ( f(A0,B1,M1,β,μ) - f(A0,B0,M1,β,μ) ) +
         (1/3) * ( f(A1,B1,M1,β,μ) - f(A1,B0,M1,β,μ) )

    ΔM = (1/3) * ( f(A0,B0,M1,β,μ) - f(A0,B0,M0,β,μ) ) +
         (1/6) * ( f(A1,B0,M1,β,μ) - f(A1,B0,M0,β,μ) ) +
         (1/6) * ( f(A0,B1,M1,β,μ) - f(A0,B1,M0,β,μ) ) +
         (1/3) * ( f(A1,B1,M1,β,μ) - f(A1,B1,M0,β,μ) )

    F0 = f(A0,B0,M0,β,μ); F1 = f(A1,B1,M1,β,μ)
    dF = F1 - F0
    Shapley3(dF, ΔA, ΔB, ΔM, dF - (ΔA + ΔB + ΔM))
end

function shapley_per_species(pool::SpeciesPool, grid::Grid, A::Matrix{Float64}, keep0::BitVector, keep1::BitVector;
                             bam::BAMParams=BAMParams(), mp::MovementParams=MovementParams())
    # baseline
    P0, B0, _, M0 = assemble_BAM(pool, grid, A, keep0; bam=bam, mp=mp)
    S0 = species_stats(pool, grid, A, keep0, P0, B0, M0, bam)

    # scenario
    P1, B1, _, M1 = assemble_BAM(pool, grid, A, keep1; bam=bam, mp=mp)
    S1 = species_stats(pool, grid, A, keep1, P1, B1, M1, bam)

    S = pool.S
    out = Vector{Shapley3{Float64}}(undef, S)
    @inbounds for s in 1:S
        out[s] = shapley3(S0.A[s], S0.B[s], S0.M[s],
                          S1.A[s], S1.B[s], S1.M[s];
                          β=bam.β, μ=bam.μ)
    end
    out, S0, S1
end

function run_sweep_and_plot(; nx=40, ny=30, S=80, basal_frac=0.45,
                             τA=0.5, τocc=0.35, α=1.0, β=1.0, μ=0.0, γ=2.0,
                             move_mode=:none, λ=0.12, T=4,
                             hl_kind=:clustered, seeds=6, front_axis=:x, front_noise=0.05,
                             loss_fracs = 0.1:0.1:0.8, seed_grid=42, seeds_pool=1:5, seeds_mask=1:5)

    grid = make_grid(nx, ny; seed=seed_grid)
    pool = build_pool(S; basal_frac=basal_frac, seed=first(seeds_pool))
    A = abiotic_matrix(pool, grid)

    bam = BAMParams(α=α, β=β, μ=μ, γ=γ, τA=τA, τocc=τocc)
    mp  = MovementParams(mode=move_mode, λ=λ, T=T)

    base_keep = trues(grid.C)  # "no loss" baseline

    xs = Float64[]; yΔ = Float64[]; yA = Float64[]; yB = Float64[]; yM = Float64[]; ySyn = Float64[]

    for f in loss_fracs
        keepfrac = 1 - f
        # ensemble over seeds
        vals = Float64[]; vA = Float64[]; vB = Float64[]; vM = Float64[]; vSyn = Float64[]
        for ps in seeds_pool, ms in seeds_mask
            pool = build_pool(S; basal_frac=basal_frac, seed=ps)
            A = abiotic_matrix(pool, grid)

            keep = begin
                if hl_kind === :random
                    random_mask(grid.C, keepfrac; seed=ms)
                elseif hl_kind === :clustered
                    clustered_mask(grid, keepfrac; nseeds=seeds, seed=ms)
                elseif hl_kind === :front
                    frontlike_mask(grid, keepfrac; axis=front_axis, noise=front_noise, seed=ms)
                else
                    error("Unknown hl_kind")
                end
            end

            shap, S0, S1 = shapley_per_species(pool, grid, A, base_keep, keep; bam=bam, mp=mp)
            cons = consumer_mask(pool)
            dF   = mean(getfield.(shap[cons], :dF))
            dA   = mean(getfield.(shap[cons], :dA))
            dB   = mean(getfield.(shap[cons], :dB))
            dM   = mean(getfield.(shap[cons], :dM))
            dSyn = mean(getfield.(shap[cons], :synergy))
            push!(vals, dF); push!(vA, dA); push!(vB, dB); push!(vM, dM); push!(vSyn, dSyn)
        end
        push!(xs, f)
        push!(yΔ,  mean(vals))
        push!(yA,  mean(vA))
        push!(yB,  mean(vB))
        push!(yM,  mean(vM))
        push!(ySyn,mean(vSyn))
    end

    # --- Makie figure (your style) ------------------------------------------
    begin
        fig = Figure(; size=(950, 420))
        ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="ΔF (consumers mean)",
                   title="Total and Shapley components — $(String(hl_kind)) HL")
        lines!(ax1, xs, yΔ, label="Total Δ")
        lines!(ax1, xs, yA, label="Abiotic")
        lines!(ax1, xs, yB, label="Biotic")
        lines!(ax1, xs, yM, label="Movement")
        lines!(ax1, xs, ySyn, label="Synergy (resid.)")
        axislegend(ax1, position=:lt, framevisible=false)
        display(fig)
    end

    return (x=xs, dF=yΔ, dA=yA, dB=yB, dM=yM, dSyn=ySyn)
end

# Example:
res = run_sweep_and_plot(
    ; nx=50, ny=40, S=100,
    α=1.0, β=1.5, μ=0.5, γ=3.0,
    move_mode=:access, λ=0.08,
    hl_kind=:front, front_axis=:x, front_noise=0.03
)

"""
Run a habitat-loss sweep, return summary statistics.
"""

using Base.Threads

function run_sweep(; nx=40, ny=30, S=80, basal_frac=0.45,
                    τA=0.5, τocc=0.35, α=1.0, β=1.0, μ=0.0, γ=2.0,
                    move_mode=:none, λ=0.12, T=4,
                    hl_kind=:clustered, seeds=6, front_axis=:x, front_noise=0.05,
                    loss_fracs=0.1:0.1:0.8, seed_grid=42, seeds_pool=1:5, seeds_mask=1:5)

    grid = make_grid(nx, ny; seed=seed_grid)
    pool = build_pool(S; basal_frac=basal_frac, seed=first(seeds_pool))
    A = abiotic_matrix(pool, grid)

    bam = BAMParams(α=α, β=β, μ=μ, γ=γ, τA=τA, τocc=τocc)
    mp  = MovementParams(mode=move_mode, λ=λ, T=T)

    base_keep = trues(grid.C)

    xs   = Float64[]
    yΔ   = Float64[]
    yA   = Float64[]
    yB   = Float64[]
    yM   = Float64[]
    ySyn = Float64[]

    for f in loss_fracs
        keepfrac = 1 - f

        # Thread-local storage
        ncomb = length(seeds_pool) * length(seeds_mask)
        dFvals   = Vector{Float64}(undef, ncomb)
        dAvals   = similar(dFvals)
        dBvals   = similar(dFvals)
        dMvals   = similar(dFvals)
        dSynvals = similar(dFvals)

        Threads.@threads for idx in 1:ncomb
            ps = seeds_pool[(idx-1) ÷ length(seeds_mask) + 1]
            ms = seeds_mask[(idx-1) % length(seeds_mask) + 1]

            pool = build_pool(S; basal_frac=basal_frac, seed=ps)
            A = abiotic_matrix(pool, grid)

            keep = if hl_kind === :random
                random_mask(grid.C, keepfrac; seed=ms)
            elseif hl_kind === :clustered
                clustered_mask(grid, keepfrac; nseeds=seeds, seed=ms)
            elseif hl_kind === :front
                frontlike_mask(grid, keepfrac; axis=front_axis, noise=front_noise, seed=ms)
            else
                error("Unknown hl_kind")
            end

            shap, S0, S1 = shapley_per_species(pool, grid, A, base_keep, keep; bam=bam, mp=mp)
            cons = consumer_mask(pool)

            dFvals[idx]   = mean(getfield.(shap[cons], :dF))
            dAvals[idx]   = mean(getfield.(shap[cons], :dA))
            dBvals[idx]   = mean(getfield.(shap[cons], :dB))
            dMvals[idx]   = mean(getfield.(shap[cons], :dM))
            dSynvals[idx] = mean(getfield.(shap[cons], :synergy))
        end

        push!(xs, f)
        push!(yΔ,   mean(dFvals))
        push!(yA,   mean(dAvals))
        push!(yB,   mean(dBvals))
        push!(yM,   mean(dMvals))
        push!(ySyn, mean(dSynvals))
    end

    return (x=xs, dF=yΔ, dA=yA, dB=yB, dM=yM, dSyn=ySyn,
            hl_kind=hl_kind)
end

"""
Plot the sweep results returned by `run_sweep`.
"""
function plot_sweep(res)
    fig = Figure(; size=(950, 420))
    ax1 = Axis(fig[1,1],
        xlabel="Area lost (fraction)",
        ylabel="ΔF (consumers mean)",
        title="Total and Shapley components — $(String(res.hl_kind)) HL")

    lines!(ax1, res.x, res.dF,   label="Total Δ")
    lines!(ax1, res.x, res.dA,   label="Abiotic")
    lines!(ax1, res.x, res.dB,   label="Biotic")
    lines!(ax1, res.x, res.dM,   label="Movement")
    lines!(ax1, res.x, res.dSyn, label="Synergy (resid.)")

    axislegend(ax1, position=:lb, framevisible=false)
    display(fig)
end

res = run_sweep(; nx=50, ny=40, S=100,
                  α=1.0, β=1.5, μ=0.5, γ=3.0,
                  move_mode=:access, λ=0.08,
                  hl_kind=:front, front_axis=:x, front_noise=0.03)

plot_sweep(res)
