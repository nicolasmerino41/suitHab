# ╭─────────────────────────────
# │ 1) Grid & helpers
# ╰─────────────────────────────
struct Grid
    nx::Int; ny::Int; C::Int
    xy::Matrix{Float64}      # 2×C
    climate::Vector{Float64} # C
end

function make_grid(nx::Int, ny::Int; seed::Int=42, texture::Float64=0.10)
    rng = MersenneTwister(seed)
    xs = range(0, 1; length=nx)
    ys = range(0, 1; length=ny)
    C = nx*ny
    xy   = Matrix{Float64}(undef, 2, C)
    clim = Vector{Float64}(undef, C)
    k=1
    for j in 1:ny, i in 1:nx
        x=xs[i]; y=ys[j]
        xy[:,k] .= (x,y)
        base   = 0.7x + 0.3y
        wiggle = texture * sin(6π*x) * sin(6π*y)
        clim[k] = base + wiggle
        k+=1
    end
    clim .-= minimum(clim); clim ./= (maximum(clim)+eps())
    Grid(nx,ny,C,xy,clim)
end

neighbors4(ix::Int, nx::Int, ny::Int) = begin
    i = ((ix - 1) % nx) + 1
    j = ((ix - 1) ÷ nx) + 1
    out = Int[]
    if i>1;  push!(out, ix-1);   end
    if i<nx; push!(out, ix+1);   end
    if j>1;  push!(out, ix-nx);  end
    if j<ny; push!(out, ix+nx);  end
    out
end

# thread-local RNG
@inline thread_rng(seed::Int, tags...) = MersenneTwister(hash((seed, tags..., Threads.threadid())))

# ╭─────────────────────────────
# │ 2) Species pool & metaweb
# ╰─────────────────────────────
struct SpeciesPool
    S::Int
    masses::Vector{Float64}
    basal::BitVector
    mu::Vector{Float64}  # niche centers
    b::Vector{Float64}   # niche breadths
    E::Vector{Vector{Int}}  # prey lists (empty for basal)
end

function build_pool(; S::Int, rng::AbstractRNG, basal_frac::Float64,
                     # metaweb knobs
                     R0_mean::Float64=12.0, R0_sd::Float64=0.5, sigma::Float64=0.5,
                     density::Float64=0.30, pmax::Float64=0.90,
                     # niche knobs
                     niche_mode::Symbol=:uniform,
                     mu_basal_centers::Tuple{Float64,Float64}=(0.25,0.75),
                     mu_basal_sd::Float64=0.05,
                     b0_basal::Float64=0.12, bspread_basal::Float64=0.05,
                     b0_cons::Float64=0.12,  bspread_cons::Float64=0.05)

    logm   = collect(range(log(1e-2), log(10.0); length=S))
    shuffle!(rng, logm)
    masses = exp.(logm)

    order = sortperm(masses)
    nB    = clamp(round(Int, basal_frac*S), 0, S)
    basal = falses(S); basal[order[1:nB]] .= true

    mu = similar(masses); b = similar(masses)
    cons_ids = findall(!, basal)
    mu[cons_ids] .= rand(rng, length(cons_ids))
    b[cons_ids]  .= b0_cons .+ bspread_cons .* rand(rng, length(cons_ids))

    bas_ids = findall(basal)
    if niche_mode === :bimodal
        nb = length(bas_ids); nb1 = nb ÷ 2; nb2 = nb - nb1
        c1, c2 = mu_basal_centers
        mu[bas_ids[1:nb1]]     .= clamp.(c1 .+ mu_basal_sd .* randn(rng, nb1), 0, 1)
        mu[bas_ids[nb1+1:end]] .= clamp.(c2 .+ mu_basal_sd .* randn(rng, nb2), 0, 1)
    else
        mu[bas_ids] .= rand(rng, length(bas_ids))
    end
    b[bas_ids] .= b0_basal .+ bspread_basal .* rand(rng, length(bas_ids))

    R0 = exp.(log(R0_mean) .+ R0_sd .* randn(rng, S))
    E = [Int[] for _ in 1:S]
    for (ii,s) in pairs(order)
        basal[s] && continue
        for jj in 1:ii-1
            q = order[jj]
            r = masses[s]/masses[q]
            z = (log(r) - log(R0[s])) / sigma
            p = pmax * exp(-0.5*z^2) * density
            if rand(rng) < p; push!(E[s], q); end
        end
        if isempty(E[s]) && ii>1
            cand = order[1:ii-1]
            target = log(masses[s]) - log(R0[s])
            qstar = cand[argmin(abs.(log.(masses[cand]) .- target))]
            push!(E[s], qstar)
        end
    end
    SpeciesPool(S, masses, basal, mu, b, E)
end

# Abiotic suitability A[s,c] (Gaussian)
function abiotic_matrix(pool::SpeciesPool, grid::Grid)
    S, C = pool.S, grid.C
    A = Matrix{Float64}(undef, S, C)
    @inbounds for s in 1:S, c in 1:C
        A[s,c] = exp(-((grid.climate[c] - pool.mu[s])^2) / (2*pool.b[s]^2 + 1e-12))
    end
    A
end

# ╭─────────────────────────────
# │ 3) Movement (component only) + masks
# ╰─────────────────────────────
struct BAMParams
    α::Float64; β::Float64; μ::Float64; γ::Float64; τA::Float64; τocc::Float64
end
BAMParams(; α=1.0, β=1.0, μ=0.0, γ=2.0, τA=0.5, τocc=0.35) = BAMParams(α,β,μ,γ,τA,τocc)

struct MovementParams
    mode::Symbol
    T::Int
end
MovementParams(; mode::Symbol=:none, T::Int=8) = MovementParams(mode,T)

# Component requirement: within kept+climate-suitable cells, M=1 in components with size≥T
function movement_component(pool::SpeciesPool, grid::Grid, A::Matrix{Float64},
                            keep::BitVector, τA::Float64, T::Int)
    S, C = pool.S, grid.C
    M = zeros(Float64, S, C)
    nx, ny = grid.nx, grid.ny
    @inbounds for s in 1:S
        ok = BitVector(undef, C)
        for i in 1:C
            ok[i] = keep[i] && (A[s,i] ≥ τA)
        end
        visited = falses(C)
        for i in 1:C
            (ok[i] && !visited[i]) || continue
            comp = Int[]; q = [i]; visited[i]=true
            while !isempty(q)
                v = popfirst!(q)
                push!(comp, v)
                for nb in neighbors4(v, nx, ny)
                    if ok[nb] && !visited[nb]
                        visited[nb]=true; push!(q, nb)
                    end
                end
            end
            if length(comp) ≥ T
                for v in comp; M[s,v] = 1.0; end
            end
        end
    end
    M
end

movement_matrix(pool::SpeciesPool, grid::Grid, A::Matrix{Float64}, keep::BitVector;
                mp::MovementParams, τA::Float64) =
    mp.mode === :none ? ones(Float64, pool.S, grid.C) :
    movement_component(pool, grid, A, keep, τA, mp.T)

# Masks (thread-safe)
function random_mask(rng::AbstractRNG, C::Int, keep_frac::Float64)
    nkeep = clamp(round(Int, keep_frac*C), 0, C)
    keep = falses(C); keep[randperm(rng, C)[1:nkeep]] .= true; keep
end
function clustered_mask(rng::AbstractRNG, grid::Grid, keep_frac::Float64; nseeds::Int=6)
    C = grid.C
    target_remove = C - clamp(round(Int, keep_frac*C), 0, C)
    removed = falses(C); q = Int[]
    seeds = randperm(rng, C)[1:min(nseeds,C)]; append!(q, seeds)
    removed_ct=0; ptr=1
    while removed_ct < target_remove
        if ptr>length(q)
            while true
                cand = rand(rng, 1:C)
                if !removed[cand]; push!(q, cand); break; end
            end
        end
        v = q[ptr]; ptr+=1
        if removed[v]; continue; end
        removed[v]=true; removed_ct+=1
        for nb in neighbors4(v, grid.nx, grid.ny)
            !removed[nb] && push!(q, nb)
        end
    end
    .!removed
end
function frontlike_mask(rng::AbstractRNG, grid::Grid, keep_frac::Float64; axis::Symbol=:x, noise::Float64=0.0)
    C = grid.C
    coord = axis === :x ? view(grid.xy, 1, :) : view(grid.xy, 2, :)
    q = quantile(coord, 1 - keep_frac)
    keep = BitVector(undef, C)
    if noise ≤ 0
        @inbounds for i in 1:C; keep[i] = coord[i] ≥ q; end
    else
        ϵ = rand(rng, C) .- 0.5
        @inbounds for i in 1:C; keep[i] = coord[i] ≥ (q + noise*ϵ[i]); end
    end
    keep
end

# ╭─────────────────────────────
# │ 4) BAM assembly (no time)
# ╰─────────────────────────────
consumer_mask(pool::SpeciesPool) = .!pool.basal

function prey_weights(pool::SpeciesPool)
    S = pool.S
    W = Vector{Vector{Float64}}(undef, S)
    for s in 1:S
        prey = pool.E[s]
        W[s] = isempty(prey) ? Float64[] : fill(1.0/length(prey), length(prey))
    end
    W
end

function assemble_BAM(pool::SpeciesPool, grid::Grid, A::Matrix{Float64}, keep::BitVector;
                      bam::BAMParams=BAMParams(), mp::MovementParams=MovementParams())
    S, C = pool.S, grid.C
    M = movement_matrix(pool, grid, A, keep; mp=mp, τA=bam.τA)
    P = falses(S, C)
    Bsup = ones(Float64, S, C)
    Score = zeros(Float64, S, C)
    order = sortperm(pool.masses)  # prey -> predators
    W = prey_weights(pool)

    @inbounds for s in order
        prey = pool.E[s]
        for i in 1:C
            keep[i] || continue
            b = (pool.basal[s] || isempty(prey)) ? 1.0 :
                (begin
                    x = 0.0
                    for (k,q) in enumerate(prey); x += W[s][k] * (P[q,i] ? 1.0 : 0.0); end
                    1.0 - exp(-bam.γ * x)
                 end)
            Bsup[s,i] = b
            Score[s,i] = (A[s,i]^bam.α) * (b^bam.β) * (M[s,i]^bam.μ)
            P[s,i] = Score[s,i] ≥ bam.τocc
        end
    end
    return P, Bsup, Score, M
end

# Species-level summaries on a mask
struct BAMStats
    A::Vector{Float64}; B::Vector{Float64}; M::Vector{Float64}; F::Vector{Float64}; BSH::Vector{Float64}
end

function species_stats(pool::SpeciesPool, grid::Grid, A::Matrix{Float64}, keep::BitVector,
                       P::BitMatrix, Bsup::Matrix{Float64}, M::Matrix{Float64}, bam::BAMParams)
    S, C = pool.S, grid.C
    As = zeros(Float64, S); Bs = zeros(Float64, S); Ms=zeros(Float64, S)
    F  = zeros(Float64, S); BSH=zeros(Float64, S)
    @inbounds for s in 1:S
        idx = Int[]
        for i in 1:C
            if keep[i] && A[s,i] ≥ bam.τA; push!(idx, i); end
        end
        if isempty(idx)
            As[s]=0; Bs[s]=0; Ms[s]=0; F[s]=0; BSH[s]=0
        else
            As[s]  = length(idx) / C
            Bs[s]  = mean(@view Bsup[s, idx])
            Ms[s]  = mean(@view M[s, idx])
            F[s]   = As[s] * (Bs[s]^bam.β) * (Ms[s]^bam.μ)
            BSH[s] = sum(@view P[s, :]) / C
        end
    end
    BAMStats(As,Bs,Ms,F,BSH)
end

# Shapley 3-factor on F(A,B,M)
struct Shapley3{T}; dF::T; dA::T; dB::T; dM::T; synergy::T; end
@inline _F(A,B,M,β,μ) = A * (B^β) * (M^μ)
function shapley3(A0,B0,M0, A1,B1,M1; β::Float64=1.0, μ::Float64=0.0)
    f = _F
    ΔA = (1/3)*(f(A1,B0,M0,β,μ)-f(A0,B0,M0,β,μ)) +
         (1/6)*(f(A1,B1,M0,β,μ)-f(A0,B1,M0,β,μ)) +
         (1/6)*(f(A1,B0,M1,β,μ)-f(A0,B0,M1,β,μ)) +
         (1/3)*(f(A1,B1,M1,β,μ)-f(A0,B1,M1,β,μ))
    ΔB = (1/3)*(f(A0,B1,M0,β,μ)-f(A0,B0,M0,β,μ)) +
         (1/6)*(f(A1,B1,M0,β,μ)-f(A1,B0,M0,β,μ)) +
         (1/6)*(f(A0,B1,M1,β,μ)-f(A0,B0,M1,β,μ)) +
         (1/3)*(f(A1,B1,M1,β,μ)-f(A1,B0,M1,β,μ))
    ΔM = (1/3)*(f(A0,B0,M1,β,μ)-f(A0,B0,M0,β,μ)) +
         (1/6)*(f(A1,B0,M1,β,μ)-f(A1,B0,M0,β,μ)) +
         (1/6)*(f(A0,B1,M1,β,μ)-f(A0,B1,M0,β,μ)) +
         (1/3)*(f(A1,B1,M1,β,μ)-f(A1,B1,M0,β,μ))
    dF = f(A1,B1,M1,β,μ) - f(A0,B0,M0,β,μ)
    Shapley3(dF, ΔA, ΔB, ΔM, dF - (ΔA+ΔB+ΔM))
end

function shapley_per_species(pool::SpeciesPool, grid::Grid, A::Matrix{Float64}, keep0::BitVector, keep1::BitVector;
                             bam::BAMParams=BAMParams(), mp::MovementParams=MovementParams())
    P0, B0, _, M0 = assemble_BAM(pool, grid, A, keep0; bam=bam, mp=mp)
    S0 = species_stats(pool, grid, A, keep0, P0, B0, M0, bam)
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

# ΔF vs ΔBSH (alignment)
function deltas_F_vs_BSH(pool::SpeciesPool, grid::Grid, A::Matrix{Float64},
                         base_keep::BitVector, keep::BitVector;
                         bam::BAMParams, mp::MovementParams)
    P0, B0, _, M0 = assemble_BAM(pool, grid, A, base_keep; bam=bam, mp=mp)
    S0 = species_stats(pool, grid, A, base_keep, P0, B0, M0, bam)
    P1, B1, _, M1 = assemble_BAM(pool, grid, A, keep; bam=bam, mp=mp)
    S1 = species_stats(pool, grid, A, keep, P1, B1, M1, bam)
    dF = S1.F .- S0.F
    dBSH = S1.BSH .- S0.BSH
    cons = consumer_mask(pool)
    (; dF, dBSH, cons)
end

# ╭─────────────────────────────
# │ 5) Axes → pool & BAM settings
# ╰─────────────────────────────
# Trim diets to reduce redundancy (used for B=:soft and B=:strong)
function _trim_diets!(rng::AbstractRNG, pool::SpeciesPool, max_prey::Int)
    for s in 1:pool.S
        pool.basal[s] && continue
        L = length(pool.E[s])
        if L > max_prey
            shuffle!(rng, pool.E[s])
            pool.E[s] = pool.E[s][1:max_prey]
        end
    end
    return pool
end

"""
Construct a pool consistent with abiotic level (A_level) and biotic level (B_level).

A_level ∈ {:neutral,:intermediate,:divergent}
B_level ∈ {:none,:soft,:strong}
"""
function build_pool_from_axes(
    rng::AbstractRNG; S::Int, basal_frac::Float64,
    A_level::Symbol, B_level::Symbol
)
    # ── Abiotic (niche) setup ──────────────────────────────────────────────
    niche_mode = :uniform
    mu_centers = (0.25, 0.75)
    mu_sd      = 0.05
    b0_basal, bspread_basal, b0_cons, bspread_cons = (0.12, 0.05, 0.12, 0.05)
    if A_level === :neutral
        niche_mode = :uniform
        b0_basal, bspread_basal = 0.22, 0.05
        b0_cons,  bspread_cons  = 0.22, 0.05
    elseif A_level === :intermediate
        niche_mode = :uniform
        b0_basal, bspread_basal = 0.12, 0.05
        b0_cons,  bspread_cons  = 0.12, 0.05
    elseif A_level === :divergent
        niche_mode = :bimodal
        mu_centers = (0.2, 0.8); mu_sd = 0.05
        b0_basal, bspread_basal = 0.08, 0.02
        b0_cons,  bspread_cons  = 0.10, 0.03
    else
        error("Unknown A_level $A_level")
    end

    # ── Biotic (metaweb) setup ─────────────────────────────────────────────
    # density ↓ ⇒ fewer links; pmax slight ↓ for strong
    density, pmax, diet_cap = 0.30, 0.90, typemax(Int)
    if B_level === :none
        density, pmax, diet_cap = 0.35, 0.90, typemax(Int)
    elseif B_level === :soft
        density, pmax, diet_cap = 0.24, 0.90, 6     # fewer links + cap 6
    elseif B_level === :strong
        density, pmax, diet_cap = 0.08, 0.85, 2     # much fewer + cap 3
    else
        error("Unknown B_level $B_level")
    end

    pool = build_pool(; S, rng, basal_frac,
                      niche_mode, mu_basal_centers=mu_centers, mu_basal_sd=mu_sd,
                      b0_basal, bspread_basal, b0_cons, bspread_cons,
                      density=density, pmax=pmax)

    # Cap diets (reduces redundancy where it matters)
    if isfinite(diet_cap)
        _trim_diets!(rng, pool, diet_cap)
    end
    return pool
end

"Translate B_level and M_level into BAM & movement parameters."
function bam_from_axes(; B_level::Symbol, M_level::Symbol,
                       α::Float64=1.0, τA::Float64=0.35, τocc::Float64=0.42)
    β, γ = 0.0, 2.0
    if B_level === :none
        β, γ = 0.0, 2.0                  # ignore B
    elseif B_level === :soft
        β, γ = 1.2, 3.0                  # matters modestly
    elseif B_level === :strong
        β, γ = 3.0, 7.0                  # sharp/strong dependence
    else
        error("Unknown B_level $B_level")
    end
    μ, mode = 0.0, :none
    if M_level === :off
        μ, mode = 0.0, :none
    elseif M_level === :on
        μ, mode = 0.8, :component        # movement contributes strongly
    else
        error("Unknown M_level $M_level")
    end
    (bam=BAMParams(; α=α, β=β, μ=μ, γ=γ, τA=τA, τocc=τocc),
     mp=MovementParams(; mode=mode, T=8))   # T overwritten later by scaling
end

# ╭─────────────────────────────
# │ 6) Threaded sweep for one HL
# ╰─────────────────────────────
function run_sweep_threaded_axes(; grid::Grid, S::Int=120, basal_frac::Float64=0.45,
        A_level::Symbol=:intermediate, B_level::Symbol=:soft, M_level::Symbol=:on,
        hl_kind::Symbol=:clustered, nseeds_cluster::Int=6,
        front_axis::Symbol=:x, front_noise::Float64=0.04,
        loss_fracs = 0.2:0.1:0.8,
        seeds_pool = 1:5, seeds_mask = 1:5,
        sim_seed::Int=1234, τA::Float64=0.5, τocc::Float64=0.35, T::Int=8)

    # BAM & movement
    pars = bam_from_axes(; B_level, M_level, τA, τocc)
    bam, mp = pars.bam, MovementParams(; mode=pars.mp.mode, T=T)

    xs = collect(loss_fracs)
    yΔ = fill(NaN, length(xs)); yA=similar(yΔ); yB=similar(yΔ); yM=similar(yΔ); ySyn=similar(yΔ)
    base_keep = trues(grid.C)

    Threads.@threads for k in eachindex(xs)
        f = xs[k]; keepfrac = 1 - f
        Δ  = Float64[]; AΔ = Float64[]; BΔ = Float64[]; MΔ = Float64[]; SΔ = Float64[]

        for ps in seeds_pool, ms in seeds_mask
            rng_pool = thread_rng(sim_seed, :pool, ps, k)
            rng_mask = thread_rng(sim_seed, :mask, ms, k, hl_kind)

            pool = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
            println("pool ", prey_stats(pool))
            A = abiotic_matrix(pool, grid)

            keep = hl_kind === :random    ? random_mask(rng_mask, grid.C, keepfrac) :
                   hl_kind === :clustered ? clustered_mask(rng_mask, grid, keepfrac; nseeds=nseeds_cluster) :
                   hl_kind === :front     ? frontlike_mask(rng_mask, grid, keepfrac; axis=front_axis, noise=front_noise) :
                   error("Unknown hl_kind = $hl_kind")

            shap, _, _ = shapley_per_species(pool, grid, A, base_keep, keep; bam=bam, mp=mp)
            cons = consumer_mask(pool)
            push!(Δ,  mean(getfield.(shap[cons], :dF)))
            push!(AΔ, mean(getfield.(shap[cons], :dA)))
            push!(BΔ, mean(getfield.(shap[cons], :dB)))
            push!(MΔ, mean(getfield.(shap[cons], :dM)))
            push!(SΔ, mean(getfield.(shap[cons], :synergy)))
        end

        yΔ[k]=mean(Δ); yA[k]=mean(AΔ); yB[k]=mean(BΔ); yM[k]=mean(MΔ); ySyn[k]=mean(SΔ)
    end

    return (x=xs, dF=yΔ, dA=yA, dB=yB, dM=yM, dSyn=ySyn,
            meta=Dict(:hl=>hl_kind, :A=>A_level, :B=>B_level, :M=>M_level, :T=>T))
end

# Shapley on F with A held at baseline (isolates B and M contributions)
function shapley_BM_condA(A0,B0,M0, A1,B1,M1; β, μ)
    f(A,B,M) = A0 * (B^β) * (M^μ)
    ΔB = (1/2)*(f(A0,B1,M0) - f(A0,B0,M0)) + (1/2)*(f(A0,B1,M1) - f(A0,B0,M1))
    ΔM = (1/2)*(f(A0,B0,M1) - f(A0,B0,M0)) + (1/2)*(f(A0,B1,M1) - f(A0,B1,M0))
    dF = (A0*(B1^β)*(M1^μ)) - (A0*(B0^β)*(M0^μ))
    (; dF, ΔB, ΔM, synergy=dF-(ΔB+ΔM))
end

# alignment for one HL (single pool+mask for clarity)
function one_alignment_axes(; grid::Grid, S::Int=120, basal_frac::Float64=0.45,
        A_level::Symbol=:intermediate, B_level::Symbol=:soft, M_level::Symbol=:on,
        hl_kind::Symbol=:clustered, keep_frac::Float64=0.5,
        sim_seed::Int=1234, pool_seed::Int=1, mask_seed::Int=1,
        τA::Float64=0.5, τocc::Float64=0.35, T::Int=8,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04)

    pars = bam_from_axes(; B_level, M_level, τA, τocc)
    bam, mp = pars.bam, MovementParams(; mode=pars.mp.mode, T=T)

    rng_pool = MersenneTwister(hash((sim_seed,:pool,pool_seed)))
    pool     = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
    A        = abiotic_matrix(pool, grid)
    base_keep = trues(grid.C)

    rng_mask = MersenneTwister(hash((sim_seed,:mask,mask_seed,hl_kind)))
    keep = hl_kind === :random    ? random_mask(rng_mask, grid.C, keep_frac) :
           hl_kind === :clustered ? clustered_mask(rng_mask, grid, keep_frac; nseeds=nseeds_cluster) :
           hl_kind === :front     ? frontlike_mask(rng_mask, grid, keep_frac; axis=front_axis, noise=front_noise) :
           error("Unknown HL $(hl_kind)")

    d = deltas_F_vs_BSH(pool, grid, A, base_keep, keep; bam=bam, mp=mp)
    x = d.dF[d.cons]; y = d.dBSH[d.cons]
    r = (std(x)>0 && std(y)>0) ? cor(x,y) : NaN
    slope = sum((x .- mean(x)).*(y .- mean(y))) / sum((x .- mean(x)).^2)
    (; x, y, r, slope)
end

# ╭─────────────────────────────
# │ 7) Plot: top (ΔF components), bottom (alignment)
# ╰─────────────────────────────
function run_combo_with_alignment(; grid::Grid, S::Int=120, basal_frac::Float64=0.45,
        A_level::Symbol=:intermediate, B_level::Symbol=:soft, M_level::Symbol=:on,
        loss_fracs = 0.2:0.1:0.8, align_keep_frac::Float64=0.5,
        seeds_pool = 1:5, seeds_mask = 1:5,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        τA::Float64=0.5, τocc::Float64=0.35, T::Int=8, sim_seed::Int=1234)

    results = Dict{Symbol,Any}()
    for hk in (:random, :clustered, :front)
        results[hk] = run_sweep_threaded_axes(; grid,S,basal_frac,A_level,B_level,M_level,
                    hl_kind=hk, nseeds_cluster, front_axis, front_noise,
                    loss_fracs, seeds_pool, seeds_mask, sim_seed, τA, τocc, T)
    end

    A_rand = one_alignment_axes(; grid,S,basal_frac,A_level,B_level,M_level,
                                hl_kind=:random, keep_frac=align_keep_frac, τA, τocc, T)
    A_cl   = one_alignment_axes(; grid,S,basal_frac,A_level,B_level,M_level,
                                hl_kind=:clustered, keep_frac=align_keep_frac, τA, τocc, T)
    A_fr   = one_alignment_axes(; grid,S,basal_frac,A_level,B_level,M_level,
                                hl_kind=:front, keep_frac=align_keep_frac, τA, τocc, T)

    # Figure
    fig = Figure(; size=(1260, 760))
    labels = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")
    Label(fig[0,1:3], "A_$(A_level), B_$(B_level), M_$(M_level)", fontsize=20)
    # top row
    for (col, hk) in enumerate((:random,:clustered,:front))
        r = results[hk]
        ax = Axis(fig[1,col], xlabel="Area lost (fraction)", ylabel="ΔF (mean cons.)",
                  title=labels[hk])
        lines!(ax, r.x, r.dF,  label="Total ΔF")
        lines!(ax, r.x, r.dA,  label="Abiotic")
        lines!(ax, r.x, r.dB,  label="Biotic")
        lines!(ax, r.x, r.dM,  label="Movement")
        lines!(ax, r.x, r.dSyn,label="Synergy")
        if col==1; axislegend(ax; position=:lt, framevisible=false); end
    end

    # bottom row
    for (col,(hk,Ahl)) in enumerate(zip((:random,:clustered,:front),(A_rand,A_cl,A_fr)))
        ax = Axis(fig[2,col], xlabel="ΔF per species", ylabel=(col==1 ? "ΔBSH per species" : ""),
                  title = labels[hk]*@sprintf("  (r=%.2f, slope=%.2f)", Ahl.r, Ahl.slope))
        scatter!(ax, Ahl.x, Ahl.y; markersize=5)
        lo = minimum(vcat(Ahl.x, Ahl.y)); hi = maximum(vcat(Ahl.x, Ahl.y))
        lines!(ax, [lo, hi], [lo, hi]; linestyle=:dash)
    end
    display(fig)

    return (res_by_hl=results, fig=fig)
end

# ╭─────────────────────────────
# │ 8) Named combination dictionary (18 combos)
# ╰─────────────────────────────
const A_LEVELS = (:neutral, :intermediate, :divergent)
const B_LEVELS = (:none, :soft, :strong)
const M_LEVELS = (:off, :on)

"Create combo names like :A_divergent__B_strong__M_on → (A,B,M)."
function build_combo_dict()
    D = Dict{Symbol,NamedTuple}()
    for A in A_LEVELS, B in B_LEVELS, M in M_LEVELS
        name = Symbol("A_", String(A), "__B_", String(B), "__M_", String(M))
        D[name] = (; A_level=A, B_level=B, M_level=M)
    end
    D
end
const COMBOS = build_combo_dict()

"Run a combo by name; change only the name to switch regimes."
function run_combo_name(name::Symbol; nx::Int=60, ny::Int=60, kwargs...)
    @assert haskey(COMBOS, name) "Unknown combo $(name). Try one of: $(collect(keys(COMBOS)))"
    grid = make_grid(nx, ny; seed=42)
    axes = COMBOS[name]
    run_combo_with_alignment(; grid, axes..., kwargs...)
end