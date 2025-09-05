# ---------- thread-local RNG ----------
@inline function thread_rng(seed::Int, tags...)::MersenneTwister
    return MersenneTwister(hash((seed, tags..., Threads.threadid())))
end

# ---------- RANDOM MASKS (thread-safe) ----------
function random_mask(rng::AbstractRNG, C::Int, keep_frac::Float64)
    nkeep = clamp(round(Int, keep_frac * C), 0, C)
    keep = falses(C)
    idx  = randperm(rng, C)[1:nkeep]
    keep[idx] .= true
    keep
end

# 4-neighborhood you already have; kept for completeness
neighbors4(ix::Int, nx::Int, ny::Int) = begin
    i = ((ix - 1) % nx) + 1
    j = ((ix - 1) ÷ nx) + 1
    out = Int[]
    if i > 1; push!(out, ix - 1); end
    if i < nx; push!(out, ix + 1); end
    if j > 1; push!(out, ix - nx); end
    if j < ny; push!(out, ix + nx); end
    out
end

function clustered_mask(rng::AbstractRNG, grid::Grid, keep_frac::Float64; nseeds::Int=6)
    C = grid.C
    target_remove = C - clamp(round(Int, keep_frac*C), 0, C)
    removed = falses(C)
    q = Int[]

    seeds = randperm(rng, C)[1:min(nseeds, C)]
    append!(q, seeds)

    removed_ct = 0
    ptr = 1
    while removed_ct < target_remove
        if ptr > length(q)
            # start a new seed in a not-removed cell
            while true
                cand = rand(rng, 1:C)
                if !removed[cand]; push!(q, cand); break; end
            end
        end
        v = q[ptr]; ptr += 1
        if removed[v]; continue; end
        removed[v] = true
        removed_ct += 1
        for nb in neighbors4(v, grid.nx, grid.ny)
            if !removed[nb]; push!(q, nb); end
        end
    end
    .!removed
end

# Front-like mask (axis-aligned with optional fuzz)
function frontlike_mask(rng::AbstractRNG, grid::Grid, keep_frac::Float64; axis::Symbol=:x, noise::Float64=0.0)
    C = grid.C
    coord = axis === :x ? view(grid.xy, 1, :) : view(grid.xy, 2, :)
    q = quantile(coord, 1 - keep_frac)
    keep = BitVector(undef, C)
    if noise ≤ 0
        @inbounds for i in 1:C
            keep[i] = coord[i] ≥ q
        end
    else
        ϵ = rand(rng, C) .- 0.5
        @inbounds for i in 1:C
            keep[i] = coord[i] ≥ (q + noise*ϵ[i])
        end
    end
    keep
end

# ---------- build_pool with rng (thread-safe) ----------
function build_pool(
    S::Int;
    rng::AbstractRNG,
    basal_frac::Float64 = 0.45,
    # diet / redundancy
    R0_mean::Float64 = 12.0,
    R0_sd::Float64   = 0.50,
    sigma::Float64   = 0.50,
    density::Float64 = 0.30,
    pmax::Float64    = 0.90,
    # climate / synchrony
    niche_mode::Symbol = :uniform,
    mu_basal_centers::Tuple{Float64,Float64} = (0.25, 0.75),
    mu_basal_sd::Float64 = 0.05,
    b0_basal::Float64 = 0.12, bspread_basal::Float64 = 0.05,
    b0_cons::Float64  = 0.12, bspread_cons::Float64  = 0.05
)
    logm   = collect(range(log(1e-2), log(10.0); length=S))
    shuffle!(rng, logm)
    masses = exp.(logm)

    order = sortperm(masses)
    nB    = clamp(round(Int, basal_frac * S), 0, S)
    basal = falses(S); basal[order[1:nB]] .= true

    mu = similar(masses)
    b  = similar(masses)

    cons_ids = findall(!, basal)
    mu[cons_ids] .= rand(rng, length(cons_ids))
    b[cons_ids]  .= b0_cons .+ bspread_cons .* rand(rng, length(cons_ids))

    bas_ids = findall(basal)
    if niche_mode === :bimodal
        nb  = length(bas_ids)
        nb1 = nb ÷ 2; nb2 = nb - nb1
        c1, c2 = mu_basal_centers
        mu[bas_ids[1:nb1]]     .= clamp.(c1 .+ mu_basal_sd .* randn(rng, nb1), 0, 1)
        mu[bas_ids[nb1+1:end]] .= clamp.(c2 .+ mu_basal_sd .* randn(rng, nb2), 0, 1)
    else
        mu[bas_ids] .= rand(rng, length(bas_ids))
    end
    b[bas_ids] .= b0_basal .+ bspread_basal .* rand(rng, length(bas_ids))

    R0 = exp.(log(R0_mean) .+ R0_sd .* randn(rng, S))

    E = [Int[] for _ in 1:S]
    for (ii, s) in pairs(order)
        basal[s] && continue
        for jj in 1:ii-1
            q = order[jj]
            r = masses[s] / masses[q]
            z = (log(r) - log(R0[s])) / sigma
            p = pmax * exp(-0.5 * z^2) * density
            if rand(rng) < p; push!(E[s], q); end
        end
        if isempty(E[s]) && ii > 1
            cand   = order[1:ii-1]
            target = log(masses[s]) - log(R0[s])
            qstar  = cand[argmin(abs.(log.(masses[cand]) .- target))]
            push!(E[s], qstar)
        end
    end

    SpeciesPool(S, masses, basal, mu, b, E)
end

"""
run_sweep_threaded(...; hl_kind=:random | :clustered | :front)

Threaded over loss_fracs. Uses per-thread RNGs for pools and masks.
Reproducible via sim_seed.
"""
function run_sweep_threaded(; grid::Grid,
        S::Int=100, basal_frac::Float64=0.45,
        τA::Float64=0.5, τocc::Float64=0.35, α::Float64=1.0, β::Float64=1.0, μ::Float64=0.0, γ::Float64=2.0,
        move_mode::Symbol=:none, λ::Float64=0.12, T::Int=4,
        hl_kind::Symbol=:clustered, nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.05,
        loss_fracs = 0.1:0.1:0.8,
        seeds_pool = 1:5, seeds_mask = 1:5,
        sim_seed::Int = 1234)

    bam = BAMParams(α=α, β=β, μ=μ, γ=γ, τA=τA, τocc=τocc)
    mp  = MovementParams(mode=move_mode, λ=λ, T=T)

    n = length(loss_fracs)
    xs   = collect(loss_fracs)
    yΔ   = fill(NaN, n); yA  = similar(yΔ); yB = similar(yΔ); yM = similar(yΔ); ySyn = similar(yΔ)

    base_keep = trues(grid.C)

    Threads.@threads for k in eachindex(xs)
        f = xs[k]; keepfrac = 1 - f
        Δ  = Float64[]; AΔ = Float64[]; BΔ = Float64[]; MΔ = Float64[]; SΔ = Float64[]

        for ps in seeds_pool, ms in seeds_mask
            rng_pool = thread_rng(sim_seed, :pool, ps, k)
            rng_mask = thread_rng(sim_seed, :mask, ms, k, hl_kind)

            pool = build_pool(S; rng=rng_pool, basal_frac=basal_frac)
            A    = abiotic_matrix(pool, grid)

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

        yΔ[k]   = mean(Δ);   yA[k] = mean(AΔ);  yB[k] = mean(BΔ);  yM[k] = mean(MΔ);  ySyn[k] = mean(SΔ)
    end

    return (x=xs, dF=yΔ, dA=yA, dB=yB, dM=yM, dSyn=ySyn,
            meta=Dict(:hl=>hl_kind, :α=>α, :β=>β, :μ=>μ, :move=>move_mode, :λ=>λ, :T=>T))
end

function plot_three_geometries(res_by_hl::Dict{Symbol,Any}; title="BAM × HL")
    begin
        fig = Figure(; size = (1200, 380))
        hl_order = (:random, :clustered, :front)
        labels   = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")

        for (col, hk) in enumerate(hl_order)
            r = res_by_hl[hk]
            ax = Axis(fig[1,col], xlabel="Area lost (fraction)", ylabel="ΔF (mean cons.)",
                      title="$(labels[hk])")
            lines!(ax, r.x, r.dF,  label="Total ΔF")
            lines!(ax, r.x, r.dA,  label="Abiotic")
            lines!(ax, r.x, r.dB,  label="Biotic")
            lines!(ax, r.x, r.dM,  label="Movement")
            lines!(ax, r.x, r.dSyn,label="Synergy")
            if col == 1
                axislegend(ax; position=:lt, framevisible=false)
            end
        end
        Label(fig[0, :], title; fontsize=16, padding=(0,0,10,0))
        display(fig)
    end
end

const BAM_COMBOS = Dict(
    :baseline => (; α=1.0, β=1.0, μ=0.0, γ=2.0, move_mode=:none,    λ=0.0,  T=4,
                   title="Baseline: β=1, μ=0, no movement"),
    :strongB_access => (; α=1.0, β=1.8, μ=0.6, γ=3.0, move_mode=:access, λ=0.08, T=4,
                         title="Strong biotic + accessibility"),
    :threshold_component => (; α=1.0, β=2.2, μ=0.8, γ=4.0, move_mode=:component, λ=0.0,  T=8,
                              title="Thresholdy biotic + component req.")
)

function run_combo_and_plot(combo::Symbol;
        nx::Int=60, ny::Int=50, S::Int=120, basal_frac::Float64=0.45,
        loss_fracs = 0.2:0.1:0.8,
        seed_grid::Int=42, seeds_pool = 1:5, seeds_mask = 1:5,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        sim_seed::Int=1234)

    @assert haskey(BAM_COMBOS, combo) "Unknown combo $(combo). Keys: $(keys(BAM_COMBOS))"
    p = BAM_COMBOS[combo]

    grid = make_grid(nx, ny; seed=seed_grid)  # your existing function (deterministic)

    res = Dict{Symbol,Any}()
    for hk in (:random, :clustered, :front)
        res[hk] = run_sweep_threaded(; grid, S, basal_frac,
            τA=0.5, τocc=0.35,
            α=p.α, β=p.β, μ=p.μ, γ=p.γ,
            move_mode=p.move_mode, λ=p.λ, T=p.T,
            hl_kind=hk, nseeds_cluster=nseeds_cluster, front_axis=front_axis, front_noise=front_noise,
            loss_fracs=loss_fracs,
            seeds_pool=seeds_pool, seeds_mask=seeds_mask,
            sim_seed=sim_seed)
    end

    plot_three_geometries(res; title=p.title)
    return res
end

# 1) Run ONE combo (and get Random vs Clustered vs Front plots)
results = run_combo_and_plot(:strongB_access;
    nx=60, ny=60, S=150,
    loss_fracs=0.2:0.1:0.8,
    seeds_pool=1:6, seeds_mask=1:6)

# 2) Switch combo by changing just the name:
results2 = run_combo_and_plot(
    :threshold_component;
    nx=60, ny=60, S=150,
    loss_fracs=0.2:0.1:0.8,
    seeds_pool=1:6, seeds_mask=1:6
)
results3 = run_combo_and_plot(
    :baseline;
    nx=60, ny=60, S=150,
    loss_fracs=0.2:0.1:0.8,
    seeds_pool=1:6, seeds_mask=1:6
)
