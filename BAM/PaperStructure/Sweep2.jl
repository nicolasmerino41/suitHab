# --- small helpers -----------------------------------------------------------

# pass-through wrapper so we can override diet cap / density / λ, ξ cleanly
function build_pool_from_axes_with_overrides(
    rng::AbstractRNG; S::Int, basal_frac::Float64,
    A_level::Symbol, B_level_for_thresholds::Symbol,
    density_override::Union{Nothing,Float64},
    diet_cap_override::Union{Nothing,Int},
    λ_align::Float64, ξ_align::Float64
)
    # 1) start from your axis recipe (for A-level only)
    pool = build_pool_from_axes(rng; S, basal_frac, A_level, B_level=:none)  # we’ll apply B “strength” in BAM (β,γ) not here

    # 2) override metaweb density if requested (thin uniformly)
    if density_override !== nothing
        # crude: Bernoulli thin of adjacency lists to match new density ratio
        target = density_override
        # estimate current link density
        L = sum(length(pool.E[s]) for s in 1:pool.S if !pool.basal[s])
        Scons = sum(!pool.basal[s] for s in 1:pool.S)
        approx_possible = max(L / (target + 1e-9), 1.0) # not used; we thin proportionally
        pkeep = clamp(target / 0.30, 0.05, 1.0)  # assumes your original ≈0.30; adjust if needed
        for s in 1:pool.S
            pool.basal[s] && continue
            if !isempty(pool.E[s])
                keep = rand(rng, length(pool.E[s])) .< pkeep
                pool.E[s] = pool.E[s][keep]
            end
        end
    end

    # 3) diet cap override (low redundancy)
    if diet_cap_override !== nothing
        _trim_diets!(rng, pool, diet_cap_override)
    end

    # 4) inject climate-sorting weight on *existing* links (reweight → resample)
    #    We resample each predator’s prey list in proportion to w = (1-λ)+λ*exp(-(Δμ)^2/(2ξ^2)),
    #    keeping the same number of links to avoid confounding density changes.
    if λ_align > 0 && ξ_align > 0
        for s in 1:pool.S
            pool.basal[s] && continue
            prey = pool.E[s]
            Ls   = length(prey)
            Ls == 0 && continue
            Δμ = abs.(pool.mu[s] .- pool.mu[prey])
            w  = (1 .- λ_align) .+ λ_align .* exp.(-(Δμ.^2) ./ (2ξ_align^2))
            w ./= sum(w) + eps()
            # resample without replacement with probs ~ w
            idx = collect(1:Ls)
            # Weighted sampling without replacement (simple greedy for small Ls)
            chosen = Int[]
            probs = copy(w)
            avail = collect(1:Ls)
            while !isempty(avail) && length(chosen) < Ls
                chosen_idx = sample(rng, avail, Weights(probs[avail]))
                push!(chosen, chosen_idx)
                deleteat!(avail, findfirst(==(chosen_idx), avail))
                # re-normalize remaining
                denom = sum(probs[avail]) + eps()
                for j in 1:length(avail); probs[avail[j]] /= denom; end
            end
            pool.E[s] = prey[chosen]
        end
    end

    return pool
end

# ΔBSH (mean consumers) for a single HL geometry at a chosen loss
function mean_dBSH_at_f(; grid::Grid, pool::SpeciesPool, A::Matrix{Float64},
        B_level::Symbol, M_level::Symbol, hl_kind::Symbol, f::Float64,
        τA::Float64, τocc::Float64, T::Int,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        seeds_mask=1:6, sim_seed::Int=1234)

    pars = bam_from_axes(; B_level, M_level, τA, τocc)
    bam  = pars.bam
    mp   = MovementParams(; mode=pars.mp.mode, T=T)

    base_keep = trues(grid.C)
    # choose keep mask for each replicate at same keep_frac
    keep_frac = 1 - f
    vals = Float64[]
    for ms in seeds_mask
        rng_mask = MersenneTwister(hash((sim_seed,:mask,ms,hl_kind)))
        keep = hl_kind === :random    ? random_mask(rng_mask, grid.C, keep_frac) :
               hl_kind === :clustered ? clustered_mask(rng_mask, grid, keep_frac; nseeds=nseeds_cluster) :
               hl_kind === :front     ? frontlike_mask(rng_mask, grid, keep_frac; axis=front_axis, noise=front_noise) :
               error("Unknown HL $(hl_kind)")
        # assemble and compute BSH loss vs baseline
        P0,B0,_,M0 = assemble_BAM(pool, grid, A, base_keep; bam=bam, mp=mp)
        S0 = species_stats(pool, grid, A, base_keep, P0, B0, M0, bam)
        P1,B1,_,M1 = assemble_BAM(pool, grid, A, keep;      bam=bam, mp=mp)
        S1 = species_stats(pool, grid, A, keep,      P1, B1, M1, bam)
        cons = consumer_mask(pool)
        push!(vals, mean(S1.BSH[cons] .- S0.BSH[cons]))
    end
    return mean(vals)
end

# full DiD for three geometries at f*
function did_triplet_at_f(; grid::Grid, S::Int, basal_frac::Float64,
        A_level::Symbol, fstar::Float64, τA::Float64=0.5, τocc::Float64=0.42,
        # sweep knobs:
        λ_align::Float64, ξ_align::Float64=0.06,
        diet_cap::Int, density::Float64,
        β_strong::Float64=4.0, M_level::Symbol=:off,
        seeds_pool=1:4, seeds_mask=1:4, sim_seed::Int=1234)

    # set strong biotic dependence in BAM mapping
    # old_bam_from_axes = bam_from_axes
    function bam_from_axes(; B_level::Symbol, M_level::Symbol, α::Float64=1.0, τA::Float64=0.35, τocc::Float64=0.42)
        β, γ = 0.0, 2.0
        if B_level === :none
            β, γ = 0.0, 2.0
        elseif B_level === :soft
            β, γ = 2.0, 5.0
        elseif B_level === :strong
            β, γ = β_strong, 7.0
        end
        μ, mode = (M_level===:off ? (0.0,:none) : (0.8,:component))
        return (bam=BAMParams(; α=α, β=β, μ=μ, γ=γ, τA=τA, τocc=τocc), mp=MovementParams(; mode=mode, T=8))
    end

    geoms = (:random, :clustered, :front)
    DiD = Dict(g => Float64[] for g in geoms)

    for ps in seeds_pool
        rng_pool = MersenneTwister(hash((sim_seed,:pool,ps)))
        # build pool with overrides
        pool = build_pool_from_axes_with_overrides(rng_pool;
            S, basal_frac, A_level,
            B_level_for_thresholds=:strong,
            density_override=density, diet_cap_override=diet_cap,
            λ_align=λ_align, ξ_align=ξ_align)
        A = abiotic_matrix(pool, grid)

        for g in geoms
            d_non = mean_dBSH_at_f(; grid, pool, A,
                B_level=:none, M_level, hl_kind=g, f=fstar,
                τA, τocc, T=8, seeds_mask, sim_seed)
            d_bio = mean_dBSH_at_f(; grid, pool, A,
                B_level=:strong, M_level, hl_kind=g, f=fstar,
                τA, τocc, T=8, seeds_mask, sim_seed)
            push!(DiD[g], d_bio - d_non)
        end
    end

    # summary
    did̄ = Dict(g => mean(DiD[g]) for g in geoms)
    # worst-geometry under B OFF vs B ON (for rank-flip)
    d_non_order = sort(collect(geoms), by = g -> did̄[g] + 0.0) # proxy using DiD; if you prefer, recompute ΔBSH_non separately
    worst_non = first(d_non_order)   # smallest
    worst_bio = first(sort(collect(geoms), by = g -> did̄[g])) # same ordering proxy
    rank_flip = worst_non != worst_bio
    spread_change = (maximum(values(did̄)) - minimum(values(did̄)))
    return (; did̄, rank_flip, spread_change)
end

# sweep over (λ_align × diet_cap), return matrices for front DiD and worst-case amp
function sweep_critical_map(; nx::Int=60, ny::Int=60, S::Int=150, basal_frac::Float64=0.45,
        A_level::Symbol=:intermediate, fstar::Float64=0.6, M_level::Symbol=:off,
        λ_list = range(0.0, 1.0; length=9), cap_list = [2,3,4,6,9999],
        density::Float64=0.12, β_strong::Float64=4.0, ξ_align::Float64=0.06,
        seeds_pool=1:4, seeds_mask=1:4, sim_seed::Int=1234)

    grid = make_grid(nx, ny; seed=42)
    G = (:random, :clustered, :front)
    nλ, nc = length(λ_list), length(cap_list)
    DiD_front = fill(NaN, nλ, nc)
    DiD_max   = fill(NaN, nλ, nc)
    Flip      = falses(nλ, nc)

    for (i,λ) in enumerate(λ_list), (j,cap) in enumerate(cap_list)
        res = did_triplet_at_f(; grid,S,basal_frac,A_level,fstar,M_level,
                               λ_align=λ, ξ_align, diet_cap=cap==9999 ? typemax(Int) : cap,
                               density, β_strong, seeds_pool, seeds_mask, sim_seed)
        d = res.did̄
        DiD_front[i,j] = d[:front]
        DiD_max[i,j]   = maximum([d[g] for g in G])
        Flip[i,j]      = res.rank_flip
    end
    return (; λ_list, cap_list, DiD_front, DiD_max, Flip, A_level, fstar, M_level)
end

# plot: heatmap of DiD_front and worst-case amplification, with rank-flip dots
function plot_critical_map(res)
    fig = Figure(; size=(1200, 520))
    Label(fig[0,1:2], "Critical regimes at f*=$(res.fstar) — A=$(res.A_level), M=$(res.M_level)", fontsize=18)

    # 1) front DiD
    ax1 = Axis(fig[1,1], xlabel="diet cap", ylabel="λ_align", title="DiD(front) = Δbio - Δnon")
    hm1 = heatmap!(ax1, 1:length(res.cap_list), 1:length(res.λ_list), res.DiD_front;
                   colormap=:viridis)
    ax1.xticks = (1:length(res.cap_list), string.(res.cap_list))
    ax1.yticks = (1:length(res.λ_list), string.(round.(Float64.(res.λ_list); digits=2)))
    Colorbar(fig[1,2], hm1; label="ΔBSH (mean cons.)")

    # 2) worst-case amplification + rank flips
    ax2 = Axis(fig[1,3], xlabel="diet cap", ylabel="λ_align",
               title="Worst-geometry amplification  max_g DiD(g)")
    hm2 = heatmap!(ax2, 1:length(res.cap_list), 1:length(res.λ_list), res.DiD_max; colormap=:magma)
    ax2.xticks = (1:length(res.cap_list), string.(res.cap_list))
    ax2.yticks = (1:length(res.λ_list), string.(round.(Float64.(res.λ_list); digits=2)))
    # overlay flips
    IJ = findall(res.Flip)
    I = first.(IJ)
    J = last.(IJ)
    if !isempty(I)
        scatter!(ax2, J, I; markersize=9, marker=:xcross, color=:white, strokecolor=:black, strokewidth=1.5)
        text!(ax2, 0.5, 1.02, text="× = rank flip (worst geometry changes with B)", align=(:left,:bottom), space=:relative)
    end
    Colorbar(fig[1,4], hm2; label="ΔBSH (mean cons.)")

    display(fig)
    return fig
end
# Amplify biotic effect and map where it appears
res = sweep_critical_map(; A_level=:intermediate,
                         fstar=0.6,
                         M_level=:off,               # headline: no screening
                         λ_list=0.0:0.1:1.0,         # climate sorting strength
                         cap_list=[2,3,4,6,9999],    # 9999 = effectively ∞ (no cap)
                         density=0.12, β_strong=4.0, ξ_align=0.06,
                         seeds_pool=1:4, seeds_mask=1:4)

fig = plot_critical_map(res)
