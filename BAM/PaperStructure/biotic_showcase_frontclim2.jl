"""
build_pool_aligned(...)

As `build_pool_from_axes`, but adds climate-assortative feeding and a hard diet cap.
Parameters:
  λ_align ∈ [0,1]  : weight of climate similarity (0=no effect, 1=full)
  ξ_align          : width of the climate kernel (smaller = more selective)
  density          : base thinning of links (as before)
  diet_cap::Int    : max prey per consumer (e.g., 2 or 3); use typemax(Int) for none
"""
function build_pool_aligned(
    rng::AbstractRNG; S::Int, basal_frac::Float64,
    A_level::Symbol,                     # use your existing niche presets
    density::Float64=0.10, pmax::Float64=0.90,
    λ_align::Float64=0.8, ξ_align::Float64=0.06,
    R0_mean::Float64=12.0, R0_sd::Float64=0.5, sigma::Float64=0.5,
    diet_cap::Int=2
)
    # 1) start from your A-level recipe
    pool = build_pool_from_axes(rng; S, basal_frac, A_level, B_level=:none)

    # 2) rebuild adjacency with climate-weighted acceptance
    order  = sortperm(pool.masses)                 # light→heavy
    R0     = exp.(log(R0_mean) .+ R0_sd .* randn(rng, S))
    Enew   = [Int[] for _ in 1:S]

    for (ii, s) in pairs(order)
        pool.basal[s] && continue
        for jj in 1:ii-1
            q = order[jj]
            # mass ratio kernel
            r = pool.masses[s] / pool.masses[q]
            z = (log(r) - log(R0[s])) / sigma
            pmass = pmax * exp(-0.5*z^2)
            # climate similarity kernel
            dμ = abs(pool.mu[s] - pool.mu[q])
            pclim = exp(-(dμ^2) / (2*ξ_align^2))^λ_align
            p = pmass * pclim * density
            if rand(rng) < p
                push!(Enew[s], q)
            end
        end
        # fallback: ensure at least 1 link
        if isempty(Enew[s]) && ii>1
            qstar = order[argmin(abs.(log.(pool.masses[order[1:ii-1]]) .- (log(pool.masses[s]) - log(R0[s]))))]
            push!(Enew[s], qstar)
        end
        # cap redundancy
        if length(Enew[s]) > diet_cap
            shuffle!(rng, Enew[s])
            Enew[s] = Enew[s][1:diet_cap]
        end
    end

    pool = SpeciesPool(pool.S, pool.masses, pool.basal, pool.mu, pool.b, Enew)
    return pool
end

"""
assemble_BAM_hard(...; bmode=:soft, θB=0.65, kB=2)

bmode = :soft      → b = 1 - exp(-γ x)               (your current)
      = :frac      → b = 1[x ≥ θB]                    (fraction of prey present)
      = :kofn      → b = 1[#present ≥ kB]             (k of n prey present)
Everything else matches `assemble_BAM`.
"""
function assemble_BAM_hard(pool::SpeciesPool, grid::Grid, A::Matrix{Float64}, keep::BitVector;
                           bam::BAMParams=BAMParams(), mp::MovementParams=MovementParams(),
                           bmode::Symbol=:soft, θB::Float64=0.65, kB::Int=2)
    S, C = pool.S, grid.C
    M = movement_matrix(pool, grid, A, keep; mp=mp, τA=bam.τA)
    P = falses(S, C)
    Bsup = ones(Float64, S, C)
    order = sortperm(pool.masses)
    W = prey_weights(pool)

    @inbounds for s in order
        prey = pool.E[s]
        for i in 1:C
            keep[i] || continue
            b = if pool.basal[s] || isempty(prey)
                1.0
            else
                # compute support
                present = 0
                for q in prey
                    present += P[q,i] ? 1 : 0
                end
                if bmode === :soft
                    x = length(prey)==0 ? 0.0 : present/length(prey)
                    1.0 - exp(-bam.γ * x)
                elseif bmode === :frac
                    x = length(prey)==0 ? 0.0 : present/length(prey)
                    x ≥ θB ? 1.0 : 0.0
                elseif bmode === :kofn
                    present ≥ kB ? 1.0 : 0.0
                else
                    error("Unknown bmode=$bmode")
                end
            end
            Bsup[s,i] = b
            score = (A[s,i]^bam.α) * (b^bam.β) * (M[s,i]^bam.μ)
            P[s,i] = score ≥ bam.τocc
        end
    end
    return P, Bsup, M
end

# Keep the warm (tail=:high) or cool (tail=:low) climate tail, with optional noise
function front_climate_mask(rng::AbstractRNG, grid::Grid, keep_frac::Float64;
                            tail::Symbol=:high, noise::Float64=0.0)
    C     = grid.C
    clim  = grid.climate
    q     = tail === :high ? quantile(clim, 1 - keep_frac) : quantile(clim, keep_frac)
    keep  = BitVector(undef, C)
    if noise ≤ 0
        if tail === :high
            @inbounds for i in 1:C; keep[i] = clim[i] ≥ q; end
        else
            @inbounds for i in 1:C; keep[i] = clim[i] ≤ q; end
        end
    else
        ϵ = noise .* (rand(rng, C) .- 0.5)
        if tail === :high
            @inbounds for i in 1:C; keep[i] = clim[i] + ϵ[i] ≥ q; end
        else
            @inbounds for i in 1:C; keep[i] = clim[i] + ϵ[i] ≤ q; end
        end
    end
    keep
end

function show_biotic_showcase(; nx=60, ny=60, S=140, basal_frac=0.45,
    A_level::Symbol=:divergent,
    λ_align=0.9, ξ_align=0.06, diet_cap=2, density=0.08,
    β=4.5, τA=0.48, τocc=0.45, M_level::Symbol=:off,
    loss_fracs=0.2:0.1:0.8, seeds_pool=1:3, seeds_mask=1:4,
    nseeds_cluster=6, front_axis=:x, front_noise=0.0, sim_seed=1234,
    bmode::Symbol=:frac, θB=0.65, kB=2)

    grid = make_grid(nx, ny; seed=42)
    rngp  = MersenneTwister(hash((sim_seed,:pool,first(seeds_pool))))

    # climate-aligned, low-redundancy pool
    pool = build_pool_aligned(rngp; S, basal_frac, A_level,
                              density, λ_align, ξ_align, diet_cap)
    A = abiotic_matrix(pool, grid)

    # helper
    function mean_dBSH(B_level, hl_kind)
        base_keep = trues(grid.C)
        pars = bam_from_axes(; B_level, M_level, τA, τocc)
        bam = BAMParams(; α=1.0, β=β, μ=pars.bam.μ, γ=pars.bam.γ, τA=τA, τocc=τocc)
        mp  = MovementParams(; mode=pars.mp.mode, T=8)

        vals = Float64[]
        for f in loss_fracs
            keepfrac = 1 - f
            Δs = Float64[]
            for ms in seeds_mask
                rng_mask = MersenneTwister(hash((sim_seed,:mask,ms,hl_kind)))
                keep = hl_kind === :random      ? random_mask(rng_mask, grid.C, keepfrac) :
                    hl_kind === :clustered   ? clustered_mask(rng_mask, grid, keepfrac; nseeds=nseeds_cluster) :
                    hl_kind === :front       ? frontlike_mask(rng_mask, grid, keepfrac; axis=front_axis, noise=front_noise) :
                    hl_kind === :front_clim  ? front_climate_mask(rng_mask, grid, keepfrac; tail=:high, noise=front_noise) :
                    error("Unknown hl_kind=$hl_kind")
                P0,B0,_ = assemble_BAM_hard(pool, grid, A, base_keep; bam=bam, mp=mp, bmode, θB, kB)
                S0 = species_stats(pool, grid, A, base_keep, P0, B0, ones(Float64,size(A)), bam)
                P1,B1,_ = assemble_BAM_hard(pool, grid, A, keep;      bam=bam, mp=mp, bmode, θB, kB)
                S1 = species_stats(pool, grid, A, keep,      P1, B1, ones(Float64,size(A)), bam)
                push!(Δs, mean(S1.BSH[.!pool.basal] .- S0.BSH[.!pool.basal]))
            end
            push!(vals, mean(Δs))
        end
        return vals
    end

    geoms = (:random,:clustered,:front)
    curves_non = Dict(g=>mean_dBSH(:none, g) for g in geoms)
    curves_bio = Dict(g=>mean_dBSH(:strong, g) for g in geoms)

    # plot
    fig = Figure(; size=(1400, 480))
    Label(fig[0,1:3], "HL effect with vs without biotic — A=$(A_level), M=$(M_level)"; fontsize=18)

    lab = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")
    cols = Dict(:random=>:dodgerblue3, :clustered=>:darkorange3, :front=>:seagreen4)

    ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="ΔBSH (mean consumers)",
               title="Climate-only (B OFF)")
    for g in geoms; lines!(ax1, collect(loss_fracs), curves_non[g], color=cols[g], label=lab[g]); end
    axislegend(ax1; position=:lb, framevisible=false)

    ax2 = Axis(fig[1,2:3], xlabel="Area lost (fraction)", ylabel="ΔBSH",
               title="Full ABM (B ON) — hardened biotic gate")
    for g in geoms; lines!(ax2, collect(loss_fracs), curves_bio[g], color=cols[g], label=lab[g]); end
    axislegend(ax2; position=:lb, framevisible=false)

    display(fig)
    return fig, curves_non, curves_bio
end

function sweep_theta_lambda(; θ_list=0.5:0.05:0.8, λ_list=0.0:0.1:1.0, fstar=0.6,
                            kwargs...)
    grid = make_grid(60,60; seed=42)
    H = fill(NaN, length(λ_list), length(θ_list))
    for (i,λ) in enumerate(λ_list), (j,θ) in enumerate(θ_list)
        fig, non, bio = show_biotic_showcase(; A_level=:intermediate, M_level=:off,
                                             λ_align=λ, θB=θ, diet_cap=2, density=0.08,
                                             loss_fracs=0.2:0.1:0.8, kwargs...)
        idx = findfirst(x->isapprox(x,fstar; atol=1e-6), 0.2:0.1:0.8)
        did_front = bio[:front][idx] - non[:front][idx]
        H[i,j] = did_front
        close(fig)
    end
    fig = Figure(; size=(780,420))
    ax  = Axis(fig[1,1], xlabel="θB (fraction gate)", ylabel="λ_align", title="DiD(front) at f*=$(fstar)")
    hm  = heatmap!(ax, 1:size(H,2), 1:size(H,1), H; colormap=:viridis)
    ax.xticks = (1:length(θ_list), round.(collect(θ_list); digits=2))
    ax.yticks = (1:length(λ_list), round.(collect(λ_list); digits=2))
    Colorbar(fig[1,2], hm; label="ΔBSH")
    display(fig)
    return fig, H
end

fig, non, bio = show_biotic_showcase(; A_level=:intermediate, M_level=:on,
    λ_align=0.9, ξ_align=0.06, diet_cap=2, density=0.08,
    bmode=:frac, θB=0.7, β=4.5, τA=0.48, τocc=0.45,
    front_axis=:x, front_noise=0.0)

