"Keep the warm (:high) or cool (:low) climate tail; :auto picks the tail that feeds more consumers."
function front_climate_mask(rng::AbstractRNG, grid::Grid, keep_frac::Float64;
                            tail::Symbol=:auto, noise::Float64=0.0,
                            prey_tail_share::Union{Nothing,Float64}=nothing)
    C     = grid.C
    clim  = grid.climate
    # if :auto we assume caller computed which tail holds more consumer prey (see diag below)
    chosen = tail
    if tail === :auto
        @assert prey_tail_share !== nothing "front_climate_mask(:auto) requires prey_tail_share"
        chosen = prey_tail_share ≥ 0.5 ? :high : :low
    end
    q = chosen === :high ? quantile(clim, 1 - keep_frac) : quantile(clim, keep_frac)
    keep = BitVector(undef, C)
    if noise ≤ 0
        if chosen === :high
            @inbounds for i in 1:C; keep[i] = clim[i] ≥ q; end
        else
            @inbounds for i in 1:C; keep[i] = clim[i] ≤ q; end
        end
    else
        ϵ = noise .* (rand(rng, C) .- 0.5)
        if chosen === :high
            @inbounds for i in 1:C; keep[i] = clim[i] + ϵ[i] ≥ q; end
        else
            @inbounds for i in 1:C; keep[i] = clim[i] + ϵ[i] ≤ q; end
        end
    end
    return keep
end

function metaweb_diagnostics(pool::SpeciesPool)
    cons = findall(!, pool.basal)
    L    = mean(length.(pool.E[cons]))                    # average number of prey
    μprey = map(s-> isempty(pool.E[s]) ? NaN : mean(pool.mu[pool.E[s]]), cons)
    ρ    = cor(filter(!isnan, pool.mu[cons]), filter(!isnan, (μprey)))   # climate assortativity
    return (; mean_links=L, assortativity=ρ)
end

"Share of each consumer's prey that sit in the warm tail (top 30% climate)."
function warm_tail_share_for_consumers(pool::SpeciesPool, grid::Grid; top_frac::Float64=0.30)
    C = grid.C
    q = quantile(grid.climate, 1 - top_frac)
    warm = findall(>=(q), grid.climate)
    cons = findall(!, pool.basal)
    # fraction of each consumer's prey that *can* be suitable in warm cells (by niche center)
    shares = Float64[]
    for s in cons
        prey = pool.E[s]
        isempty(prey) && (push!(shares, 0.0); continue)
        push!(shares, mean(pool.mu[prey] .>= 0.5))  # coarse: prey biased to warm niche centers
    end
    # return mean share across consumers (rough proxy)
    return mean(shares)
end

# # after building pool:
# diag = metaweb_diagnostics(pool)
# println("links≈", round(diag.mean_links, digits=2), "  assort≈", round(diag.assortativity, digits=2))

# warm_share = warm_tail_share_for_consumers(pool, grid)  # ≳0.6 means consumers lean warm
# println("consumer prey warm-tail share ≈ ", round(warm_share, digits=2))

# # when making the front mask at keepfrac = 1 - f
# keep = hl_kind === :front_clim ?
#     front_climate_mask(rng_mask, grid, keepfrac; tail=:auto, prey_tail_share=warm_share) :
#         frontlike_mask(rng_mask, grid, keepfrac; axis=front_axis, noise=front_noise) :
#             error("Unknown hl_kind=$hl_kind")

function biotic_showcase(; nx=60, ny=60, S=140, basal_frac=0.45,
    A_level=:divergent,                     # keeps a strong climate axis
    density=0.055, diet_cap=2,              # low redundancy
    β=6.0, τA=0.50, τocc=0.50,              # stricter suitability & occupancy
    bmode=:frac, θB=0.75, kB=2,             # hard biotic gate
    M_level=:off,                           # no movement buffering
    loss_fracs=0.2:0.1:0.8, seeds_mask=1:6, sim_seed=1234)

    grid = make_grid(nx,ny; seed=42)
    rngp = MersenneTwister(hash((sim_seed,:pool,1)))
    # aligned metaweb with small niche jitter and diet cap
    pool = build_pool_aligned(rngp; S, basal_frac, A_level,
                              density=density, pmax=0.9,
                              λ_align=1.0, ξ_align=0.04, diet_cap=diet_cap)
    println(metaweb_diagnostics(pool))
    A    = abiotic_matrix(pool, grid)
    base = trues(grid.C)

    warm_share = warm_tail_share_for_consumers(pool, grid)   # tells front_clim which tail to cut

    function meanΔBSH(B_on::Bool, hl_kind::Symbol)
        pars = bam_from_axes(; B_level=(B_on ? :strong : :none), M_level, τA, τocc)
        bam  = BAMParams(; α=1.0, β=β, μ=pars.bam.μ, γ=pars.bam.γ, τA=τA, τocc=τocc)
        mp   = MovementParams(; mode=pars.mp.mode, T=8)

        vals = Float64[]
        for f in loss_fracs
            keepfrac = 1 - f
            Δ = Float64[]
            for ms in seeds_mask
                rngm = MersenneTwister(hash((sim_seed,:mask,ms,hl_kind)))
                keep = hl_kind === :random     ? random_mask(rngm, grid.C, keepfrac) :
                       hl_kind === :clustered  ? clustered_mask(rngm, grid, keepfrac; nseeds=6) :
                       hl_kind === :front_clim ? front_climate_mask(rngm, grid, keepfrac; tail=:auto,
                                                                    prey_tail_share=warm_share) :
                       error("hl_kind")
                P0,B0,_ = assemble_BAM_hard(pool, grid, A, base; bam=bam, mp=mp,
                                            bmode=bmode, θB=θB, kB=kB)
                S0 = species_stats(pool, grid, A, base, P0, B0, ones(Float64,size(A)), bam)
                P1,B1,_ = assemble_BAM_hard(pool, grid, A, keep; bam=bam, mp=mp,
                                            bmode=bmode, θB=θB, kB=kB)
                S1 = species_stats(pool, grid, A, keep, P1, B1, ones(Float64,size(A)), bam)
                push!(Δ, mean(S1.BSH[.!pool.basal] .- S0.BSH[.!pool.basal]))
            end
            push!(vals, mean(Δ))
        end
        return vals
    end

    geoms = (:random, :clustered, :front_clim)
    Δnon  = Dict(g=>meanΔBSH(false,g) for g in geoms)
    Δbio  = Dict(g=>meanΔBSH(true, g) for g in geoms)

    # quick plot (Makie)
    fig = Figure(; size=(1380, 500))
    Label(fig[0,1:2], "HL effect with vs without biotic — A=$(A_level), M=$(M_level)"; fontsize=18)
    lab  = Dict(:random=>"Random", :clustered=>"Clustered", :front_clim=>"Front (climate-tail)")
    col  = Dict(:random=>:dodgerblue3, :clustered=>:darkorange3, :front_clim=>:seagreen4)

    ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="ΔBSH (mean cons.)", title="Climate-only (B OFF)")
    for g in geoms; lines!(ax1, collect(loss_fracs), Δnon[g], color=col[g], label=lab[g]); end
    axislegend(ax1; position=:lb, framevisible=false)

    ax2 = Axis(fig[1,2], xlabel="Area lost (fraction)", ylabel="ΔBSH",
               title="Full ABM (B ON) — diet cap 2 + θB=0.75")
    for g in geoms; lines!(ax2, collect(loss_fracs), Δbio[g], color=col[g], label=lab[g]); end
    axislegend(ax2; position=:lb, framevisible=false)

    display(fig)
    return (; fig, Δnon, Δbio, warm_share)
end

biotic_showcase(; A_level=:divergent, M_level=:off,
    density=0.055, diet_cap=2, β=6.0, τA=0.50, τocc=0.50,
    bmode=:frac, θB=0.75, kB=2)
