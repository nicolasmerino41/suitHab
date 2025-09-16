function biotic_showcase_frontclim(; nx=60, ny=60, S=140, basal_frac=0.45,
    A_level=:divergent,                           # or :intermediate
    λ_align=1.0, ξ_align=0.04, diet_cap=2, density=0.06,
    β=5.0, τA=0.48, τocc=0.46, M_level=:off,
    loss_fracs=0.2:0.1:0.8, seeds_mask=1:6, sim_seed=1234,
    bmode=:frac, θB=0.70, kB=2)

    grid = make_grid(nx,ny; seed=42)
    rngp = MersenneTwister(hash((sim_seed,:pool,1)))

    # metaweb: low redundancy + strong climate assortativity
    pool = build_pool_aligned(rngp; S, basal_frac, A_level,
                              density=density, λ_align=λ_align, ξ_align=ξ_align,
                              diet_cap=diet_cap)
    @show metaweb_diagnostics(pool)

    A    = abiotic_matrix(pool, grid)
    base = trues(grid.C)

    # helper for one geometry
    function mean_curve(B_on::Bool, hl_kind::Symbol)
        pars = bam_from_axes(; B_level = (B_on ? :strong : :none), M_level, τA, τocc)
        bam  = BAMParams(; α=1.0, β=β, μ=pars.bam.μ, γ=pars.bam.γ, τA=τA, τocc=τocc)
        mp   = MovementParams(; mode=pars.mp.mode, T=8)

        vals = Float64[]
        for f in loss_fracs
            keepfrac = 1 - f
            Δs = Float64[]
            for ms in seeds_mask
                rngm = MersenneTwister(hash((sim_seed,:mask,ms,hl_kind)))
                keep = hl_kind === :random     ? random_mask(rngm, grid.C, keepfrac) :
                       hl_kind === :clustered  ? clustered_mask(rngm, grid, keepfrac; nseeds=6) :
                       hl_kind === :front_clim ? front_climate_mask(rngm, grid, keepfrac; tail=:high, noise=0.0) :
                       error("hl_kind must be :random/:clustered/:front_clim")

                P0,B0,_ = assemble_BAM_hard(pool, grid, A, base; bam=bam, mp=mp, bmode=bmode, θB=θB, kB=kB)
                S0 = species_stats(pool, grid, A, base, P0, B0, ones(Float64, size(A)), bam)
                P1,B1,_ = assemble_BAM_hard(pool, grid, A, keep; bam=bam, mp=mp, bmode=bmode, θB=θB, kB=kB)
                S1 = species_stats(pool, grid, A, keep, P1, B1, ones(Float64, size(A)), bam)
                push!(Δs, mean(S1.BSH[.!pool.basal] .- S0.BSH[.!pool.basal]))
            end
            push!(vals, mean(Δs))
        end
        return vals
    end

    geoms = (:random, :clustered, :front_clim)
    curves_non = Dict(g=>mean_curve(false, g) for g in geoms)
    curves_bio = Dict(g=>mean_curve(true,  g) for g in geoms)

    # --- plot
    fig = Figure(; size=(1400, 480))
    Label(fig[0,1:3], "HL effect with vs without biotic — A=$(A_level), M=$(M_level)"; fontsize=18)

    lab  = Dict(:random=>"Random", :clustered=>"Clustered", :front_clim=>"Front (climate-tail)")
    cols = Dict(:random=>:dodgerblue3, :clustered=>:darkorange3, :front_clim=>:seagreen4)

    ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="ΔBSH (mean consumers)",
               title="Climate-only (B OFF)")
    for g in geoms; lines!(ax1, collect(loss_fracs), curves_non[g], color=cols[g], label=lab[g]); end
    axislegend(ax1; position=:lb, framevisible=false)

    ax2 = Axis(fig[1,2:3], xlabel="Area lost (fraction)", ylabel="ΔBSH",
               title="Full ABM (B ON) — diet cap 2 + threshold gate")
    for g in geoms; lines!(ax2, collect(loss_fracs), curves_bio[g], color=cols[g], label=lab[g]); end
    axislegend(ax2; position=:lb, framevisible=false)

    display(fig)
    return fig, curves_non, curves_bio
end
function metaweb_diagnostics(pool::SpeciesPool)
    cons = findall(!, pool.basal)
    L    = mean(length.(pool.E[cons]))
    μprey = map(s-> isempty(pool.E[s]) ? NaN : mean(pool.mu[pool.E[s]]), cons)
    ρ    = cor(filter(!isnan, pool.mu[cons]), filter(!isnan, (μprey)))  # climate assortativity
    return (; mean_links=L, assortativity=ρ)
end
# Expect: mean_links ≈ 2, assortativity ρ ≳ 0.6

biotic_showcase_frontclim(; A_level=:divergent, M_level=:on,
    λ_align=1.0, ξ_align=0.04, diet_cap=2, density=0.06,
    bmode=:frac, θB=0.70, β=5.0, τA=0.48, τocc=0.46)
