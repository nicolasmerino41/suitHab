# Re-run a small sweep across R95 at one σ to extract composition effects
σ_fixed = 0.14
sweepR  = NamedTuple[(; R95=r) for r in 1:10]
fixedR  = (; Cgrid=Cgrid, align=default_align, σ=σ_fixed, R95=default_R95,
           C=default_C, S=S, basal_frac=basal_frac, τA=τA, kreq=kreq)

function once_with_breakdown(rng, pars, fixed)
    v = run_once(MersenneTwister(rand(rng,1:10^9)); Cgrid = fixed.Cgrid,
                 align=fixed.align, σ=fixed.σ,
                 R95=get(pars,:R95,fixed.R95), connectance=fixed.C,
                 motif_mix=get(pars,:motif,:mixed),
                 S=fixed.S, basal_frac=fixed.basal_frac, τA=fixed.τA, kreq=fixed.kreq)
    # We need species-level ΔA to split by basal/consumer
    mw = MetaWeb.build_metaweb(MersenneTwister(42); S=fixed.S, basal_frac=fixed.basal_frac,
                               connectance=fixed.C, R95=get(pars,:R95,fixed.R95), motif_mix=:mixed)
    # recompute AM/BAM maps deterministically for the same R95 (cheap once here)
    μ, σi = Niches.make_niches(MersenneTwister(43), fixed.S; align=fixed.align, σ=fixed.σ, basal_frac=fixed.basal_frac)
    out = BAM.compute_AM_BAM(MersenneTwister(44), mw, fixed.Cgrid, μ, σi, BAM.Params(;τA=fixed.τA, kreq=fixed.kreq))
    am = BAM.species_areas(out[:AM_maps]); bm = BAM.species_areas(out[:BAM_maps])
    is_cons = mw.trophic_role .!= :basal
    Δ_cons = mean(am[is_cons] .- bm[is_cons])
    Δ_basal = mean(am[.!is_cons] .- bm[.!is_cons])  # ~0
    frac_cons = mean(is_cons)
    return (; R95=get(pars,:R95,fixed.R95), Δ_all=v.ΔA, Δ_cons, Δ_basal, frac_cons, P=v.Psuff, πA=v.πA)
end

Rres = DataFrame([once_with_breakdown(rng, p, fixedR) for p in sweepR])

begin
    fig = Figure(; size=(1100,420))
    ax1 = Axis(fig[1,1], xlabel="R95", ylabel="ΔArea", title="ΔArea: all vs consumers")
    lines!(ax1, Rres.R95, Rres.Δ_all; label="mean over all species")
    lines!(ax1, Rres.R95, Rres.Δ_cons; linestyle=:dash, label="consumers only")
    axislegend(ax1, position=:rb, framevisible=false)

    ax2 = Axis(fig[1,2], xlabel="R95", ylabel="fraction consumers", title="Composition shifts")
    lines!(ax2, Rres.R95, Rres.frac_cons)

    ax3 = Axis(fig[1,3], xlabel="R95", ylabel="", title="Mechanism checks")
    lines!(ax3, Rres.R95, Rres.P; label="P_suff")
    lines!(ax3, Rres.R95, Rres.πA; label="π_A")
    axislegend(ax3, position=:lb, framevisible=false)

    display(fig)
    save(joinpath(@__DIR__, "data", "figs", "R95_decomposition.png"), fig)
end
