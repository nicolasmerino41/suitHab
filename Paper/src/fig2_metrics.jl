module Fig2Metrics

using CairoMakie, Statistics, Printf
using ..BSH
using ..HL
using ..Metrics

export fig2_metrics_block, per_species_failure_cdf, penalty_vs_diet

"""
Make a compact table: rows = {KS, tailΔ, GiniΔ}, cols = {random, clustered, front}.
- tailΔ = P_BAM(loss>tail_cut) - P_AM(loss>tail_cut)
- GiniΔ  = Gini_BAM - Gini_AM
Saves figure to `out`.
"""
function fig2_metrics_block(; rng, pool, grid, pars, fstar::Float64,
                             tail_cut::Float64=0.8, seed_A::Int=1,
                             out::AbstractString="figs/fig2_metrics/Fig2_metrics.png")
    geoms = (:random, :clustered, :front)
    stats = Dict{Symbol,NamedTuple}()

    for g in geoms
        rAM, rBAM = BSH.per_species_relative_loss(rng, pool, grid, pars;
                                                  fstar=fstar, geometry=g, seed_A=seed_A)
        stats[g] = Metrics.dist_metrics(rAM, rBAM; tail_cut=tail_cut)
    end

    # unwrap: take Δ (last element) if tuple, otherwise the value
    unwrap(x) = x isa Tuple ? x[end] : x

    fig = Figure(; size=(900, 330))
    ax = Axis(
        fig[1,1],
        xticks=(1:3, collect(String.(geoms))),  # now both parts are vectors
        yticks=(1:3, ["KS", "tailΔ", "GiniΔ"]),
        title=@sprintf("Dist. metrics at f* = %.2f (tail=%.2f)", fstar, tail_cut)
    )



    hidedecorations!(ax, grid=false); hidespines!(ax)

    vals = [
        stats[:random].KS       stats[:clustered].KS       stats[:front].KS;
        stats[:random].tailDiff stats[:clustered].tailDiff stats[:front].tailDiff;
        stats[:random].giniDiff stats[:clustered].giniDiff stats[:front].giniDiff
    ]

    for i in 1:3, j in 1:3
        text!(ax, j, i; text=@sprintf("%.3f", vals[i,j]),
            align=(:center,:center), fontsize=22)
    end

    save(out, fig); display(fig)
    # return (; stats, fig)
end

"""
Per-species distribution of biotic gate failure among A-suitable kept cells
at f* (one geometry). Saves a one-panel CDF.
"""
function per_species_failure_cdf(; rng, pool, grid, pars, fstar::Float64,
                                 geometry::Symbol=:front, seed_A::Int=1,
                                 out::AbstractString="figs/Fig2_failure_cdf.png")
    A = BSH.abiotic_matrix(pool, grid; seed=seed_A)
    keepfrac = 1 - fstar
    keep = geometry===:random    ? HL.random_mask(rng, grid.C, keepfrac) :
           geometry===:clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                   HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)

    bam = BSH.assemble_BAM(pool, grid, A, keep; pars=pars)

    cons = [s for s in 1:pool.S if !pool.basal[s]]
    fracfail = Float64[]
    for s in cons
        den = 0; num = 0
        @inbounds for i in 1:grid.C
            if keep[i] & (A[s,i] ≥ pars.τA)
                den += 1
                (bam.preyfrac[s,i] < pars.τB) && (num += 1)
            end
        end
        den>0 && push!(fracfail, num/den)
    end

    ss = sort(fracfail); n = length(ss); yy = (1:n)./n
    fig = Figure(; size=(650, 430))
    ax  = Axis(fig[1,1], xlabel="P(failure | A-suitable & kept)",
                         ylabel="fraction of consumers",
                         title="Per-species biotic failure — $(String(geometry)) at f*=$(fstar)")
    stairs!(ax, ss, yy)
    save(out, fig); display(fig)
    return (; fracfail, fig)
end

"""
Species-level 'BAM penalty' vs diet size at f* (one geometry).
Penalty = (relBAM - relAM) per species (consumers only).
Also prints Spearman ρ.
"""
function penalty_vs_diet(; rng, pool, grid, pars, fstar::Float64,
                         geometry::Symbol=:front, seed_A::Int=1,
                         out::AbstractString="figs/penalty_vs_diet.png")
    relAM, relBAM = BSH.per_species_relative_loss(rng, pool, grid, pars;
                                                  fstar=fstar, geometry=geometry, seed_A=seed_A)
    cons = [s for s in 1:pool.S if !pool.basal[s]]
    pen  = relBAM[cons] .- relAM[cons]
    diets = [length(pool.prey[s]) for s in cons]

    # Spearman by rank
    r_pen   = sortperm(pen);  r_diet = sortperm(diets)
    ρ = cor(r_pen, r_diet)     # Spearman (rank–rank Pearson)

    fig = Figure(; size=(650, 430))
    ax  = Axis(fig[1,1], xlabel="diet size (k)", ylabel="BAM−AM (per species)",
               title="Penalty vs diet — $(String(geometry)) at f*=$(fstar)\nSpearman ρ=$(round(ρ,digits=3))")
    scatter!(ax, diets, pen, markersize=6)
    save(out, fig); display(fig)
    return (; ρ, fig)
end

end # module
