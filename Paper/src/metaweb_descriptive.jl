# src/metaweb_descriptive.jl
# Descriptive plots for a metaweb: adjacency + diet-size distribution.
# Uses only Metawebs.SpeciesPool fields (masses, basal, prey).

module MetawebDescriptive

using Statistics, Printf
using CairoMakie
using ..Metawebs

export plot_metaweb_descriptive, plot_metaweb_spectrum

# adjacency as Bool matrix (pred rows, prey cols)
function _adjacency(pool::Metawebs.SpeciesPool)
    S = pool.S
    A = falses(S, S)
    @inbounds for s in 1:S
        pool.basal[s] && continue
        for q in pool.prey[s]
            A[s,q] = true
        end
    end
    A
end

diet_sizes(pool::Metawebs.SpeciesPool) = [length(pool.prey[s]) for s in 1:pool.S if !pool.basal[s]]

function connectance(pool::Metawebs.SpeciesPool)
    cons = findall(!, pool.basal)
    E = sum(length.(pool.prey[cons]))
    denom = length(cons) * pool.S
    denom == 0 ? 0.0 : E/denom
end

"""
    plot_metaweb_descriptive(pool; name="metaweb", outdir="Paper/figs/metawebs")

Saves `<outdir>/<name>_desc.png` with:
- adjacency matrix (species ordered by mass),
- diet-size histogram for consumers,
- mass rank line for sanity.
Returns a NamedTuple with (sizes, conn, fig, out).
"""
function plot_metaweb_descriptive(pool::Metawebs.SpeciesPool; name::String="metaweb", outdir::String="Paper/figs/metawebs")
    isdir(outdir) || mkpath(outdir)

    ord = sortperm(pool.masses)   # nicer block structure
    A = _adjacency(pool)[ord, ord]
    sizes = diet_sizes(pool)
    conn  = connectance(pool)

    fig = Figure(;size=(1100,380))

    ax1 = Axis(fig[1,1], title="$name — adjacency (pred rows, prey cols)")
    image!(ax1, Float64.(A)'; interpolate=false, colormap=[:white, :black])
    hidespines!(ax1); hidedecorations!(ax1, grid=false)

    ax2 = Axis(fig[1,2], title="Diet sizes (consumers)", xlabel="k", ylabel="count")
    if isempty(sizes)
        text!(ax2, 0.5, 0.5, text="no consumers", align=(:center,:center))
    else
        kmax = maximum(sizes)
        hist!(ax2, sizes; bins=0:(kmax+1), color=:steelblue)
        text!(ax2, 0.02, 0.92, text=@sprintf("mean=%.2f  med=%d  conn=%.3f",
                                             mean(sizes), median(sizes), conn),
              space=:relative, align=(:left,:top))
    end

    ax3 = Axis(fig[1,3], title="Mass order", xlabel="rank", ylabel="mass")
    lines!(ax3, 1:pool.S, sort(pool.masses); color=:gray)

    out = joinpath(outdir, "$(name)_desc.png")
    save(out, fig)
    display(fig)
    return (; sizes, conn, fig, out)
end

"""
    plot_metaweb_spectrum(pools; names=[], outdir="Paper/figs/metawebs")

Convenience to describe several metawebs; `pools` is a Vector{SpeciesPool}.
"""
function plot_metaweb_spectrum(pools::Vector{Metawebs.SpeciesPool}; names=String[], outdir::String="Paper/figs/metawebs")
    results = NamedTuple[]
    for (i, pool) in enumerate(pools)
        nm = isempty(names) ? "pool$(i)" : names[i]
        push!(results, plot_metaweb_descriptive(pool; name=nm, outdir=outdir))
    end
    results
end

end # module
