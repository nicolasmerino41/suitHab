module MetawebDescriptive

using CairoMakie, Statistics, Printf
using ..Metawebs

export plot_metaweb_spectrum

# -- helpers -------------------------------------------------------------------

"Integer diet sizes for consumers only."
_consumer_diets(pool::Metawebs.SpeciesPool) =
    [length(pool.prey[s]) for s in 1:pool.S if !pool.basal[s]]

"Adjacency matrix A[s,q] = 1 if predator s eats prey q (entire S×S)."
function _adjacency_matrix(pool::Metawebs.SpeciesPool)
    S = pool.S
    A = zeros(Int, S, S)
    @inbounds for s in 1:S
        for q in pool.prey[s]
            (1 ≤ q ≤ S) && (A[s,q] = 1)
        end
    end
    A
end

"Simple bar counts for integer k."
function _int_hist_counts(ks::Vector{Int})
    isempty(ks) && return (Int[], Int[])
    kmin, kmax = minimum(ks), maximum(ks)
    xs = collect(kmin:kmax)
    ys = [count(==(x), ks) for x in xs]
    xs, ys
end

"One figure: adjacency (left), diet histogram (middle), mass curve (right)."
function _plot_one_metaweb(pool::Metawebs.SpeciesPool;
        name::AbstractString="metaweb", outpath::AbstractString="figs/metawebs/metaweb.png")

    S = pool.S
    A = _adjacency_matrix(pool)
    diets = _consumer_diets(pool)
    masses_sorted = sort(pool.masses)

    # diagnostics from your own helper
    d = Metawebs.metaweb_diagnostics(pool)

    # ---- layout
    fig = Figure(; size=(1900, 650))

    # 1) adjacency (pred rows, prey cols)
    ax1 = Axis(fig[1,1]; title="$name — adjacency (pred rows, prey cols)")
    # draw as scattered pixels for speed
    I = Int[]; J = Int[]
    @inbounds for i in 1:S, j in 1:S
        if A[i,j] == 1
            push!(J, j); push!(I, i)
        end
    end
    scatter!(ax1, J, I; markersize=2.2, color=:black)
    xlims!(ax1, 0.5, S+0.5); ylims!(ax1, 0.5, S+0.5)
    # yflip!(ax1, true)
    hidespines!(ax1, :t, :r)

    # 2) diet histogram (consumers)
    ax2 = Axis(fig[1,2]; title="Diet sizes (consumers)", xlabel="k", ylabel="count")
    xs, ys = _int_hist_counts(Int.(diets))
    if !isempty(xs)
        barplot!(ax2, xs, ys)
        txt = @sprintf("mean=%.2f   med=%d   conn=%.3f", mean(diets), Int(median(diets)), d.C)
        # put label near top-left of plot area
        xmax = maximum(xs); ymax = maximum(ys)
        text!(ax2, minimum(xs) + 0.02*(xmax - minimum(xs)),
                   ymax - 0.05*(ymax>0 ? ymax : 1.0);
                   text=txt, align=(:left,:top), color=:dimgray)
    end
    hidespines!(ax2, :t, :r)

    # 3) mass order
    ax3 = Axis(fig[1,3]; title="Mass order", xlabel="rank", ylabel="mass", yscale=log10)
    lines!(ax3, 1:S, masses_sorted; color=:gray25, linewidth=3)
    hidespines!(ax3, :t, :r)

    mkpath(dirname(outpath))
    save(outpath, fig)
    fig
end

# -- public --------------------------------------------------------------------

"""
plot_metaweb_spectrum(pools; names, outdir="figs/metawebs")

Create one descriptive figure per pool:
- left: adjacency raster (pred rows × prey cols)
- middle: histogram of consumer diet sizes (with mean / median / connectance)
- right: body-mass vs rank (log-y)

`pools` and `names` must have the same length.
"""
function plot_metaweb_spectrum(pools::Vector{Metawebs.SpeciesPool};
        names::Vector{<:AbstractString},
        outdir::AbstractString="figs/metawebs")

    @assert length(pools) == length(names) "pools and names must match in length"
    mkpath(outdir)
    figs = Any[]
    for (pool, nm) in zip(pools, names)
        fig = _plot_one_metaweb(pool; name=string(nm),
                                outpath=joinpath(outdir, "metaweb_desc_"*string(nm)*".png"))
        push!(figs, fig)
        display(fig)   # ✅ show it immediately
    end
    figs
end

end # module
