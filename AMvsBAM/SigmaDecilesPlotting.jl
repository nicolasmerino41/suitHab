# ---------------------------
# ΔArea vs R95 conditioned on σ deciles
# ---------------------------
using Statistics, DataFrames, CairoMakie
const Mke = CairoMakie

function plot_deltaarea_vs_R95_by_sigma_quantiles(
    res2::DataFrame; 
    nquant::Int=10, savepath::Union{Nothing,String}=nothing
)
    # @assert :R95 ∈ names(res2) && :σ ∈ names(res2) && :ΔAmean ∈ names(res2)

    # 1) Build σ quantile bins (deciles by default)
    qs     = range(0, 1; length=nquant+1) |> collect
    edges  = quantile(res2.σ, qs)                         # length = nquant+1
    # guard against identical edges (can happen if σ grid is coarse)
    for i in 2:length(edges)
        if edges[i] ≤ edges[i-1]
            edges[i] = nextfloat(edges[i-1])
        end
    end

    # Assign each row to its σ-decile (1..nquant)
    decile = similar(res2.σ, Int)
    for (i, s) in pairs(res2.σ)
        # searchsortedlast gives index in 1..nquant+1; clamp to 1..nquant
        d = searchsortedlast(edges, s)
        decile[i] = clamp(d, 1, nquant)
    end
    res2_with = copy(res2)
    res2_with.decile = decile

    # 2) Aggregate: for each (R95, decile) compute mean ΔAmean
    g = combine(groupby(res2_with, [:R95, :decile]), :ΔAmean => mean => :ΔAmean_mean)

    # For labeling/legend: compute σ-range per decile
    σlo = [minimum(res2_with.σ[res2_with.decile .== d]) for d in 1:nquant]
    σhi = [maximum(res2_with.σ[res2_with.decile .== d]) for d in 1:nquant]
    σlab = ["σ∈[$(@sprintf("%.3f", σlo[d])),$(@sprintf("%.3f", σhi[d]))]" for d in 1:nquant]

    # 3) Plot: one line per decile (all deciles; color by decile)
    fig = Mke.Figure(; size=(950, 500))
    ax  = Mke.Axis(fig[1,1];
        xlabel = "R95 (diet redundancy)",
        ylabel = "ΔArea (AM − BAM)",
        title  = "ΔArea vs R95 conditioned on σ deciles"
    )

    # consistent R95 ordering on the x-axis
    xs = sort(unique(g.R95))

    cmap = cgrad(:viridis, nquant)
    for d in 1:nquant
        sub = g[g.decile .== d, :]
        # ensure y follows xs order
        y = [first(sub.ΔAmean_mean[(sub.R95 .== x)]) for x in xs]
        Mke.lines!(ax, xs, y; color=cmap[d], label=σlab[d])
    end
    Mke.axislegend(ax; position=:rb, framevisible=false, nbanks=2)
    display(fig)

    if savepath !== nothing
        save(savepath, fig)
    end
    return fig
end

# Example usage (saves PNG next to your other figs)
mkpath(joinpath(@__DIR__, "..", "data", "figs"))
plot_deltaarea_vs_R95_by_sigma_quantiles(
    res2;
    nquant=10,
    savepath=joinpath(@__DIR__, "data", "figs", "Lines_DeltaArea_vs_R95_bySigmaDeciles.png")
)
