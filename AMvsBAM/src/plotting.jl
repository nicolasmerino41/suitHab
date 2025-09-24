module Plotting

using CairoMakie
const Mke = CairoMakie
export verrorbars!, heatmap_with_colorbar!, ecdf_plot!

# vertical error bars at x with mean y spanning [lo, hi]
function verrorbars!(ax, xs::AbstractVector, y::AbstractVector, lo::AbstractVector, hi::AbstractVector;
                     whisker::Float64=0.03, color=:black, lw::Real=1.5)
    @assert length(xs)==length(y)==length(lo)==length(hi)
    for (x, m, l, h) in zip(xs,y,lo,hi)
        lines!(ax, [x,x], [l,h]; color=color, linewidth=lw)
        lines!(ax, [x - whisker, x + whisker], [l,l]; color=color, linewidth=lw)
        lines!(ax, [x - whisker, x + whisker], [h,h]; color=color, linewidth=lw)
    end
end

function heatmap_with_colorbar!(fig, pos, X, Y, Z; xlabel="", ylabel="", title="")
    ax = Axis(fig[pos...], xlabel=xlabel, ylabel=ylabel, title=title)
    hm = heatmap!(ax, X, Y, Z)
    Colorbar(fig[pos[1], pos[2]+1], hm)
    return ax
end

function ecdf_plot!(ax, data::AbstractVector{<:Real}; label="")
    s = sort(data)
    n = length(s)
    x = s
    y = range(1/n, 1; length=n)
    lines!(ax, x, y; label=label)
end

end # module
