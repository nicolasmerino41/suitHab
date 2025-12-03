module Plotting

export plot_grid

using CairoMakie

function plot_grid(G; title="grid")
    fig = Figure(; size = (1100, 600))
    ax = Axis(fig[1,1]; title=title)
    heatmap!(ax, G; colormap=:viridis, colorrange=(0, 1))
    display(fig)
end

end
