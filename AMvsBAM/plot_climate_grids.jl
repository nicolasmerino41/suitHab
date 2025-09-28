function plot_climate_grids(grid_grad, grid_ridge, grid_fractal; outpath=nothing, clim=nothing, titles=nothing)
    # Adapter: accept matrices or structs with values/nx/ny
    to_matrix(g) = g isa AbstractMatrix ? g :
                   reshape(getfield(g, :values), getfield(g, :nx), getfield(g, :ny))

    Gs = (to_matrix(grid_grad), to_matrix(grid_ridge), to_matrix(grid_fractal))

    # Determine color limits if not provided
    if clim === nothing
        gmin = minimum(minimum, Gs)
        gmax = maximum(maximum, Gs)
        clim = (float(gmin), float(gmax))
    end

    # Titles
    titles === nothing && (titles = ["Gradient","Ridge","Fractal"])

    fig = Figure(; size=(1200, 420))

    for (j, (g, title)) in enumerate(zip(Gs, titles))
        ax = Axis(fig[1, j], title=title, xlabel="x", ylabel="y")
        heatmap!(ax, g; colormap=:viridis, colorrange=clim)
    end

    Colorbar(fig[:, end + 1], label="Climate value", colorrange=clim)
    fig[0, :] = Label(fig, "Areas of study (climate)", fontsize=22, tellwidth=false)

    if outpath !== nothing
        save(outpath, fig)
    end
    display(fig)
    return fig
end

# Example (note the variable name fix):
grid_gradient = Climate.make_climate_grid(nx, ny; kind=:gradient, seed=11)
grid_ridge = Climate.make_climate_grid(nx, ny; kind=:ridge, seed=11)
grid_fractal = Climate.make_climate_grid(nx, ny; kind=:fractal, seed=11)
plot_climate_grids(grid_gradient, grid_ridge, grid_fractal; outpath="AMvsBAM/data/figs/final/climate_grids.png")

