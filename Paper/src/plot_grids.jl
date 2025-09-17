# plot_grids.jl
# Visualize several grid types (climate maps) side-by-side.
using Random, Statistics, Printf
using CairoMakie

function plot_grid!(ax, grid::Grid; title::String="")
    nx, ny = grid.nx, grid.ny
    # reshape climate to (ny,nx) for image()
    Z = reshape(grid.climate, (nx, ny))'  # rows=y, cols=x
    image!(ax, 1:nx, 1:ny, Z; interpolate=false, colormap=:viridis)
    hidespines!(ax); hidedecorations!(ax, grid=false)
    ax.title = title
end

function plot_all_grids(; nx=60, ny=60, out="out/grids.png")
    gridA = make_grid_gradient(nx,ny; seed=42, texture=0.10)
    gridB = make_grid_patchy(nx,ny; seed=43, octaves=3, amp=0.5)
    gridC = make_grid_ridge(nx,ny; seed=44)

    fig = Figure(resolution=(1200, 380))
    plot_grid!(Axis(fig[1,1]), gridA; title="Gradient (smooth+wiggle)")
    plot_grid!(Axis(fig[1,2]), gridB; title="Patchy (Gaussian bumps)")
    plot_grid!(Axis(fig[1,3]), gridC; title="Ridge (oblique band)")
    Colorbar(fig[1,4], limits=(0,1), colormap=:viridis, label="climate")

    # save(out, fig); 
    display(fig)
end
