begin
    # assume you already have `grid`; else make a quick one
    if !@isdefined(grid)
        grid = make_grid(60, 60; seed=11, texture=0.10)
    end

    clim = reshape(grid.climate, grid.nx, grid.ny)'  # matrix for plotting

    fig = Figure(; size = (560, 480))
    ax  = Axis(fig[1,1], title = "Grid climate", xlabel = "x", ylabel = "y")
    hm  = heatmap!(ax, clim)
    Colorbar(fig[1,2], hm, label = "climate")
    display(fig)
end
