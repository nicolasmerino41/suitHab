module TernaryPlots
using CairoMakie

export ternary_coords, ternary_scatter!

"Convert (wA,wB,wM), sum=1, to 2D ternary coordinates."
ternary_coords(wA, wB, wM) = (x = 0.5*(2wM + wB), y = (√3/2)*wB)

"Plot colored points on an existing Makie axis as a ternary triangle."
function ternary_scatter!(ax, W::Vector{NTuple{3,Float64}}, z::AbstractVector;
        colormap=:viridis, markersize=12)

    # triangle frame A(0,0), M(1,0), B(0.5,√3/2)
    lines!(ax, [0.0, 1.0], [0.0, 0.0], color=:black)
    lines!(ax, [1.0, 0.5], [0.0, √3/2], color=:black)
    lines!(ax, [0.5, 0.0], [√3/2, 0.0], color=:black)
    text!(ax, "A", position=(0.0,-0.03), align=(:left,:top))
    text!(ax, "M", position=(1.0,-0.03), align=(:right,:top))
    text!(ax, "B", position=(0.5, √3/2+0.03), align=(:center,:bottom))

    xs = Float64[]; ys = Float64[]
    for (a,b,m) in W
        c = ternary_coords(a,b,m)
        push!(xs, c.x); push!(ys, c.y)
    end
    sc = scatter!(ax, xs, ys; color=z, colormap=colormap, markersize=markersize)
    return sc
end

end # module
