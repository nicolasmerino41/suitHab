using CairoMakie

"""
stackedbars_poly!(ax, x, Y; width=0.8, colors=[:dodgerblue,:orange,:forestgreen])

- x :: Vector{<:Real}           positions (one per regime)
- Y :: Matrix{<:Real}           N×K matrix; columns = components (Climate, Interaction, Synergy)
- width :: Real                 bar width
- colors :: Vector{Colorant}    color per component (length K)
"""
function stackedbars_poly!(ax::Axis, x::AbstractVector, Y::AbstractMatrix;
                           width=0.8, colors=[:dodgerblue, :orange, :forestgreen])
    n, k = size(Y)
    @assert length(x) == n "x and Y size mismatch"
    @assert length(colors) == k "colors length must equal number of components (columns of Y)"

    halfw = width/2
    for i in 1:n
        pos = 0.0
        neg = 0.0
        for j in 1:k
            h = Y[i, j]
            x0 = x[i] - halfw
            x1 = x[i] + halfw
            if h ≥ 0
                y0, y1 = pos, pos + h
                poly!(ax, Point2f[(x0, y0), (x1, y0), (x1, y1), (x0, y1)];
                      color = colors[j])
                pos = y1
            else
                y0, y1 = neg + h, neg
                poly!(ax, Point2f[(x0, y0), (x1, y0), (x1, y1), (x0, y1)];
                      color = colors[j])
                neg = y0
            end
        end
    end
    return ax
end
