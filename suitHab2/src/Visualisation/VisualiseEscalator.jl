module VisualiseEscalator

export visualise_escalator

using LinearAlgebra, Measures
using Distributions
using Graphs
using GraphPlot
using NetworkLayout, GraphMakie
using CairoMakie
"""
    visualise_escalator(
        A, s; B,
        xnoise = 0.25,
        figres = (1600,1200)
    )

Plots the PPM network with:
- x wrapping every B species
- exact TL on y-axis
- jitter on x-axis for non-basal nodes
- arrow edges
- color bins by TL (basal green, TL∈[1,2)=blue, [2,3)=orange, [3,4)=red)
"""
function visualise_escalator(A, s; B, xnoise=0.75, figres=(1100,620))

    g = DiGraph(A)
    S = length(s)

    # --- Base x wrapping ---
    x = Float64.([(i-1) % B for i in 1:S])
    y = s

    # --- Add jitter ONLY to non-basals ---
    for i in 1:S
        if s[i] > 1.0
            x[i] += (rand() - 0.5) * xnoise   # small shift
        end
    end

    # --- Color nodes by TL bins ---
    colors = Vector{Symbol}(undef, S)
    for i in 1:S
        if s[i] == 1
            colors[i] = :green
        elseif 1 < s[i] <= 2
            colors[i] = :blue
        elseif 2 <= s[i] < 3
            colors[i] = :orange
        elseif 3 <= s[i] < 4
            colors[i] = :red
        elseif 4 <= s[i] < 5
            colors[i] = :purple
        else
            colors[i] = :brown
        end
    end

    # --- Create figure ---
    fig = Figure(; size=figres)
    ax = Axis(fig[1,1];
        title = "Escalator Trophic Layout",
        xlabel = "X (wrapped every B species)",
        ylabel = "Trophic Level",
    )

    for prey in 1:S
        for pred in inneighbors(g, prey)

            # line from prey to predator
            lines!(ax,
                [x[prey], x[pred]],
                [y[prey], y[pred]];
                color=:black, linewidth=0.5
            )

            # arrowhead at predator
            dx = x[pred] - x[prey]
            dy = y[pred] - y[prey]
            θ = atan(dy, dx)

            arrow_x = x[pred] - 0.08*cos(θ)
            arrow_y = y[pred] - 0.08*sin(θ)

            poly!(ax,
                Point2f[
                    (arrow_x, arrow_y),
                    (arrow_x - 0.05*cos(θ + 0.3),
                    arrow_y - 0.05*sin(θ + 0.3)),
                    (arrow_x - 0.05*cos(θ - 0.3),
                    arrow_y - 0.05*sin(θ - 0.3))
                ],
                color=:black)
        end
    end


    # --- Draw nodes ---
    scatter!(ax, x, y; color=colors, markersize=10)

    # --- Axis limits ---
    xlims!(ax, -1, B)
    ylims!(ax, minimum(y)-0.2, maximum(y)+0.2)

    display(fig)
end

# visualise_escalator(metaweb.A, metaweb.s; B=B)

end