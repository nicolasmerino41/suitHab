module Niches

using Random, Statistics, Distributions

export make_niches

# """
#     make_climate_grid(nx, ny; kind=:gradient)

# Return an `nx × ny` grid of a 1-D climate variable in [0,1].
# `kind` can be `:gradient` or `:ridge` (bimodal/patchy).
# """
# function make_climate_grid(nx::Int, ny::Int; kind::Symbol=:gradient, seed::Int=1)
#     rng = MersenneTwister(seed)
#     if kind == :gradient
#         x = range(0, 1; length=nx)
#         y = range(0, 1; length=ny)
#         C = [xi for xi in x, yi in y]  # simple x-gradient
#     else
#         # ridge: two bands
#         C = zeros(nx, ny)
#         for i in 1:nx, j in 1:ny
#             C[i,j] = rand(rng) < 0.5 ? 0.25 + 0.05*randn(rng) : 0.75 + 0.05*randn(rng)
#         end
#         C .= clamp.(C, 0, 1)
#     end
#     return C
# end

"""
    make_niches(rng, S; align=0.5, σ=0.12, basal_frac=0.3)

Return species climatic optima `μ` in [0,1] and niche widths `σi`.
`align` controls similarity between consumer and prey optima:
- align=1: consumers share prey climates
- align=0: independent draws
We simply implement via adding correlated noise around prey pool mean.
"""
function make_niches(rng::AbstractRNG, S::Int; align::Float64=0.5, σ::Float64=0.12, basal_frac::Float64=0.3)
    nb = max(1, round(Int, basal_frac*S))
    μ = rand(rng, S)  # start random
    # force consumers to be correlated with lower-index (prey-dominated) set
    baseμ = μ[1:nb]
    for j in nb+1:S
        μ[j] = clamp(align*mean(baseμ) + (1-align)*rand(rng), 0, 1)
    end
    σi = fill(σ, S)
    return μ, σi
end

end # module
