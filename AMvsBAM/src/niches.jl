module Niches

using Random, Statistics, Distributions

export make_niches

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
