module Niches

export make_niches

using Random, Distributions

function make_niches(rng::AbstractRNG, S::Int; σ=0.12)
    μ = rand(rng, S)
    σi = fill(σ, S)
    μ, σi
end

end
