module Biotic

export build_basal_biotic, propagate_biotic

using Random
using ..SpatialUtils: gaussian_blur, norm01

function build_basal_biotic(nb, nx, ny; σ=2.0, seed=1)
    rng = MersenneTwister(seed)
    B = Vector{Matrix{Bool}}(undef, nb)

    for k in 1:nb
        noise = rand(rng, nx, ny)
        smooth = gaussian_blur(noise, σ)
        field = smooth .> 0.5
        B[k] = field
    end
    B
end

function propagate_biotic(metaweb, B_basal)
    A = metaweb.A
    S = size(A,1)
    nx, ny = size(B_basal[1])
    B = [zeros(Bool, nx, ny) for _ in 1:S]

    for (k, idx) in enumerate(metaweb.basal)
        B[idx] = B_basal[k]
    end

    order = sortperm(metaweb.s)

    for i in order
        if i ∈ metaweb.basal
            continue
        end
        prey = findall(A[i, :] .== 1)
        tmp = zeros(Bool, nx, ny)
        for p in prey
            tmp .|= B[p]
        end
        B[i] = tmp
    end

    B
end

end
