module Suitability

export abiotic_suitability, biotic_suitability, combine_suitability, movement_filter

using ..Movement: apply_movement

function abiotic_suitability(climate, μ, σi)
    nx, ny = size(climate)
    S = length(μ)
    A = [zeros(Bool, nx, ny) for _ in 1:S]

    for i in 1:S
        diff = abs.(climate .- μ[i])
        A[i] = diff .<= σi[i]
    end
    A
end

biotic_suitability(B) = B

function combine_suitability(A, B)
    S = length(A)
    nx, ny = size(A[1])
    C = [A[i] .& B[i] for i in 1:S]
    C
end

function movement_filter(C, k)
    S = length(C)
    [apply_movement(C[i], k) for i in 1:S]
end

end
