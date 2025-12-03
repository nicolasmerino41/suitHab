module MetawebUtils

export trophic_levels, trophic_incoherence

using LinearAlgebra, Statistics

function trophic_levels(A::Matrix{Int})
    S = size(A,1)
    kin = sum(A, dims=2)
    v = [max(kin[i],1) for i in 1:S]
    Λ = Diagonal(v) - A
    Λ \ v
end

function trophic_incoherence(A, s)
    xs = Float64[]
    for i in 1:length(s), j in 1:length(s)
        if A[i,j] == 1
            push!(xs, s[i] - s[j])
        end
    end
    sqrt(mean(xs.^2) - 1)
end

end
