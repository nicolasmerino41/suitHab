module Metrics

using Statistics

export mean_area, gini, delta_area, delta_gini, qband

mean_area(v::AbstractVector{<:Real}) = mean(v)

function gini(x::AbstractVector{<:Real})
    n = length(x)
    n==0 && return 0.0
    s = sort(x)
    c = cumsum(s)
    G = (n + 1 - 2*sum(c) / c[end]) / n
    return G
end

delta_area(am::AbstractVector{<:Real}, bam::AbstractVector{<:Real}) = mean(am) - mean(bam)
delta_gini(am::AbstractVector{<:Real}, bam::AbstractVector{<:Real}) = gini(am) - gini(bam)

qband(xs::AbstractVector{<:Real}) = (quantile(xs,0.1), mean(xs), quantile(xs,0.9))

end # module
