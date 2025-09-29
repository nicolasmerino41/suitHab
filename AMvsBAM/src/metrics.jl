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

# Quantile band that ignores NaNs/missings; returns (NaN,NaN,NaN) if nothing left
function qband(xs::AbstractVector{<:Real}; lo=0.1, hi=0.9)
    v = Float64[]
    @inbounds for x in xs
        if x !== missing
            y = float(x)
            if isfinite(y) && !isnan(y)
                push!(v, y)
            end
        end
    end
    if isempty(v)
        return (NaN, NaN, NaN)
    end
    return (quantile(v, lo), median(v), quantile(v, hi))
end

end # module
