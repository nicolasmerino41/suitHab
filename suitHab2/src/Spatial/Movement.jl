module Movement

export apply_movement

using ..SpatialUtils: label_components

function apply_movement(S::AbstractMatrix{Bool}, k::Int)
    labels, ncomp = label_components(S)
    keep = zeros(Bool, size(S))
    for comp_id in 1:ncomp
        mask = labels .== comp_id
        if sum(mask) ≥ k
            keep .|= mask
        end
    end
    keep
end

end
