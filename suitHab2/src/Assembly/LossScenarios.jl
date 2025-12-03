module LossScenarios

export apply_loss

using Random

function apply_loss(F, frac; kind=:random, seed=1)
    rng = MersenneTwister(seed)
    nx, ny = size(F[1])
    total = nx * ny
    nremove = round(Int, frac * total)

    mask = ones(Bool, nx, ny)

    if kind == :random
        inds = shuffle(rng, collect(1:total))[1:nremove]
        for idx in inds
            i = fld(idx-1, ny) + 1
            j = mod(idx-1, ny) + 1
            mask[i,j] = false
        end
    end

    [F[i] .& mask for i in 1:length(F)]
end

end
