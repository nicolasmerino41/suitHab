module HL

export random_mask, clustered_mask, front_mask, neighbors4

using Random, StatsBase

neighbors4(ix::Int, nx::Int, ny::Int) = begin
    i = ((ix - 1) % nx) + 1
    j = ((ix - 1) ÷ nx) + 1
    out = Int[]
    (i>1)  && push!(out, ix-1)
    (i<nx) && push!(out, ix+1)
    (j>1)  && push!(out, ix-nx)
    (j<ny) && push!(out, ix+nx)
    out
end

function random_mask(rng::AbstractRNG, C::Int, keep_frac::Float64)
    nkeep = clamp(round(Int, keep_frac*C), 0, C)
    keep  = falses(C); keep[randperm(rng, C)[1:nkeep]] .= true; keep
end

function clustered_mask(rng::AbstractRNG, nx::Int, ny::Int, keep_frac::Float64; nseeds::Int=8)
    C = nx*ny; target_remove = C - clamp(round(Int, keep_frac*C), 0, C)
    removed = falses(C); q = Int[]
    seeds = randperm(rng, C)[1:min(nseeds,C)]; append!(q, seeds)
    removed_ct=0; ptr=1
    while removed_ct < target_remove
        if ptr>length(q); push!(q, rand(rng, 1:C)); end
        v = q[ptr]; ptr+=1
        removed[v] && continue
        removed[v]=true; removed_ct+=1
        for nb in neighbors4(v, nx, ny)
            !removed[nb] && push!(q, nb)
        end
    end
    .!removed
end

function front_mask(rng::AbstractRNG, xy::Matrix{Float64}, keep_frac::Float64; axis::Symbol=:x, noise::Float64=0.0)
    C = size(xy,2)
    coord = axis === :x ? view(xy,1,:) : view(xy,2,:)
    q = quantile(coord, 1 - keep_frac)
    keep = BitVector(undef, C)
    if noise ≤ 0
        @inbounds for i in 1:C; keep[i] = coord[i] ≥ q; end
    else
        ϵ = rand(rng, C) .- 0.5
        @inbounds for i in 1:C; keep[i] = coord[i] ≥ (q + noise*ϵ[i]); end
    end
    keep
end

end # module
