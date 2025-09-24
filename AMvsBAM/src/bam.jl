module BAM

using Random, Statistics
using ..MetaWeb: Metaweb
using ..Niches

export Params, compute_AM_BAM, prey_sufficiency, K_spectrum, species_areas

"""
    Params(; τA=0.5, kreq=1)

Minimal parameter set for BAM without movement / habitat loss.
"""
Base.@kwdef struct Params
    τA::Float64 = 0.5
    kreq::Int   = 1
end

# Gaussian niche response at climate value c
_gauss(c, μ, σ) = exp(-0.5*((c-μ)/σ)^2)

"""
    compute_AM_BAM(rng, mw, Cgrid, μ, σi, pars)

Return Dict with:
- `AM_maps` and `BAM_maps`: S-element vectors of Bool grids (nx×ny).
- `prey_present`: S-element vector of prey-presence grids (for consumers).
"""
function compute_AM_BAM(rng::AbstractRNG, mw::Metaweb, Cgrid::AbstractMatrix,
                        μ::Vector{Float64}, σi::Vector{Float64}, pars::Params)
    S = mw.S
    nx, ny = size(Cgrid)
    # abiotic pass per species
    Akeep = [trues(nx,ny) for _ in 1:S]
    for i in 1:S
        for x in 1:nx, y in 1:ny
            Akeep[i][x,y] = _gauss(Cgrid[x,y], μ[i], σi[i]) >= pars.τA
        end
    end
    # AM maps
    AM_maps = Akeep

    # prey presence per consumer (prey that also pass A in cell)
    prey_present = [falses(nx,ny) for _ in 1:S]
    for j in 1:S
        # prey indices for j
        prey = findall(q->mw.A[j,q], 1:S)
        if isempty(prey)
            prey_present[j] .= false
        else
            tmp = falses(nx,ny)
            for q in prey
                tmp .|= Akeep[q]
            end
            prey_present[j] .= tmp
        end
    end

    # BAM maps: consumer cell is kept if Akeep and at least kreq prey present
    BAM_maps = [falses(nx,ny) for _ in 1:S]
    for j in 1:S
        prey = findall(q->mw.A[j,q], 1:S)
        if isempty(prey)
            BAM_maps[j] .= false
            continue
        end
        # K(x,y): number of prey with Akeep at (x,y)
        K = fill(0, nx, ny)
        for q in prey
            K .+= Int.(Akeep[q])
        end
        BAM_maps[j] .= Akeep[j] .& (K .>= pars.kreq)
    end

    return Dict(:AM_maps=>AM_maps, :BAM_maps=>BAM_maps, :prey_present=>prey_present, :Akeep=>Akeep)
end

"""
    species_areas(maps) -> area vector (S,)

Fraction of grid cells kept for each species (mean across cells).
"""
function species_areas(maps::Vector{BitMatrix})
    nx, ny = size(maps[1])
    S = length(maps)
    out = zeros(Float64, S)
    for i in 1:S
        out[i] = count(maps[i]) / (nx*ny)
    end
    return out
end

"""
    prey_sufficiency(mw, Akeep; kreq=1)

P_suff = P( >=k prey present | A-kept cell for the consumer ), averaged over consumers.
"""
function prey_sufficiency(mw::Metaweb, Akeep::Vector{BitMatrix}; kreq::Int=1)
    S = mw.S
    nx, ny = size(Akeep[1])
    p = 0.0
    n = 0
    for j in 1:S
        prey = findall(q->mw.A[j,q], 1:S)
        isempty(prey) && continue
        K = fill(0, nx, ny)
        for q in prey
            K .+= Int.(Akeep[q])
        end
        kept = Akeep[j]
        num = count((x)->x, (K .>= kreq) .& kept)
        den = count(kept)
        if den>0
            p += num/den
            n += 1
        end
    end
    return n==0 ? 0.0 : p/n
end

"""
    K_spectrum(mw, Akeep)

Return a vector of K values pooled over all consumers and A-eligible cells.
"""
function K_spectrum(mw::Metaweb, Akeep::Vector{BitMatrix})
    S = mw.S
    nx, ny = size(Akeep[1])
    Ks = Int[]
    for j in 1:S
        prey = findall(q->mw.A[j,q], 1:S)
        isempty(prey) && continue
        K = fill(0, nx, ny)
        for q in prey
            K .+= Int.(Akeep[q])
        end
        for x in 1:nx, y in 1:ny
            if Akeep[j][x,y]
                push!(Ks, K[x,y])
            end
        end
    end
    return Ks
end

end # module
