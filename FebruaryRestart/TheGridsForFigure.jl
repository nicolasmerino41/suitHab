#!/usr/bin/env julia
# Julia 1.11
# Standalone script: plot representative Random vs Autocorrelated grids
using Random
using Statistics
using CairoMakie

# ============================================================
# PARAMETERS (match your main script)
# ============================================================
const NX = 60
const NY = 60
const NCELLS = NX * NY

const E_MIN = 0.0
const E_MAX = 100.0

const AUTOCORR_ITERS = 18
const AUTOCORR_ALPHA = 0.55

Random.seed!(42)

# ============================================================
# Index helpers
# ============================================================

@inline linidx(x::Int, y::Int) = (y - 1) * NX + x
@inline x_of(i::Int) = ((i - 1) % NX) + 1
@inline y_of(i::Int) = ((i - 1) ÷ NX) + 1

function build_neighbors_4()
    neigh = Vector{Vector{Int}}(undef, NCELLS)
    for i in 1:NCELLS
        x = x_of(i); y = y_of(i)
        nb = Int[]
        if x > 1;  push!(nb, linidx(x-1,y)); end
        if x < NX; push!(nb, linidx(x+1,y)); end
        if y > 1;  push!(nb, linidx(x,y-1)); end
        if y < NY; push!(nb, linidx(x,y+1)); end
        neigh[i] = nb
    end
    return neigh
end

const NEIGH_4 = build_neighbors_4()

# ============================================================
# Rescale helper
# ============================================================

function rescale_to_range!(v::Vector{Float64}, lo::Float64, hi::Float64)
    mn = minimum(v)
    mx = maximum(v)
    if mx == mn
        fill!(v, 0.5*(lo+hi))
        return v
    end
    @inbounds for i in eachindex(v)
        v[i] = lo + (hi-lo) * (v[i] - mn) / (mx - mn)
    end
    return v
end

# ============================================================
# Environmental field generator (EXACTLY your logic)
# ============================================================

function smooth_field_once!(E::Vector{Float64}, tmp::Vector{Float64}, α::Float64)
    @inbounds for i in 1:NCELLS
        s = E[i]
        nb = NEIGH_4[i]
        m = s
        for j in nb
            m += E[j]
        end
        m /= (1 + length(nb))
        tmp[i] = (1-α)*E[i] + α*m
    end
    copyto!(E, tmp)
    return E
end

function make_environment(kind::Symbol)
    E = randn(NCELLS)

    if kind == :autocorr
        tmp = similar(E)
        for _ in 1:AUTOCORR_ITERS
            smooth_field_once!(E, tmp, AUTOCORR_ALPHA)
        end
    elseif kind != :random
        error("Unknown kind: $kind")
    end

    rescale_to_range!(E, E_MIN, E_MAX)

    return reshape(E, NX, NY)
end

# ============================================================
# Generate fields
# ============================================================

E_random = make_environment(:random)
E_auto   = make_environment(:autocorr)

# ============================================================
# Plot (same Makie heatmap style)
# ============================================================
begin
    

    f = Figure(size = (1000, 450))

    ax1 = Axis(f[1,1],
        title = "Autocorrelated",
        xticksvisible = false,
        yticksvisible = false
    )

    ax2 = Axis(f[1,2],
        title = "Random",
        xticksvisible = false,
        yticksvisible = false
    )

    h1 = heatmap!(ax1, E_auto)
    h2 = heatmap!(ax2, E_random)

    # Shared color range
    # colorrange!(h1, 0, 100)
    # colorrange!(h2, 0, 100)

    Colorbar(f[1,3], h1, label = "Environmental value")
    display(f)
    save("SI_environmental_grids_julia.png", f)

end

println("Saved: SI_environmental_grids_julia.png")