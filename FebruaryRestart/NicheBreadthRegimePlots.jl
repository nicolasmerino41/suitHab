#!/usr/bin/env julia
using CairoMakie
using Random

# ============================================================
# Parameters matching your pipeline
# ============================================================

const E_MIN = 0.0
const E_MAX = 100.0
const S_SHOW = 10          # how many species to draw per panel
const RNG_SEED = 123

# Same regimes as your model
struct BreadthRegime
    name::String
    meanσ::Float64
    logsd::Float64
end

const REGIMESS = [
    BreadthRegime("Narrow + LowVar",  7.5, 0.20),
    BreadthRegime("Narrow + HighVar", 7.5, 0.55),
    BreadthRegime("Broad + LowVar",   16.0, 0.20),
    BreadthRegime("Broad + HighVar",  16.0, 0.55),
]

# ============================================================
# Helpers
# ============================================================

const SQRT2PI = sqrt(2π)

@inline function normal_pdf(x, μ, σ)
    invσ = 1 / max(σ, 1e-6)
    return (invσ / SQRT2PI) .* exp.(-0.5 .* ((x .- μ) .* invσ) .^ 2)
end

function draw_sigmas(rng::AbstractRNG, regime::BreadthRegime, n::Int)
    σ = Vector{Float64}(undef, n)
    for i in 1:n
        val = exp(log(regime.meanσ) + regime.logsd * randn(rng))
        σ[i] = clamp(val, 1.5, 45.0)   # same clamp as pipeline
    end
    return σ
end

# ============================================================
# Global x (environmental gradient)
# ============================================================

x = range(E_MIN, E_MAX; length=2000)

# ============================================================
# Plot
# ============================================================

rng = MersenneTwister(RNG_SEED)

f = Figure(size=(1200, 900))
Label(f[0, :], "Niche breadth regimes (μ uniform, σ lognormal) — faithful to pipeline", fontsize=20)

axes = Axis[]

for (k, reg) in enumerate(REGIMESS)
    r = (k - 1) ÷ 2 + 1
    c = (k - 1) % 2 + 1
    ax = Axis(f[r, c])
    push!(axes, ax)

    # μ uniform across full domain (as in simulate_one!)
    μ = rand(rng, S_SHOW) .* (E_MAX - E_MIN) .+ E_MIN

    # σ from regime
    σ = draw_sigmas(rng, reg, S_SHOW)

    # plot S_SHOW species
    for i in 1:S_SHOW
        y = normal_pdf(x, μ[i], σ[i])
        lines!(ax, x, y)
    end

    ax.title = reg.name
    ax.xlabel = "Environment"
    ax.ylabel = "Suitability density"
end

# ============================================================
# Force identical x only (y free)
# ============================================================

for ax in axes
    xlims!(ax, E_MIN, E_MAX)
end
linkxaxes!(axes...)

save("breadth_regimes_realistic.png", f)
display(f)

println("Saved: breadth_regimes_realistic.png")
