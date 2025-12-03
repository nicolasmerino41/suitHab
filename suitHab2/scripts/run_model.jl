include("../../SetUp.jl");
include("../src/suitHab.jl")
using .suitHab

rng = MersenneTwister(1)

S = 500
B = 24
L = 20
T = 1.2

metaweb = build_ppm(S, B, L, T)
total_links_possible = S * (S - 1)
realized_links = count(==(1), metaweb.A)
connectance = realized_links / total_links_possible

nx = 60
ny = 60
climate = make_climate_grid(nx, ny; kind=:ridge)

μ, σi = make_niches(rng, S; σ=0.1)

A = abiotic_suitability(climate, μ, σi)

B_basal = build_basal_biotic(length(metaweb.basal), nx, ny)
B_layer = propagate_biotic(metaweb, B_basal)

S_layer = combine_suitability(A, B_layer)

F = movement_filter(S_layer, 4)

R = realized_metaweb(metaweb, F)

println("Realized metaweb links: ", sum(R))

for i in 24:30
    plot_grid(B_layer[i])
end