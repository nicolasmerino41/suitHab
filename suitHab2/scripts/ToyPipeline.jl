S = 100
B = 10
L = 20
T = 0.01

nx, ny = 30, 30

metaweb = build_ppm(S, B, L, T)
A, s = metaweb.A, metaweb.s
visualise_escalator(A, s; B=B)

μ, σ = make_niches(rng, S; σ=0.12)

clim = make_climate_grid(nx, ny; kind=:fractal)
A_layer = abiotic_suitability(clim, μ, σ)
for i in rand(1:S, 5)
    plot_grid(A_layer[i]; title="Abiotic suitability")
end

B0 = build_basal_biotic(B, nx, ny)
for i in rand(1:B, 5)
    plot_grid(B0[i]; title="Biotic suitability")
end

B_layer = propagate_biotic(metaweb, B0)

for i in 31:100
    plot_grid(B_layer[i]; title="Biotic suitability for consumer $i")
end