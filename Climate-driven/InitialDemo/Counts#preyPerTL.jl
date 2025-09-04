# metaweb trophic levels (before climate)
tl = trophic_level_metaweb(pool)
using StatsBase
println("Metaweb TL counts (consumers): ", countmap(tl[.!pool.basal]))

# which consumer TLs actually have ANY realized presence before loss?
Z0 = climate_pass(pool, grid; τ=τ)
P0 = assemble(Z0, pool)
present = vec(sum(P0; dims=2) .> 0)
println("Realized TL counts (pre-loss): ",
        countmap(tl[.!pool.basal .& present]))
