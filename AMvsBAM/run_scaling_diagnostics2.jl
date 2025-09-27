# =====================================================================
# build_scaling_diagnostics.jl
# Produces: data/scaling_grid.csv and data/scaling_species.csv
# Requires your existing src/ modules (MetaWeb, Niches, BAM, Metrics, Climate).
# =====================================================================

include("../SetUp.jl");
include("src/metaweb.jl");   using .MetaWeb
include("src/niches.jl");    using .Niches
include("src/bam.jl");       using .BAM
include("src/metrics.jl");   using .Metrics
include("src/climate.jl");   using .Climate

mkpath(joinpath(@__DIR__, "data"))

rng = MersenneTwister(42)

# baseline params
τA       = 0.5
kreq     = 1
C        = 0.10
R95      = 5
σ        = 0.12
align    = 0.4
basal_fr = 0.30
replicates = 12

# ---- one replicate run ----
function run_once(; rng, nx, ny, S, gridkind=:gradient)
    Cgrid = Climate.make_climate_grid(nx, ny; kind=gridkind, seed=11)
    mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_fr, connectance=C, R95=R95, motif_mix=:mixed)
    μ, σi = Niches.make_niches(rng, S; align=align, σ=σ, basal_frac=basal_fr)
    pars  = BAM.Params(; τA=τA, kreq=kreq)
    out   = BAM.compute_AM_BAM(rng, mw, Cgrid, μ, σi, pars)
    am = BAM.species_areas(out[:AM_maps])
    bm = BAM.species_areas(out[:BAM_maps])
    ΔA = Metrics.delta_area(am, bm)
    ΔG = Metrics.delta_gini(am, bm)
    Ps = BAM.prey_sufficiency(mw, out[:Akeep]; kreq=kreq)
    (; ΔA, ΔG, Ps)
end

# ---- grid-size scaling (S fixed) ----
function scaling_grid(; S=175, gridkind=:gradient)
    sizes = [(50,50), (80,80), (160,160)]
    rows = DataFrame(nx=Int[], ny=Int[], ncells=Int[], ΔA=Float64[], ΔG=Float64[], P_suff=Float64[])
    for (nx,ny) in sizes
        vals = [run_once(; rng=MersenneTwister(rand(rng, 1:10^9)), nx, ny, S, gridkind) for _ in 1:replicates]
        push!(rows, (nx=nx, ny=ny, ncells=nx*ny,
                     ΔA=mean(getfield.(vals, :ΔA)),
                     ΔG=mean(getfield.(vals, :ΔG)),
                     P_suff=mean(getfield.(vals, :Ps))))
    end
    rows
end

# ---- species-number scaling (grid fixed) ----
function scaling_species(; nx=100, ny=100, gridkind=:gradient)
    Svals = [100, 175, 300, 600]
    rows = DataFrame(S=Int[], ΔA=Float64[], ΔG=Float64[], P_suff=Float64[])
    for S in Svals
        vals = [run_once(; rng=MersenneTwister(rand(rng, 1:10^9)), nx, ny, S, gridkind) for _ in 1:replicates]
        push!(rows, (S=S,
                     ΔA=mean(getfield.(vals, :ΔA)),
                     ΔG=mean(getfield.(vals, :ΔG)),
                     P_suff=mean(getfield.(vals, :Ps))))
    end
    rows
end

sg = scaling_grid(; S=175, gridkind=:ridge)  # pick whichever grid you want to showcase
ss = scaling_species(; nx=120, ny=120, gridkind=:ridge)

CSV.write(joinpath(@__DIR__, "data", "scaling_grid_ridge.csv"), sg)
CSV.write(joinpath(@__DIR__, "data", "scaling_species_ridge.csv"), ss)
println("Wrote data/scaling_grid.csv and data/scaling_species.csv")
