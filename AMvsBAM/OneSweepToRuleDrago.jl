using Pkg
cd(pwd())

dir = pwd()

using Random, StatsBase, Statistics, DataFrames, CSV, Printf, GLM
include("src/metaweb.jl");   using .MetaWeb
# include("src/metaweb_old_approach.jl");   using .MetaWeb
include("src/niches.jl");   using .Niches
include("src/bam.jl");      using .BAM
include("src/metrics.jl");  using .Metrics
# include("src/plotting.jl"); using .Plotting
include("src/climate.jl");     using .Climate

# ---------------------------
# CONFIG
# ---------------------------
mkpath(joinpath(@__DIR__, "data"))
mkpath(joinpath(@__DIR__, "data", "figs"))

rng           = MersenneTwister(1)
S             = 175
basal_frac    = 0.30
nx, ny        = 40, 40
tauA          = 0.5
kreq          = 1
replicates    = 10 # 20 is the original
 
# Fixed defaults for slices
default_R95   = 5
default_C     = 0.10
default_sigma = 0.12
default_align = 0.4

Cs     = range(0.005, 0.1; length=10)
Aligns = range(0.0, 1.0; length=10)
R95s   = Int.(range(1.0, 10.0; length=10))
Sigmas = range(0.02, 0.3; length=10)
kreqs  = Int.(range(1.0, 3; length=3))

# climate grid (choose gradient here)
grid_type = "gradient"
Cgrid = Climate.make_climate_grid(nx, ny; kind=Symbol(grid_type), seed=11)

# ---------------------------
# helpers
# ---------------------------
function run_once(rng; Cgrid, align, sigma, R95, motif_mix=:mixed,
                  S=175, basal_frac=0.3, tauA=0.5, kreq=1, connectance=0.10)
    mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                           connectance=connectance, R95=R95,
                           motif_mix=motif_mix, align=align)
    mu, sigma_i = Niches.make_niches(rng, S; align=align, sigma=sigma, basal_frac=basal_frac)
    pars  = BAM.Params(; tauA=tauA, kreq=kreq)
    out   = BAM.compute_AM_BAM(rng, mw, Cgrid, mu, sigma_i, pars)

    am = BAM.species_areas(out[:AM_maps])
    bm = BAM.species_areas(out[:BAM_maps])

    # average abiotic admission among consumers: pi_A
    is_cons   = mw.trophic_role .!= :basal
    vpi       = filter(isfinite, am[is_cons])
    piA_cons  = isempty(vpi) ? NaN : mean(vpi)

    # prey sufficiency (may be NaN if denominator hits zero for some runs)
    Psuff = BAM.prey_sufficiency(mw, out[:Akeep]; kreq=kreq)

    # K spectrum (empty -> NaN quantile)
    Ks    = BAM.K_spectrum(mw, out[:Akeep])
    vK    = filter(isfinite, Ks)
    Kp90  = isempty(vK) ? NaN : quantile(vK, 0.9)

    dA    = Metrics.delta_area(am, bm)
    dG    = Metrics.delta_gini(am, bm)

    return (dA=dA, dG=dG, Psuff=Psuff, Kmean=mean(Ks), Kp90=Kp90, piA=piA_cons)
end

function replicate_sweep(rng, sweep; fixed::NamedTuple, replicates::Int=20)
    n = length(sweep)
    # Pre-initialize with empty vectors for each thread
    results_threads = [NamedTuple[] for _ in 1:Threads.nthreads()]

    Threads.@threads for idx in 1:n
        tid = Threads.threadid()

        pars = sweep[idx]
        # independent RNG per thread, re-seeded for reproducibility
        local_rng = MersenneTwister(rand(rng, 1:10^9))

        vals = [run_once(MersenneTwister(rand(local_rng, 1:10^9));
                         Cgrid = fixed.Cgrid,
                         align = get(pars, :align, fixed.align),
                         sigma = get(pars, :sigma, fixed.sigma),
                         R95 = get(pars, :R95, fixed.R95),
                         connectance = get(pars, :C, fixed.C),
                         motif_mix = get(pars, :motif, :mixed),
                         S=fixed.S, basal_frac=fixed.basal_frac,
                         tauA=fixed.tauA, kreq=fixed.kreq)
                for _ in 1:replicates]

        dAs   = [v.dA for v in vals];  dAlo, dAmean, dAhi = Metrics.qband(dAs)
        dGs   = [v.dG for v in vals];  dGlo, dGmean, dGhi = Metrics.qband(dGs)
        Ps    = [v.Psuff for v in vals]; Plo, Pmean, Phi  = Metrics.qband(Ps)
        K90s  = [v.Kp90 for v in vals]; Klo, Kmean, Khi   = Metrics.qband(K90s)

        pis   = [v.piA for v in vals];  piAlo, piAmean, piAhi = Metrics.qband(pis)

        row = merge(pars, (; dAlo, dAmean, dAhi,
                             dGlo, dGmean, dGhi,
                             Plo, Pmean, Phi,
                             Klo, Kmean, Khi,
                             piAlo, piAmean, piAhi))

        push!(results_threads[tid], row)
    end

    # flatten thread-local results
    flat = reduce(vcat, results_threads)
    return DataFrame(flat)
end

# ---------------------------
# Sweep 1: (C, align)
# ---------------------------
# sweep1 = NamedTuple[(; C=c, align=a) for c in Cs for a in Aligns]
sweep1 = NamedTuple[(; C=c, sigma=s, align=a, R95=R, kreq=k) for c in Cs for s in Sigmas, a in Aligns, R in R95s, k in kreqs]
fixed1 = (; Cgrid=Cgrid, align=default_align, sigma=default_sigma, R95=default_R95,
           C=default_C, S=S, basal_frac=basal_frac, tauA=tauA, kreq=kreq)
# if true # Only run if needed
res1 = replicate_sweep(rng, sweep1; fixed=fixed1)
CSV.write(joinpath(@__DIR__, "data", "OneSweepToRule_$(grid_type)_mixed.csv"), res1)
