# ==========================================================
# Run baseline AM vs BAM (no habitat loss / no movement)
# Produces CSVs and PNGs in data/ and data/figs/
# Plots use Makie; each plot is wrapped in begin...display(fig) blocks
# ==========================================================
include("../SetUp.jl");
# include("src/metaweb.jl");   using .MetaWeb
include("src/metaweb_old_approach.jl");   using .MetaWeb
include("src/niches.jl");   using .Niches
include("src/bam.jl");      using .BAM
include("src/metrics.jl");  using .Metrics
include("src/plotting.jl"); using .Plotting
include("src/climate.jl");     using .Climate
const Mke = CairoMakie
# ---------------------------
# CONFIG
# ---------------------------
mkpath(joinpath(@__DIR__, "data"))
mkpath(joinpath(@__DIR__, "data", "figs"))

rng           = MersenneTwister(1)
S             = 175
basal_frac    = 0.30
nx, ny        = 40, 40
τA            = 0.5
kreq          = 1
replicates    = 10 # 20 is the original
 
# Fixed defaults for slices
default_R95   = 5
default_C     = 0.10
default_σ     = 0.12
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
function run_once(rng; Cgrid, align, σ, R95, motif_mix=:mixed,
                  S=175, basal_frac=0.3, τA=0.5, kreq=1, connectance=0.10)
    mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                               connectance=connectance, R95=R95, motif_mix=motif_mix)
    μ, σi = Niches.make_niches(rng, S; align=align, σ=σ, basal_frac=basal_frac)
    pars  = BAM.Params(; τA=τA, kreq=kreq)
    out   = BAM.compute_AM_BAM(rng, mw, Cgrid, μ, σi, pars)

    am = BAM.species_areas(out[:AM_maps])
    bm = BAM.species_areas(out[:BAM_maps])

    # average abiotic admission among consumers: π_A
    is_cons   = mw.trophic_role .!= :basal
    vπ        = filter(isfinite, am[is_cons])
    piA_cons  = isempty(vπ) ? NaN : mean(vπ)

    # prey sufficiency (may be NaN if denominator hits zero for some runs)
    Psuff = BAM.prey_sufficiency(mw, out[:Akeep]; kreq=kreq)

    # K spectrum (empty -> NaN quantile)
    Ks    = BAM.K_spectrum(mw, out[:Akeep])
    vK    = filter(isfinite, Ks)
    Kp90  = isempty(vK) ? NaN : quantile(vK, 0.9)

    ΔA    = Metrics.delta_area(am, bm)
    ΔG    = Metrics.delta_gini(am, bm)

    return (ΔA=ΔA, ΔG=ΔG, Psuff=Psuff, Kmean=mean(Ks), Kp90=Kp90, πA=piA_cons)
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
                         σ = get(pars, :σ, fixed.σ),
                         R95 = get(pars, :R95, fixed.R95),
                         connectance = get(pars, :C, fixed.C),
                         motif_mix = get(pars, :motif, :mixed),
                         S=fixed.S, basal_frac=fixed.basal_frac,
                         τA=fixed.τA, 
                         kreq=get(pars, :kreq, fixed.kreq)
                         )
                for _ in 1:replicates]

        ΔAs   = [v.ΔA for v in vals];  ΔAlo, ΔAmean, ΔAhi = Metrics.qband(ΔAs)
        ΔGs   = [v.ΔG for v in vals];  ΔGlo, ΔGmean, ΔGhi = Metrics.qband(ΔGs)
        Ps    = [v.Psuff for v in vals]; Plo, Pmean, Phi  = Metrics.qband(Ps)
        K90s  = [v.Kp90 for v in vals]; Klo, Kmean, Khi   = Metrics.qband(K90s)

        πAs   = [v.πA for v in vals];  πAlo, πAmean, πAhi = Metrics.qband(πAs)

        row = merge(pars, (; ΔAlo, ΔAmean, ΔAhi,
                             ΔGlo, ΔGmean, ΔGhi,
                             Plo, Pmean, Phi,
                             Klo, Kmean, Khi,
                             πAlo, πAmean, πAhi))

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
sweep1 = NamedTuple[(; C=c, σ=s, align=a, R95=R, kreq=k) for c in Cs for s in Sigmas, a in Aligns, R in R95s, k in kreqs]
fixed1 = (; Cgrid=Cgrid, align=default_align, σ=default_σ, R95=default_R95,
           C=default_C, S=S, basal_frac=basal_frac, τA=τA, kreq=kreq)
# if true # Only run if needed
res1 = replicate_sweep(rng, sweep1; fixed=fixed1)
CSV.write(joinpath(@__DIR__, "data", "OneSweepToRule_$(grid_type)_mixed.csv"), res1)
# end
res1 = CSV.read(joinpath(@__DIR__, "data", "OneSweepToRule_$(grid_type)_mixed.csv"), DataFrame)

