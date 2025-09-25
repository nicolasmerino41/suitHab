# ==========================================================
# run_scaling_diagnostics.jl
#
# Purpose
# -------
# Find scale-robust settings for AM vs BAM by:
#   1) adding a "scaleless" fractal climate grid,
#   2) estimating the climate correlation length L,
#   3) working with dimensionless controls:
#        σ̂ = σ / L            (niche breadth vs climate scale)
#        κ  = kreq / R95       (prey requirement vs diet redundancy)
#        λ  = C * S            (mean links per species)
#        ρ  = S / Ncells       (species-to-cell ratio)
#   4) running three decisive diagnostics:
#        A) finite-size convergence (grid resolution),
#        B) species-number scaling (S) at fixed λ and ρ,
#        C) mechanism heatmap ΔArea over (σ̂, κ)
#
# ---- your AM/BAM code ---
include("../SetUp.jl");-
include("src/metaweb.jl");    using .MetaWeb
include("src/niches.jl");     using .Niches
include("src/bam.jl");        using .BAM
include("src/metrics.jl");    using .Metrics
include("src/climate.jl");    using .climate
const Mke = CairoMakie

# ----------------------------------------------------------
# General run helper (one replicate)
# ----------------------------------------------------------
function run_once(rng::AbstractRNG, Cgrid;
                  S::Int, basal_frac::Float64,
                  connectance::Float64, R95::Real,
                  align::Real, σ::Real, τA::Real, kreq::Int,
                  motif_mix::Symbol=:mixed)

    # metaweb + niches
    mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                               connectance=connectance, R95=R95,
                               motif_mix=motif_mix)
    μ, σi = Niches.make_niches(rng, S; align=align, σ=σ, basal_frac=basal_frac)
    pars  = BAM.Params(; τA=τA, kreq=kreq)

    out   = BAM.compute_AM_BAM(rng, mw, Cgrid, μ, σi, pars)
    
    am = BAM.species_areas(out[:AM_maps])
    bm = BAM.species_areas(out[:BAM_maps])

    # who's a consumer?
    is_cons = mw.trophic_role .== :consumer
    # average abiotic admission among consumers: π_A
    piA_cons = mean(am[is_cons])

    ΔA    = Metrics.delta_area(am, bm)
    ΔG    = Metrics.delta_gini(am, bm)
    Psuff = BAM.prey_sufficiency(mw, out[:Akeep]; kreq=kreq)   # mediation metric
    Ks    = BAM.K_spectrum(mw, out[:Akeep])                    # prey co-retention spectrum

    return (ΔA=ΔA, ΔG=ΔG, Psuff=Psuff, Kmean=mean(Ks), Kp90=quantile(Ks,0.9), πA=piA_cons)
end

# summarize over replicates
qband(v) = (quantile(v, 0.05), mean(v), quantile(v, 0.95))

# ----------------------------------------------------------
# CONFIG (baseline)
# ----------------------------------------------------------
mkpath(joinpath(@__DIR__, "data"))
mkpath(joinpath(@__DIR__, "data", "figs"))

rng           = MersenneTwister(11)
S0            = 175
basal_frac    = 0.30
τA            = 0.50
kreq          = 1
replicates    = 20

# default “shape” of the metaweb & niches for slices
C0            = 0.10      # connectance
R95_0         = 5         # diet redundancy
σ0            = 0.12
align0        = 0.5

# grid choices
nx0, ny0      = 40, 40
grid_kinds    = [:gradient, :ridge, :fractal]  # we will use :fractal for “scaleless” diagnostics
Cgrid0        = make_climate_grid(nx0, ny0; kind=:fractal, seed=99)  # base field
L0            = estimate_corr_length(Cgrid0)   # correlation length (cells)
σ̂0           = σ0 / L0

# ----------------------------------------------------------
# A) FINITE-SIZE CONVERGENCE (same physical extent, ↑ resolution)
# ----------------------------------------------------------
res_fs = DataFrame(Ncells=Int[], nx=Int[], S=Int[], C=Float64[], R95=Float64[], σ=Float64[],
                   align=Float64[], L=Float64[], σhat=Float64[], λ=Float64[], ρ=Float64[],
                   ΔAlo=Float64[], ΔAmean=Float64[], ΔAhi=Float64[],
                   ΔGlo=Float64[], ΔGmean=Float64[], ΔGhi=Float64[],
                   Plo=Float64[], Pmean=Float64[], Phi=Float64[])

for (nx, ny) in ((40,40), (80,80), (160,160))
    Cgrid = resample_to(Cgrid0, nx, ny)  # same extent, finer pixels
    L     = estimate_corr_length(Cgrid)
    σhat  = σ0 / L
    λ     = C0 * S0
    ρ     = S0 / (nx*ny)

    vals = [run_once(MersenneTwister(rand(rng, 1:10^9)), Cgrid;
                     S=S0, basal_frac=basal_frac,
                     connectance=C0, R95=R95_0,
                     align=align0, σ=σ0, τA=τA, kreq=kreq)
            for _ in 1:replicates]
    ΔAs = [v.ΔA for v in vals]; ΔGs = [v.ΔG for v in vals]; Ps = [v.Psuff for v in vals]
    ΔAlo, ΔAmean, ΔAhi = qband(ΔAs)
    ΔGlo, ΔGmean, ΔGhi = qband(ΔGs)
    Plo,  Pmean,  Phi  = qband(Ps)

    push!(res_fs, (; Ncells=nx*ny, nx, S=S0, C=C0, R95=R95_0, σ=σ0, align=align0,
                    L, σhat, λ, ρ, ΔAlo, ΔAmean, ΔAhi, ΔGlo, ΔGmean, ΔGhi, Plo, Pmean, Phi))
end

CSV.write(joinpath(@__DIR__, "data", "finite_size.csv"), res_fs)

# plot: ΔArea, ΔGini, P_suff vs Ncells
begin
    fig = Figure(; size=(1100, 340))
    ax1 = Axis(
        fig[1,1],
        xlabel="number of cells", ylabel="ΔArea",
        title="Finite-size convergence (ΔArea)",
        xticklabelsize=10, yticklabelsize=10
    )
    lines!(ax1, res_fs.Ncells, res_fs.ΔAmean); scatter!(ax1, res_fs.Ncells, res_fs.ΔAmean)
    ax2 = Axis(
        fig[1,2],
        xlabel="number of cells", ylabel="ΔGini",
        title="Finite-size convergence (ΔGini)",
        xticklabelsize=10, yticklabelsize=10
        )
    lines!(ax2, res_fs.Ncells, res_fs.ΔGmean); scatter!(ax2, res_fs.Ncells, res_fs.ΔGmean)
    ax3 = Axis(
        fig[1,3],
        xlabel="number of cells", ylabel="P_suff",
        title="Finite-size convergence (P_suff)",
        xticklabelsize=10, yticklabelsize=10
        )
    lines!(ax3, res_fs.Ncells, res_fs.Pmean);  scatter!(ax3, res_fs.Ncells, res_fs.Pmean)
    display(fig)
    save(joinpath(@__DIR__, "data", "figs", "A_FiniteSize_convergence.png"), fig)
end

# ----------------------------------------------------------
# B) SPECIES-NUMBER SCALING (vary S; hold λ and ρ fixed)
# ----------------------------------------------------------
λ_target = C0 * S0                      # keep mean degree
ρ_target = S0 / (nx0 * ny0)             # keep crowding per area

S_list   = [75, 150, 300, 600]          # will adapt grid and C to match targets

res_S = DataFrame(S=Int[], nx=Int[], ny=Int[], Ncells=Int[], C=Float64[],
                  R95=Float64[], σ=Float64[], align=Float64[],
                  L=Float64[], σhat=Float64[], λ=Float64[], ρ=Float64[],
                  ΔAlo=Float64[], ΔAmean=Float64[], ΔAhi=Float64[],
                  ΔGlo=Float64[], ΔGmean=Float64[], ΔGhi=Float64[],
                  Plo=Float64[], Pmean=Float64[], Phi=Float64[])

for S in S_list
    # choose grid so S/Ncells ≈ ρ_target
    Ncells = max(100, round(Int, S / ρ_target))
    nside  = round(Int, sqrt(Ncells))
    nx = nside; ny = nside; Ncells = nx*ny

    # adapt connectance to keep λ = C*S ≈ λ_target
    C = clamp(λ_target / S, 0.01, 0.5)

    # use the same base fractal field at highest resolution; resample to (nx,ny)
    Cgrid = resample_to(Cgrid0, nx, ny)
    L     = estimate_corr_length(Cgrid)
    σhat  = σ0 / L
    λ     = C * S
    ρ     = S / Ncells

    vals = [run_once(MersenneTwister(rand(rng, 1:10^9)), Cgrid;
                     S=S, basal_frac=basal_frac,
                     connectance=C, R95=R95_0,
                     align=align0, σ=σ0, τA=τA, kreq=kreq)
            for _ in 1:replicates]
    ΔAs = [v.ΔA for v in vals]; ΔGs = [v.ΔG for v in vals]; Ps = [v.Psuff for v in vals]
    ΔAlo, ΔAmean, ΔAhi = qband(ΔAs)
    ΔGlo, ΔGmean, ΔGhi = qband(ΔGs)
    Plo,  Pmean,  Phi  = qband(Ps)

    push!(res_S, (; S, nx, ny, Ncells, C, R95=R95_0, σ=σ0, align=align0,
                   L, σhat, λ, ρ, ΔAlo, ΔAmean, ΔAhi, ΔGlo, ΔGmean, ΔGhi, Plo, Pmean, Phi))
end

CSV.write(joinpath(@__DIR__, "data", "species_scaling.csv"), res_S)

# plot: ΔArea, ΔGini, P_suff vs S  (markers show λ & ρ ~ constant)
begin
    fig = Figure(; size=(1100, 340))
    ax1 = Axis(fig[1,1], xlabel="S", ylabel="ΔArea", title="Species scaling (ΔArea)")
    lines!(ax1, res_S.S, res_S.ΔAmean); scatter!(ax1, res_S.S, res_S.ΔAmean)
    ax2 = Axis(fig[1,2], xlabel="S", ylabel="ΔGini", title="Species scaling (ΔGini)")
    lines!(ax2, res_S.S, res_S.ΔGmean); scatter!(ax2, res_S.S, res_S.ΔGmean)
    ax3 = Axis(fig[1,3], xlabel="S", ylabel="P_suff", title="Species scaling (P_suff)")
    lines!(ax3, res_S.S, res_S.Pmean);  scatter!(ax3, res_S.S, res_S.Pmean)
    display(fig)
    save(joinpath(@__DIR__, "data", "figs", "B_Species_scaling.png"), fig)
end

# ----------------------------------------------------------
# C) MECHANISM HEATMAP: ΔArea over (σ̂, κ) at fixed λ
# ----------------------------------------------------------
# Sweep σ and R95, compute σ̂ and κ = kreq / R95; keep S and C so that λ≈λ_target.
Sigmas = collect(range(0.06, 0.24; length=10))
R95s   = Int.(collect(range(1.0, 10.0; length=10)))

res_mech = DataFrame(σ=Float64[], R95=Float64[], σhat=Float64[], κ=Float64[],
                     ΔAlo=Float64[], ΔAmean=Float64[], ΔAhi=Float64[],
                     Plo=Float64[], Pmean=Float64[], Phi=Float64[])

# keep grid fixed (scaleless fractal) to isolate mechanism
Cgrid_mech = resample_to(Cgrid0, 80, 80)
L_mech     = estimate_corr_length(Cgrid_mech)

for σ in Sigmas, R95 in R95s
    σhat = σ / L_mech
    κ    = kreq / R95
    C    = clamp(λ_target / S0, 0.01, 0.5)  # keep λ≈λ_target at S0
    vals = [run_once(MersenneTwister(rand(rng, 1:10^9)), Cgrid_mech;
                     S=S0, basal_frac=basal_frac,
                     connectance=C, R95=R95,
                     align=align0, σ=σ, τA=τA, kreq=kreq)
            for _ in 1:replicates]
    ΔAs = [v.ΔA for v in vals]; Ps = [v.Psuff for v in vals]
    ΔAlo, ΔAmean, ΔAhi = qband(ΔAs)
    Plo,  Pmean,  Phi  = qband(Ps)
    push!(res_mech, (; σ, R95, σhat, κ, ΔAlo, ΔAmean, ΔAhi, Plo, Pmean, Phi))
end

CSV.write(joinpath(@__DIR__, "data", "mechanism_sigmahat_kappa.csv"), res_mech)

# Heatmap ΔArea over (σ̂, κ)
begin
    fig = Figure(; size=(900, 520))
    ax  = Axis(fig[1,1], xlabel="σ̂ = σ / L", ylabel="κ = kreq / R95",
               title="ΔArea (AM − BAM) over (σ̂, κ)")
    xs = sort(unique(res_mech.σhat)); ys = sort(unique(res_mech.κ))
    Z  = [first(res_mech.ΔAmean[(isapprox.(res_mech.σhat, x; atol=1e-8)) .&
                                (isapprox.(res_mech.κ,    y; atol=1e-8))]) for x in xs, y in ys]
    heatmap!(ax, xs, ys, Z')
    Colorbar(fig[1,2], label="ΔArea")
    display(fig)
    save(joinpath(@__DIR__, "data", "figs", "C_HM_DeltaArea_sigmahat_kappa.png"), fig)
end

# Also show mediation: ΔArea vs P_suff for these runs
begin
    fig = Figure(; size=(800, 520))
    ax  = Axis(fig[1,1], xlabel="P_suff", ylabel="ΔArea",
               title="Mediation: ΔArea vs P_suff (σ̂, κ sweep)")
    scatter!(ax, res_mech.Pmean, res_mech.ΔAmean)
    display(fig)
    save(joinpath(@__DIR__, "data", "figs", "C_Scatter_DeltaArea_vs_Psuff.png"), fig)
end

println("Done. Wrote:")
println("  data/finite_size.csv, data/species_scaling.csv, data/mechanism_sigmahat_kappa.csv")
println("  data/figs/A_FiniteSize_convergence.png")
println("  data/figs/B_Species_scaling.png")
println("  data/figs/C_HM_DeltaArea_sigmahat_kappa.png")
println("  data/figs/C_Scatter_DeltaArea_vs_Psuff.png")
