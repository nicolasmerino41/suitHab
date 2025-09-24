# ==========================================================
# scripts/run_diagnostics.jl
# Diagnostics for AM vs BAM (no HL / no movement)
# - K ECDF across A-eligible cells
# - P_suff bars with error bars
# - Per-species area ECDFs (AM vs BAM)
# - AM vs BAM scatter for baseline
# Figures -> data/figs/
# ==========================================================
using Random, Statistics
using CSV, DataFrames
using CairoMakie

const Mke = CairoMakie
include("../SetUp.jl");
include("src/metaweb.jl");   using .MetaWeb
include("src/niches.jl");   using .Niches
include("src/bam.jl");      using .BAM
include("src/metrics.jl");  using .Metrics
include("src/plotting.jl"); using .Plotting

# ---------------------------
# CONFIG
# ---------------------------
mkpath(joinpath(@__DIR__, "..", "data", "figs"))

rng           = MersenneTwister(2025)
S             = 175
basal_frac    = 0.30
nx, ny        = 40, 40

# abiotic and prey gate
τA            = 0.50
kreq          = 1

# climate grid (keep fixed; no HL)
Cgrid = Niches.make_climate_grid(nx, ny; kind=:nongradient, seed=99)

# Replicates to stabilize diagnostics
replicates    = 20

# Three contrasting scenarios
scenarios = [
    (name = "baseline",  C = 0.10, align = 0.40, R95 = 5, σ = 0.12, motif = :mixed),
    (name = "lowR-lowA", C = 0.06, align = 0.00, R95 = 2, σ = 0.08, motif = :chains),
    (name = "highR-hiA", C = 0.20, align = 0.80, R95 = 8, σ = 0.20, motif = :omnivory),
]

# ---------------------------
# helpers
# ---------------------------
function run_once(rng, scen; Cgrid, S, basal_frac, τA, kreq)
    mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                               connectance=scen.C, R95=scen.R95, motif_mix=scen.motif)
    μ, σi = Niches.make_niches(rng, S; align=scen.align, σ=scen.σ, basal_frac=basal_frac)
    pars  = BAM.Params(; τA=τA, kreq=kreq)
    out   = BAM.compute_AM_BAM(rng, mw, Cgrid, μ, σi, pars)

    am = BAM.species_areas(out[:AM_maps])
    bm = BAM.species_areas(out[:BAM_maps])

    Ps  = BAM.prey_sufficiency(mw, out[:Akeep]; kreq=kreq)
    Ks  = BAM.K_spectrum(mw, out[:Akeep])

    return (; am, bm, Ps, Ks)
end

function collect_diag(rng, scen; replicates, Cgrid, S, basal_frac, τA, kreq)
    all_am = Float64[]; all_bm = Float64[]
    all_K  = Int[];     all_Ps = Float64[]
    for r in 1:replicates
        seed = rand(rng, 1:10^9)
        v = run_once(MersenneTwister(seed), scen; Cgrid=Cgrid, S=S, basal_frac=basal_frac, τA=τA, kreq=kreq)
        append!(all_am, v.am)
        append!(all_bm, v.bm)
        append!(all_K,  v.Ks)
        push!(all_Ps,   v.Ps)
    end
    return (; am = all_am, bm = all_bm, K = all_K, Ps = all_Ps)
end

# ---------------------------
# RUN
# ---------------------------
diag = Dict{String,Any}()
for scen in scenarios
    diag[scen.name] = collect_diag(rng, scen; replicates, Cgrid, S, basal_frac, τA, kreq)
end

# ---------------------------
# PLOTS
# ---------------------------

# 1) ECDF of K across scenarios
begin
    fig = Figure(; size=(900, 520))
    ax  = Axis(fig[1,1], xlabel = "K = # prey present in A-eligible cell",
                        ylabel = "ECDF",
                        title  = "Prey co-retention (K) across scenarios")

    for (i, scen) in enumerate(scenarios)
        d = diag[scen.name]
        Plotting.ecdf_plot!(ax, Float64.(d.K); label=scen.name)
    end
    Mke.axislegend(ax, position=:rb, framevisible=false)
    display(fig)
    save(joinpath(@__DIR__, "..", "data", "figs", "Diag_K_ECDF.png"), fig)
end

# 2) Bars of P_suff with error bars
begin
    fig = Figure(; size=(900, 520))
    ax  = Axis(fig[1,1], ylabel = "P_suff = P(≥k prey | A-kept cell)",
                        xticks = (1:length(scenarios), [s.name for s in scenarios]),
                        title  = "Prey sufficiency among A-eligible cells")

    means = [mean(diag[s.name].Ps) for s in scenarios]
    los   = [quantile(diag[s.name].Ps, 0.1) for s in scenarios]
    his   = [quantile(diag[s.name].Ps, 0.9) for s in scenarios]

    barplot!(ax, 1:length(means), means)
    Plotting.verrorbars!(ax, 1:length(means), means, los, his; whisker=0.25)
    display(fig)
    save(joinpath(@__DIR__, "..", "data", "figs", "Diag_Psuff_bars.png"), fig)
end

# 3) ECDFs of per-species areas (AM vs BAM) for each scenario
begin
    fig = Figure(; size=(1200, 420))
    for (col, scen) in enumerate(scenarios)
        ax = Axis(fig[1, col], xlabel = "per-species suitable area", ylabel = col==1 ? "ECDF" : "",
                            title  = scen.name)
        Plotting.ecdf_plot!(ax, diag[scen.name].am; label="AM")
        Plotting.ecdf_plot!(ax, diag[scen.name].bm; label="BAM")
        Mke.axislegend(ax, position=:rb, framevisible=false)
    end
    display(fig)
    save(joinpath(@__DIR__, "..", "data", "figs", "Diag_Area_ECDFs.png"), fig)
end

# 4) AM vs BAM per-species scatter (baseline)
begin
    base = diag["baseline"]
    fig  = Figure(; size=(700, 600))
    ax   = Axis(fig[1,1], xlabel="AM per-species area", ylabel="BAM per-species area",
                         title="AM vs BAM per-species area (baseline)")
    scatter!(ax, base.am, base.bm)
    # 1:1 guide
    m = max(maximum(base.am), maximum(base.bm))
    lines!(ax, [0.0, m], [0.0, m]; linestyle=:dash, color=:gray)
    display(fig)
    save(joinpath(@__DIR__, "..", "data", "figs", "Diag_AM_vs_BAM_scatter_baseline.png"), fig)
end

println("Diagnostics written to data/figs/:")
println(" - Diag_K_ECDF.png")
println(" - Diag_Psuff_bars.png")
println(" - Diag_Area_ECDFs.png")
println(" - Diag_AM_vs_BAM_scatter_baseline.png")
