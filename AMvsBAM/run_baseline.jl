# ==========================================================
# Run baseline AM vs BAM (no habitat loss / no movement)
# Produces CSVs and PNGs in data/ and data/figs/
# Plots use Makie; each plot is wrapped in begin...display(fig) blocks
# ==========================================================

include("../SetUp.jl");
include("src/metaweb.jl");   using .MetaWeb
include("src/niches.jl");   using .Niches
include("src/bam.jl");      using .BAM
include("src/metrics.jl");  using .Metrics
include("src/plotting.jl"); using .Plotting
const Mke = CairoMakie
# ---------------------------
# CONFIG
# ---------------------------
mkpath(joinpath(@__DIR__, "..", "data"))
mkpath(joinpath(@__DIR__, "..", "data", "figs"))

rng           = MersenneTwister(1)
S             = 175
basal_frac    = 0.30
nx, ny        = 40, 40
τA            = 0.5
kreq          = 1
replicates    = 20

# Fixed defaults for slices
default_R95   = 5
default_C     = 0.10
default_σ     = 0.12
default_align = 0.4

Cs     = range(0.01, 0.20; length=30)
Aligns = range(0.0, 1.0; length=30)
R95s   = Int.(range(1.0, 10.0; length=10))
Sigmas = range(0.06, 0.24; length=24)


# climate grid (choose gradient here)
Cgrid = Niches.make_climate_grid(nx, ny; kind=:nongradient, seed=11)

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

    ΔA    = Metrics.delta_area(am, bm)
    ΔG    = Metrics.delta_gini(am, bm)
    Psuff = BAM.prey_sufficiency(mw, out[:Akeep]; kreq=kreq)
    Ks    = BAM.K_spectrum(mw, out[:Akeep])

    return (ΔA=ΔA, ΔG=ΔG, Psuff=Psuff, Kmean=mean(Ks), Kp90=quantile(Ks, 0.9))
end

function replicate_sweep(rng, sweep; fixed::NamedTuple)
    results = DataFrame()
    for pars in sweep
        vals = [run_once(MersenneTwister(rand(rng, 1:10^9));
                    Cgrid = fixed.Cgrid,
                    align = get(pars, :align, fixed.align),
                    σ = get(pars, :σ, fixed.σ),
                    R95 = get(pars, :R95, fixed.R95),
                    connectance = get(pars, :C, fixed.C),
                    motif_mix = get(pars, :motif, :mixed),
                    S=fixed.S, basal_frac=fixed.basal_frac, τA=fixed.τA, kreq=fixed.kreq)
                for _ in 1:replicates]
        # summarize
        ΔAs   = [v.ΔA for v in vals];  ΔAlo, ΔAmean, ΔAhi = Metrics.qband(ΔAs)
        ΔGs   = [v.ΔG for v in vals];  ΔGlo, ΔGmean, ΔGhi = Metrics.qband(ΔGs)
        Ps    = [v.Psuff for v in vals]; Plo, Pmean, Phi  = Metrics.qband(Ps)
        K90s  = [v.Kp90 for v in vals]; Klo, Kmean, Khi   = Metrics.qband(K90s)

        row = merge(pars, (; ΔAlo, ΔAmean, ΔAhi, ΔGlo, ΔGmean, ΔGhi, Plo, Pmean, Phi, Klo, Kmean, Khi))
        push!(results, row, promote=true)
    end
    return results
end

# ---------------------------
# Sweep 1: (C, align)
# ---------------------------
sweep1 = NamedTuple[(; C=c, align=a) for c in Cs for a in Aligns]
fixed1 = (; Cgrid=Cgrid, align=default_align, σ=default_σ, R95=default_R95,
           C=default_C, S=S, basal_frac=basal_frac, τA=τA, kreq=kreq)
res1 = replicate_sweep(rng, sweep1; fixed=fixed1)
CSV.write(joinpath(@__DIR__, "..", "data", "sweep_C_align.csv"), res1)

# Heatmap of ΔArea(C, align)
begin
    fig = Figure(; size=(900,500))
    ax  = Axis(fig[1,1], xlabel="connectance C", ylabel="alignment", title="ΔArea (AM-BAM) — C vs align")
    # grid
    xs = sort(unique(res1.C)); ys = sort(unique(res1.align))
    Z  = [first(res1.ΔAmean[(res1.C.==x) .& (res1.align.==y)]) for x in xs, y in ys]
    heatmap!(ax, xs, ys, Z')
    cb = Colorbar(fig[1,2], label="ΔArea")
    display(fig)
    save(joinpath(@__DIR__, "..", "data", "figs", "HM_DeltaArea_C_align.png"), fig)
end

# ΔArea vs Psuff scatter (mediation)
begin
    fig = Figure(; size=(800,500))
    ax  = Axis(fig[1,1], xlabel="P_suff", ylabel="ΔArea", title="Mediation: ΔArea vs P_suff (C, align sweep)")
    scatter!(ax, res1.Pmean, res1.ΔAmean)
    display(fig)
    save(joinpath(@__DIR__, "..", "data", "figs", "Scatter_DeltaArea_vs_Psuff.png"), fig)
end

# ---------------------------
# Sweep 2: (R95, σ)
# ---------------------------
sweep2 = NamedTuple[(; R95=r, σ=s) for r in R95s for s in Sigmas]
fixed2 = (; Cgrid=Cgrid, align=default_align, σ=default_σ, R95=default_R95,
           C=default_C, S=S, basal_frac=basal_frac, τA=τA, kreq=kreq)
res2 = replicate_sweep(rng, sweep2; fixed=fixed2)
CSV.write(joinpath(@__DIR__, "..", "data", "sweep_R_sigma.csv"), res2)

# Heatmap of ΔArea(R95, σ)
begin
    fig = Figure(; size=(900,500))
    ax  = Axis(fig[1,1], xlabel="R95 (diet redundancy)", ylabel="niche breadth σ", title="ΔArea (AM-BAM) — R95 vs σ")
    xs = sort(unique(res2.R95)); ys = sort(unique(res2.σ))
    Z  = [first(res2.ΔAmean[(res2.R95.==x) .& (res2.σ.==y)]) for x in xs, y in ys]
    heatmap!(ax, xs, ys, Z')
    cb = Colorbar(fig[1,2], label="ΔArea")
    display(fig)
    save(joinpath(@__DIR__, "..", "data", "figs", "HM_DeltaArea_R_sigma.png"), fig)
end

# ΔGini bars at four illustrative corners
corners = [(; R95=min(R95s...), σ=min(Sigmas...)),
           (; R95=min(R95s...), σ=max(Sigmas...)),
           (; R95=max(R95s...), σ=min(Sigmas...)),
           (; R95=max(R95s...), σ=max(Sigmas...))]

labels = ["low R, narrow σ", "low R, broad σ", "high R, narrow σ", "high R, broad σ"]
vals   = Float64[]; los = Float64[]; his = Float64[]

for c in corners
    row = res2[(isapprox.(res2.R95, c.R95; atol=1e-8)) .& 
               (isapprox.(res2.σ,  c.σ;  atol=1e-8)), :]
    if nrow(row) == 0
        @warn "No match found for corner R95=$(c.R95), σ=$(c.σ)"
        continue
    end
    push!(vals, first(row.ΔGmean))
    push!(los,  first(row.ΔGlo))
    push!(his,  first(row.ΔGhi))
end

begin
    fig = Figure(; size=(900,500))
    ax  = Axis(fig[1,1], ylabel="ΔGini (AM − BAM)", xticks=(1:4, labels), title="Inequality shift across corners")
    barplot!(ax, 1:4, vals)
    Plotting.verrorbars!(ax, 1:4, vals, los, his; whisker=0.2)
    hlines!(ax, 0.0; color=:gray, linestyle=:dash)
    display(fig)
    save(joinpath(@__DIR__, "..", "data", "figs", "Bars_DeltaGini_corners.png"), fig)
end

println("Done. CSVs and figures written to data/ and data/figs/")
