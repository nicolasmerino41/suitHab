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
Cgrid = Climate.make_climate_grid(nx, ny; kind=:ridge, seed=11)

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

    # who's a consumer?
    is_cons = mw.trophic_role .!= :basal
    # average abiotic admission among consumers: π_A
    piA_cons = mean(am[is_cons])

    ΔA    = Metrics.delta_area(am, bm)
    ΔG    = Metrics.delta_gini(am, bm)
    Psuff = BAM.prey_sufficiency(mw, out[:Akeep]; kreq=kreq)
    Ks    = BAM.K_spectrum(mw, out[:Akeep])

    return (ΔA=ΔA, ΔG=ΔG, Psuff=Psuff, Kmean=mean(Ks), Kp90=quantile(Ks,0.9), πA=piA_cons)
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

        πAs = [v.πA for v in vals]; πAlo, πAmean, πAhi = Metrics.qband(πAs)
        row = merge(pars, (; ΔAlo, ΔAmean, ΔAhi, ΔGlo, ΔGmean, ΔGhi,
                            Plo, Pmean, Phi, Klo, Kmean, Khi, πAlo, πAmean, πAhi))

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
if true # Only run if needed
    res1 = replicate_sweep(rng, sweep1; fixed=fixed1)
    CSV.write(joinpath(@__DIR__, "data", "sweep_C_align_ridge.csv"), res1)
end
res1 = CSV.read(joinpath(@__DIR__, "data", "sweep_C_align.csv"), DataFrame)

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
    save(joinpath(@__DIR__, "data", "figs", "HM_DeltaArea_C_align_ridge.png"), fig)
end

# ΔArea vs Psuff scatter (mediation)
begin
    fig = Figure(; size=(800,500))
    ax  = Axis(fig[1,1], xlabel="P_suff", ylabel="ΔArea", title="Mediation: ΔArea vs P_suff (C, align sweep)")
    scatter!(ax, res1.Pmean, res1.ΔAmean)
    display(fig)
    save(joinpath(@__DIR__, "data", "figs", "Scatter_DeltaArea_vs_Psuff_ridge.png"), fig)
end

# ---------------------------
# Sweep 2: (R95, σ)
# ---------------------------
sweep2 = NamedTuple[(; R95=r, σ=s) for r in R95s for s in Sigmas]
fixed2 = (; Cgrid=Cgrid, align=default_align, σ=default_σ, R95=default_R95,
           C=default_C, S=S, basal_frac=basal_frac, τA=τA, kreq=kreq)
if true # Only run if you changed something
    res2 = replicate_sweep(rng, sweep2; fixed=fixed2)
    CSV.write(joinpath(@__DIR__, "data", "sweep_R_sigma_ridge.csv"), res2)
end
res2 = CSV.read(joinpath(@__DIR__, "data", "sweep_R_sigma.csv"), DataFrame)

# Heatmap of ΔArea(R95, σ)
begin
    fig = Figure(; size=(900,500))
    ax  = Axis(fig[1,1], xlabel="R95 (diet redundancy)", ylabel="niche breadth σ", title="ΔArea (AM-BAM) — R95 vs σ")
    xs = sort(unique(res2.R95)); ys = sort(unique(res2.σ))
    Z  = [first(res2.ΔAmean[(res2.R95.==x) .& (res2.σ.==y)]) for x in xs, y in ys]
    heatmap!(ax, xs, ys, Z')
    cb = Colorbar(fig[1,2], label="ΔArea")
    display(fig)
    save(joinpath(@__DIR__, "data", "figs", "HM_DeltaArea_R_sigma_ridge.png"), fig)
end

# ---- ΔGini bars at four illustrative corners (robust) ----
# Helper: return the dataframe row in `df` nearest to (Rtarget, σtarget)
nearest_row(df, Rtarget, σtarget) = begin
    d2 = (df.R95 .- Rtarget).^2 .+ (df.σ .- σtarget).^2
    df[argmin(d2), :]
end

# Define the *targets* using min/max of your parameter ranges
corner_targets = [
    (; R95 = minimum(R95s), σ = minimum(Sigmas)),  # low R, narrow σ
    (; R95 = minimum(R95s), σ = maximum(Sigmas)),  # low R, broad σ
    (; R95 = maximum(R95s), σ = minimum(Sigmas)),  # high R, narrow σ
    (; R95 = maximum(R95s), σ = maximum(Sigmas)),  # high R, broad σ
]

labels = ["low R, narrow σ", "low R, broad σ", "high R, narrow σ", "high R, broad σ"]

vals = Float64[]; los = Float64[]; his = Float64[]
for c in corner_targets
    row = nearest_row(res2, c.R95, c.σ)
    push!(vals, row.ΔGmean)
    push!(los,  row.ΔGlo)
    push!(his,  row.ΔGhi)
end

begin
    fig = Figure(; size=(900, 500))
    ax  = Axis(
        fig[1,1],
        ylabel = "ΔGini (BAM − AM)", xticks = (1:4, labels),
        title  = "ΔGini at illustrative corners"
    )
    barplot!(ax, 1:4, vals)
    Plotting.verrorbars!(ax, 1:4, vals, los, his; whisker=0.25)
    hlines!(ax, 0.0; color=:gray, linestyle=:dash)
    ax.xticklabelrotation = π/10  # optional: avoid label overlap
    display(fig)
    save(joinpath(@__DIR__, "data", "figs", "Bars_DeltaGini_corners_ridge.png"), fig)
end

println("Done. CSVs and figures written to data/ and data/figs/")

# --- π_A and P_suff vs alignment at a fixed connectance ---
C_fixed = 0.06  # pick something in your “ridge” zone

# helper: nearest row at (C, align)
function nearest_row_at(df, Ctarget, atarget)
    d2 = (df.C .- Ctarget).^2 .+ (df.align .- atarget).^2
    df[argmin(d2), :]
end

xs   = sort(unique(res1.align))
πA_y = Float64[]; Ps_y = Float64[]
for a in xs
    row = nearest_row_at(res1, C_fixed, a)
    push!(πA_y, row.πAmean)
    push!(Ps_y, row.Pmean)
end
Δ̂_y = πA_y .* (1 .- Ps_y)

begin
    fig = Figure(; size=(1100,420))
    ax1 = Axis(fig[1,1], xlabel="alignment", ylabel="π_A (consumers)", title="Abiotic admission vs alignment")
    lines!(ax1, xs, πA_y)

    ax2 = Axis(fig[1,2], xlabel="alignment", ylabel="P_suff", title="Prey sufficiency vs alignment")
    lines!(ax2, xs, Ps_y)

    ax3 = Axis(fig[1,3], xlabel="alignment", ylabel="ΔArea ≈ π_A × (1−P_suff)", title="Approx ΔArea from mediation")
    lines!(ax3, xs, Δ̂_y)

    display(fig)
    save(joinpath(@__DIR__, "data", "figs", "Mediation_curves_ridge.png"), fig)
end

