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
replicates    = 20 # 20 is the original
 
# Fixed defaults for slices
default_R95   = 5
default_C     = 0.10
default_σ     = 0.12
default_align = 0.4

Cs     = range(0.01, 0.20; length=20)
Aligns = range(0.0, 1.0; length=20)
R95s   = Int.(range(1.0, 10.0; length=10))
Sigmas = range(0.06, 0.24; length=24)

# climate grid (choose gradient here)
Cgrid = Climate.make_climate_grid(nx, ny; kind=:fractal, seed=11)

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
if false # Only run if needed
    res1 = replicate_sweep(rng, sweep1; fixed=fixed1)
    CSV.write(joinpath(@__DIR__, "data", "sweep_C_align_gradient.csv"), res1)
end
res1 = CSV.read(joinpath(@__DIR__, "data", "sweep_C_align_fractal.csv"), DataFrame)

# Heatmap of ΔArea(C, align)
begin
    fig = Figure(; size=(900,500))
    ax  = Axis(fig[1,1], xlabel="connectance C", ylabel="alignment", title="ΔArea (AM-BAM) — C vs align")
    # grid
    xs = sort(unique(res1.C)); ys = sort(unique(res1.align))
    Z  = [first(res1.ΔAmean[(res1.C.==x) .& (res1.align.==y)]) for x in xs, y in ys]
    lo, hi = extrema(Z)
    heatmap!(ax, xs, ys, Z')
    cb = Colorbar(fig[1,2], label="ΔArea", colorrange=(lo, hi))
    display(fig)
    # save(joinpath(@__DIR__, "data", "figs", "HM_DeltaArea_C_align_gradient.png"), fig)
end

# ΔArea vs Psuff scatter (mediation)
begin
    fig = Figure(; size=(800,500))
    ax  = Axis(fig[1,1], xlabel="P_suff", ylabel="ΔArea", title="Mediation: ΔArea vs P_suff (C, align sweep)")
    scatter!(ax, res1.Pmean, res1.ΔAmean)
    display(fig)
    # save(joinpath(@__DIR__, "data", "figs", "Scatter_DeltaArea_vs_Psuff_gradient.png"), fig)
end

# ---------------------------
# Sweep 2: (R95, σ)
# ---------------------------
sweep2 = NamedTuple[(; R95=r, σ=s) for r in R95s for s in Sigmas]
fixed2 = (; Cgrid=Cgrid, align=default_align, σ=default_σ, R95=default_R95,
           C=default_C, S=S, basal_frac=basal_frac, τA=τA, kreq=kreq)
if false # Only run if you changed something
    res2 = replicate_sweep(rng, sweep2; fixed=fixed2)
    CSV.write(joinpath(@__DIR__, "data", "sweep_R_sigma_gradient.csv"), res2)
end
res2 = CSV.read(joinpath(@__DIR__, "data", "sweep_R_sigma_ridge.csv"), DataFrame)
############## TEMPORARY: just for testing
# --- SAFE R95–σ HEATMAP BUILDER ---------------------------------------------
# 1) one row per (R95, σ) with mean ΔA (even if you had replicates written to the CSV)
res2_clean = combine(groupby(res2, [:R95, :σ]),
                     :ΔAmean => mean => :ΔA)

# 2) sorted axes
xs = sort(unique(res2_clean.R95))            # Int
ys = sort(unique(res2_clean.σ))              # Float64

# 3) build a dense matrix with no silent fallbacks
Z = fill(NaN, length(xs), length(ys))
let lut = Dict( (res2_clean.R95[i], res2_clean.σ[i]) => res2_clean.ΔA[i]
                for i in 1:nrow(res2_clean) )
    for (ix, x) in pairs(xs), (iy, y) in pairs(ys)
        v = get(lut, (x, y), NaN)
        Z[ix, iy] = v
    end
end

# 4) checks: no gaps, unique combos
@assert all(!isnan, Z) "Missing (R95, σ) combinations in res2."
@assert nrow(res2_clean) == length(xs)*length(ys) "Duplicate or missing grid points."

# 5) plot (note: Makie expects Z' if you want x across columns, y up rows)
begin
    fig = Figure(; size=(980,520))
    ax  = Axis(fig[1,1],
        xlabel = "R95 (diet redundancy)",
        ylabel = "niche breadth σ",
        title  = "ΔArea (AM − BAM) — R95 vs σ")

    heatmap!(ax, xs, ys, permutedims(Z))   # <- transpose for x by y layout
    Colorbar(fig[1,2], label="ΔArea")
    display(fig)
    # save(joinpath(@__DIR__, "data", "figs", "HM_DeltaArea_R_sigma_FIXED.png"), fig)
end

using GLM, StatsModels
m = lm(@formula(ΔA ~ 1 + R95 + σ), res2_clean)   # or add R95*σ if you want the interaction
coef(m)  # look at the sign and magnitude on R95

# after you build xs, ys, Z exactly as you posted
@show size(Z)             # should be (length(xs), length(ys))
@assert size(Z) == (length(xs), length(ys))

# Column-wise means "over σ", i.e., trend along R95
col_means = vec(mean(Z, dims=2))
@show xs
@show col_means           # should be ≈ monotone DOWN (small slope)

# Row-wise means "over R95", i.e., trend along σ
row_means = vec(mean(Z, dims=1))
@show ys
@show row_means           # should be clearly increasing with σ
begin
    lo, hi = extrema(Z)                 # e.g., (0.03, 0.20)

    fig = Figure(; size=(980,520))
    ax  = Axis(
        fig[1,1],
        xlabel = "R95 (diet redundancy)",
        ylabel = "niche breadth σ",
        title  = "ΔArea (AM − BAM) — R95 vs σ"
    )

    hm = heatmap!(ax, xs, ys, Z; colorrange=(lo, hi))   # <-- key line
    Colorbar(fig[1,2], hm, label="ΔArea", ticks=round.(range(lo, hi; length=6), digits=2))
    display(fig)

end

##########################################
# Heatmap of ΔArea(R95, σ)
begin
    fig = Figure(; size=(900,500))
    ax  = Axis(fig[1,1], xlabel="R95 (diet redundancy)", ylabel="niche breadth σ", title="ΔArea (AM-BAM) — R95 vs σ")
    xs = sort(unique(res2.R95)); ys = sort(unique(res2.σ))
    Z  = [first(res2.ΔAmean[(res2.R95.==x) .& (res2.σ.==y)]) for x in xs, y in ys]
    heatmap!(ax, xs, ys, Z')
    cb = Colorbar(fig[1,2], label="ΔArea")
    display(fig)
    # save(joinpath(@__DIR__, "data", "figs", "HM_DeltaArea_R_sigma_gradient.png"), fig)
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
    # save(joinpath(@__DIR__, "data", "figs", "Bars_DeltaGini_corners_gradient.png"), fig)
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
    # save(joinpath(@__DIR__, "data", "figs", "Mediation_curves_gradient.png"), fig)
end

# --- realized diet breadth vs R95 ---
function mean_indegree_vs_R95(; rng=MersenneTwister(1), Rs=1:10, S=175, basal_frac=0.3, C=0.10)
    means = Float64[]
    for r in Rs
        mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac, connectance=C, R95=r, motif_mix=:mixed)
        indeg = [sum(mw.A[:,j]) for j in 1:size(mw.A,2)]               # in-degree per consumer column
        is_cons = mw.trophic_role .!= :basal
        push!(means, mean(indeg[is_cons]))
    end
    Rs, means
end

Rs, mdeg = mean_indegree_vs_R95()
@info "Realized mean indegree by R95" (Rs=Rs, mean_indegree=mdeg)

# build a lookup Dict keyed by (R95, σ)
keyys = Tuple{Int,Float64}[(res2.R95[i], res2.σ[i]) for i in 1:nrow(res2)]
lut  = Dict( keyys[i] => res2.ΔAmean[i] for i in 1:nrow(res2) )

xs = sort(unique(res2.R95))
ys = sort(unique(res2.σ))
Z  = [get(lut, (x,y), NaN) for x in xs, y in ys]   # no silent fallbacks

# optional: assert no NaNs
@assert all(!isnan, Z) "Missing (R95, σ) combos in res2"