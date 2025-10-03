# ==========================================================
# Run baseline AM vs BAM (no habitat loss / no movement)
# Produces CSVs and PNGs in data/ and data/figs/
# Plots use Makie; each plot is wrapped in begin...display(fig) blocks
# ==========================================================
include("../SetUp.jl");
include("src/metaweb.jl");   using .MetaWeb
# include("src/metaweb_old_approach.jl");   using .MetaWeb
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
S             = 200
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

Cs     = range(0.001, 0.1; length=15)
Aligns = range(0.0, 1.0; length=15)
R95s   = Int.(range(4.0, 40.0; length=10))
Sigmas = range(0.02, 0.3; length=15)

# climate grid (choose gradient here)
grid_type = "fractal"
Cgrid = Climate.make_climate_grid(nx, ny; kind=Symbol(grid_type), seed=11)

# ---------------------------
# helpers
# ---------------------------
function run_once(rng; Cgrid, align, σ, R95,
                  motif_mix=:mixed, S=175, basal_frac=0.3,
                  τA=0.5, kreq=1, connectance=0.10)

    # decide which is free by looking at what the caller passed:
    control = (:R95 => R95, :C => connectance)
    mode = (!isnothing(R95) && (isnothing(connectance))) ? :R95 : :C

    mw = MetaWeb.build_metaweb(rng; S=S, basal_frac=basal_frac,
                               control=mode,
                               connectance=(mode==:C ? connectance : nothing),
                               R95=(mode==:R95 ? R95 : nothing),
                               motif_mix=motif_mix)

    μ, σi = Niches.make_niches(rng, S; align=align, σ=σ, basal_frac=basal_frac)
    pars  = BAM.Params(; τA=τA, kreq=kreq)
    out   = BAM.compute_AM_BAM(rng, mw, Cgrid, μ, σi, pars)

    am = BAM.species_areas(out[:AM_maps])
    bm = BAM.species_areas(out[:BAM_maps])

    is_cons  = mw.trophic_role .!= :basal
    vπ       = filter(isfinite, am[is_cons])
    piA_cons = isempty(vπ) ? NaN : mean(vπ)

    Psuff = BAM.prey_sufficiency(mw, out[:Akeep]; kreq=kreq)
    Ks    = BAM.K_spectrum(mw, out[:Akeep])
    vK    = filter(isfinite, Ks)
    Kp90  = isempty(vK) ? NaN : quantile(vK, 0.9)

    ΔA    = Metrics.delta_area(am, bm)
    ΔG    = Metrics.delta_gini(am, bm)

    # realized stats (useful for sanity)
    Creal = MetaWeb.global_connectance(mw.A)
    nb    = round(Int, basal_frac*S)
    outdeg = sum(mw.A; dims=2)[:]
    R95real = quantile(outdeg[(nb+1):end], 0.95)

    return (ΔA=ΔA, ΔG=ΔG, Psuff=Psuff, Kmean=mean(Ks),
            Kp90=Kp90, πA=piA_cons, Creal=Creal, R95real=R95real)
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
                         τA=fixed.τA, kreq=fixed.kreq)
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
# Sweep 1: (C, align)  → control = :C (R95 free)
# ---------------------------
sweep1 = NamedTuple[(; C=c, align=a) for c in Cs for a in Aligns]

fixed1 = (; Cgrid=Cgrid,
          align=default_align, σ=default_σ,
          R95=default_R95,          # ignored when C is present in the sweep
          C=default_C,              # default used only if a point doesn't provide :C
          S=S, basal_frac=basal_frac, τA=τA, kreq=kreq)

# run + persist (toggle the if you want to reuse previous output)
res1 = replicate_sweep(rng, sweep1; fixed=fixed1, replicates=replicates)
CSV.write(joinpath(@__DIR__, "data", "sweep_C_align_$grid_type.csv"), res1)
# res1 = CSV.read(joinpath(@__DIR__, "data", "sweep_C_align_$grid_type.csv"), DataFrame)

# sanity: this sweep should have used control=:C everywhere
@assert all(res1.control .== :C) "Expected control=:C for (C, align) sweep."

# ---------------------------
# Heatmap of ΔArea(C, align)
# ---------------------------
begin
    # aggregate in case you later change replicate_sweep to return multiple rows per cell
    res1_clean = combine(groupby(res1, [:C, :align]),
                         :ΔAmean => mean => :ΔA)

    xs = sort(unique(res1_clean.C))
    ys = sort(unique(res1_clean.align))

    # build dense matrix safely
    Z = fill(NaN, length(xs), length(ys))
    let lut = Dict( (res1_clean.C[i], res1_clean.align[i]) => res1_clean.ΔA[i]
                    for i in 1:nrow(res1_clean) )
        for (ix, x) in pairs(xs), (iy, y) in pairs(ys)
            v = get(lut, (x, y), NaN)
            Z[ix, iy] = v
        end
    end
    @assert all(!isnan, Z) "Missing (C, align) combinations in res1."

    lo, hi = extrema(Z)

    fig = Figure(; size=(900,500))
    ax  = Axis(fig[1,1], xlabel="connectance C", ylabel="alignment",
               title="ΔArea (AM–BAM) — C vs align $grid_type")
    hm  = heatmap!(ax, xs, ys, permutedims(Z); colorrange=(lo,hi))
    Colorbar(fig[1,2], hm, label="ΔArea")
    display(fig)
    save(joinpath(@__DIR__, "data", "figs",
         "HM_DeltaArea_C_align_$grid_type.png"), fig)
end

# ---------------------------
# ΔArea vs Psuff scatter (mediation)
# ---------------------------
begin
    fig = Figure(; size=(800,500))
    ax  = Axis(fig[1,1], xlabel="P_suff (mean)", ylabel="ΔArea (mean)",
               title="ΔArea vs P_suff — (C, align) sweep $grid_type")
    scatter!(ax, res1.Pmean, res1.ΔAmean)
    display(fig)
    save(joinpath(@__DIR__, "data", "figs",
         "Scatter_DeltaArea_vs_Psuff_$grid_type.png"), fig)
end

# ---------------------------
# helpers for plotting
# ---------------------------
# Build a dense matrix Z[x,y] from tall results; value_col is a Symbol in res
function _pivot_for_heatmap(res::DataFrame, xcol::Symbol, ycol::Symbol; value_col::Symbol=:ΔAmean)
    # aggregate in case replicate_sweep ever returns multiple rows per (x,y)
    tbl = combine(groupby(res, [xcol, ycol]), value_col => mean => :val)

    xs = sort(unique(tbl[!, xcol]))
    ys = sort(unique(tbl[!, ycol]))

    Z = fill(NaN, length(xs), length(ys))
    lut = Dict( (tbl[i, xcol], tbl[i, ycol]) => tbl[i, :val] for i in 1:nrow(tbl) )
    for (ix, x) in pairs(xs), (iy, y) in pairs(ys)
        Z[ix, iy] = get(lut, (x, y), NaN)
    end
    @assert all(!isnan, Z) "Missing grid points for $(xcol)×$(ycol)."
    return xs, ys, Z
end

function _heatmap!(xs, ys, Z; xlabel::AbstractString, ylabel::AbstractString, title::AbstractString)
    lo, hi = extrema(Z)
    ax = Axis(current_figure()[1,1], xlabel=xlabel, ylabel=ylabel, title=title)
    hm = heatmap!(ax, xs, ys, permutedims(Z); colorrange=(lo, hi))
    Colorbar(current_figure()[1,2], hm, label="ΔArea")
    return ax
end

# ---------------------------
# Sweep 1.2: (C, σ)  → control=:C (R95 free)
# ---------------------------
sweep1_2 = NamedTuple[(; C=c, σ=s) for c in Cs for s in Sigmas]

fixed1_2 = (; Cgrid=Cgrid, align=default_align, σ=default_σ, R95=default_R95, # R95 ignored here
             C=default_C, S=S, basal_frac=basal_frac, τA=τA, kreq=kreq)

res1_2 = replicate_sweep(rng, sweep1_2; fixed=fixed1_2, replicates=replicates)
@assert all(res1_2.control .== :C) "Expected control=:C for (C, σ) sweep."

begin
    fig = Figure(size=(900,500))
    xs, ys, Z = _pivot_for_heatmap(res1_2, :C, :σ; value_col=:ΔAmean)
    _heatmap!(xs, ys, Z; xlabel="connectance C", ylabel="niche breadth σ",
              title="ΔArea (AM–BAM) — C vs σ $grid_type")
    display(fig)
    save(joinpath(@__DIR__, "data", "figs", "HM_DeltaArea_C_sigma_$grid_type.png"), fig)
end

# ---------------------------
# Sweep 2: (R95, σ) → control=:R95 (C free)
# ---------------------------
sweep2 = NamedTuple[(; R95=r, σ=s) for r in R95s for s in Sigmas]

fixed2 = (; Cgrid=Cgrid, align=default_align, σ=default_σ, R95=default_R95, # default only
           C=default_C, S=S, basal_frac=basal_frac, τA=τA, kreq=kreq)

res2 = replicate_sweep(rng, sweep2; fixed=fixed2, replicates=replicates)
@assert all(res2.control .== :R95) "Expected control=:R95 for (R95, σ) sweep."

begin
    fig = Figure(size=(980,520))
    xs, ys, Z = _pivot_for_heatmap(res2, :R95, :σ; value_col=:ΔAmean)
    _heatmap!(xs, ys, Z; xlabel="R95 (95th out-degree)", ylabel="niche breadth σ",
              title="ΔArea (AM–BAM) — R95 vs σ $grid_type")
    display(fig)
    save(joinpath(@__DIR__, "data", "figs", "HM_DeltaArea_R_sigma_$grid_type.png"), fig)
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
    save(joinpath(@__DIR__, "data", "figs", "Bars_DeltaGini_corners_$grid_type.png"), fig)
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
    save(joinpath(@__DIR__, "data", "figs", "Mediation_curves_$grid_type.png"), fig)
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