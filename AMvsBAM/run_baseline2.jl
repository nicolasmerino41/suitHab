############## run_baseline2.jl — clean recording of free variables ##############

include("../SetUp.jl")
include("src/metaweb.jl");   using .MetaWeb
include("src/niches.jl");    using .Niches
include("src/bam.jl");       using .BAM
include("src/metrics.jl");   using .Metrics
include("src/climate.jl");   using .Climate

const Mke = CairoMakie

# ---------------------------
# CONFIG
# ---------------------------
mkpath(joinpath(@__DIR__, "data", "UNI"))
mkpath(joinpath(@__DIR__, "data", "UNI", "figs"))

rng           = MersenneTwister(1)
S             = 200
basal_frac    = 0.30
nx, ny        = 40, 40
τA            = 0.5
kreq          = 1
replicates    = 10

# Defaults for slices
default_R95   = 5
default_C     = 0.10
default_σ     = 0.4
default_align = 0.4

Cs     = range(0.001, 0.10; length=16)
Sigmas = range(0.02, 0.30; length=16)

# climate grid
grid_type = "gradient"   # try "gradient" or "fractal" too
Cgrid = Climate.make_climate_grid(nx, ny; kind=Symbol(grid_type), seed=11)

# ---------------------------
# CORE: one run
# ---------------------------
"""
run_once(; ...)

Decides the control regime automatically:
- if R95 is `nothing` and C is provided  → control=:C (C fixed, R95 free)
- if C is `nothing` and R95 is provided  → control=:R95 (R95 fixed, C free)

Returns ΔA, ΔG, Psuff, K stats, πA plus Creal & R95real and the control used.
"""
function run_once(rng; Cgrid, align, σ, R95,
                  motif_mix=:mixed, S=175, basal_frac=0.3,
                  τA=0.5, kreq=1, connectance=0.10)

    # decide which is free
    mode = (isnothing(connectance) && !isnothing(R95)) ? :R95 : :C

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

    # consumers only
    is_cons   = mw.trophic_role .!= :basal
    vπ        = filter(isfinite, am[is_cons])
    piA_cons  = isempty(vπ) ? NaN : mean(vπ)

    Psuff = BAM.prey_sufficiency(mw, out[:Akeep]; kreq=kreq)

    Ks    = BAM.K_spectrum(mw, out[:Akeep])
    vK    = filter(isfinite, Ks)
    Kp90  = isempty(vK) ? NaN : quantile(vK, 0.9)

    ΔA    = Metrics.delta_area(am, bm)
    ΔG    = Metrics.delta_gini(am, bm)

    # realized structural stats
    Creal = MetaWeb.global_connectance(mw.A)
    nb    = round(Int, basal_frac*S)
    outdeg = sum(mw.A; dims=2)[:]
    R95real = quantile(outdeg[(nb+1):end], 0.95)

    return (ΔA=ΔA, ΔG=ΔG, Psuff=Psuff, Kmean=mean(Ks),
            Kp90=Kp90, πA=piA_cons, Creal=Creal, R95real=R95real, control=mode)
end

# ---------------------------
# sweep helper
# ---------------------------
"""
replicate_sweep(rng, sweep; fixed, replicates)

- Automatically routes each grid point to the right control regime based on
  whether the point carries :C or :R95.
- Aggregates lo/mean/hi for all outcome metrics.
- Also aggregates the **realized** free variables:
    - for control=:C, reports R95real_lo/mean/hi and Creal_mean (sanity)
    - for control=:R95, reports Creal_lo/mean/hi and R95real_mean (sanity)
- Returns a DataFrame with one row per point in `sweep`.
"""
function replicate_sweep(rng, sweep; fixed::NamedTuple, replicates::Int=20)
    n = length(sweep)
    results_threads = [NamedTuple[] for _ in 1:Threads.nthreads()]

    Threads.@threads for idx in 1:n
        tid  = Threads.threadid()
        pars = sweep[idx]

        # independent RNG per thread, re-seeded for reproducibility
        local_rng = MersenneTwister(rand(Int))

        # detect control regime for this point:
        want_C   = haskey(pars, :C)   && !haskey(pars, :R95)
        want_R95 = haskey(pars, :R95) && !haskey(pars, :C)
        mode = want_R95 ? :R95 : :C

        vals = [ run_once(MersenneTwister(rand(Int));
                  Cgrid=fixed.Cgrid,
                  align=get(pars, :align, fixed.align),
                  σ=get(pars, :σ, fixed.σ),
                  R95=(mode==:R95 ? get(pars, :R95, fixed.R95) : nothing),
                  connectance=(mode==:C ? get(pars, :C, fixed.C) : nothing),
                  motif_mix=get(pars, :motif, :mixed),
                  S=fixed.S, basal_frac=fixed.basal_frac,
                  τA=fixed.τA, kreq=fixed.kreq)
        for _ in 1:replicates ]

        # helpers
        band(x) = Metrics.qband(x)  # returns (lo, mean, hi)

        ΔAlo, ΔAmean, ΔAhi = band([v.ΔA     for v in vals])
        ΔGlo, ΔGmean, ΔGhi = band([v.ΔG     for v in vals])
        Plo,  Pmean,  Phi  = band([v.Psuff  for v in vals])
        Klo,  Kmean,  Khi  = band([v.Kp90   for v in vals])
        πAlo, πAmean, πAhi = band([v.πA     for v in vals])

        Creals  = [v.Creal   for v in vals]
        R95reals= [v.R95real for v in vals]

        if mode == :C
            Rlo, Rmean, Rhi = band(R95reals)
            # keep Creal mean for sanity (should be ~target C)
            row = merge(pars, (; control=mode,
                                 ΔAlo, ΔAmean, ΔAhi,
                                 ΔGlo, ΔGmean, ΔGhi,
                                 Plo, Pmean, Phi,
                                 Klo, Kmean, Khi,
                                 πAlo, πAmean, πAhi,
                                 Creal_mean = mean(Creals),
                                 R95real_lo = Rlo, R95real_mean = Rmean, R95real_hi = Rhi))
            push!(results_threads[tid], row)

        else # mode == :R95
            Clo, Cmean, Chi = band(Creals)
            row = merge(pars, (; control=mode,
                                 ΔAlo, ΔAmean, ΔAhi,
                                 ΔGlo, ΔGmean, ΔGhi,
                                 Plo, Pmean, Phi,
                                 Klo, Kmean, Khi,
                                 πAlo, πAmean, πAhi,
                                 R95real_mean = mean(R95reals),
                                 Creal_lo = Clo, Creal_mean = Cmean, Creal_hi = Chi))
            push!(results_threads[tid], row)
        end
    end

    flat = reduce(vcat, results_threads)
    return DataFrame(flat)
end

# ---------------------------
# plotting helpers
# ---------------------------
# Dense heatmap (ΔAmean by x,y)
function heatmap_from_results(res::DataFrame, xcol::Symbol, ycol::Symbol; title::AbstractString)
    tbl = combine(groupby(res, [xcol, ycol]), :ΔAmean => mean => :val)
    xs = sort(unique(tbl[!, xcol]))
    ys = sort(unique(tbl[!, ycol]))
    Z = fill(NaN, length(xs), length(ys))
    lut = Dict( (tbl[i, xcol], tbl[i, ycol]) => tbl[i, :val] for i in 1:nrow(tbl) )
    for (ix,x) in pairs(xs), (iy,y) in pairs(ys)
        Z[ix,iy] = get(lut, (x,y), NaN)
    end
    @assert all(!isnan, Z) "Missing points in $(xcol)×$(ycol)."

    lo, hi = extrema(Z)
    fig = Figure(; size=(960,520))
    ax  = Axis(fig[1,1], xlabel=string(xcol), ylabel=string(ycol), title=title)
    hm  = heatmap!(ax, xs, ys, permutedims(Z); colorrange=(lo,hi))
    Colorbar(fig[1,2], hm, label="ΔArea")
    # display(fig)
    return fig
end

# Univariate ΔA vs swept var, colored by the **realized free variable**
function plot_univariate_with_free(res::DataFrame, xcol::Symbol)
    @assert :control ∈ propertynames(res)
    ctrl = first(res.control)
    free_sym = (ctrl == :C) ? :R95real_mean : :Creal_mean
    fig = Figure(; size=(840,480))
    ax  = Axis(fig[1,1], xlabel=string(xcol), ylabel="ΔArea (mean)",
               title="ΔArea vs $(xcol) (colored by $(free_sym)) $(grid_type)")
    sc  = scatter!(ax, res[!, xcol], res.ΔAmean; color=res[!, free_sym], colormap=:viridis)
    Colorbar(fig[1,2], sc, label=string(free_sym))
    # display(fig)
    return fig
end

# ---------------------------
# EXAMPLE SWEEP A: (C, σ)  → control=:C, R95 free
# ---------------------------
sweep_Cσ = NamedTuple[(; C=c, σ=s) for c in Cs for s in Sigmas]
fixed_Cσ = (; Cgrid=Cgrid, align=default_align, σ=default_σ, R95=default_R95,
            C=default_C, S=S, basal_frac=basal_frac, τA=τA, kreq=kreq)

res_Cσ = replicate_sweep(rng, sweep_Cσ; fixed=fixed_Cσ, replicates=replicates)
@assert all(res_Cσ.control .== :C)

# Heatmap ΔA(C, σ)
fig_HM_Cσ = heatmap_from_results(res_Cσ, :C, :σ; title="ΔArea (AM–BAM) — C vs σ $(grid_type)")

# Univariate: ΔA vs C, colored by realized R95
# (use the mid-σ slice to keep it 1D or just show all points)
# res_midσ = filter(row -> isapprox(row.σ, median(Sigmas); atol=1e-8), res_Cσ)
fig_uni_C = plot_univariate_with_free(res_Cσ, :C)
save(joinpath(@__DIR__, "data/UNI", "figs", "UNI_Delta_vs_C_col_by_R95_$grid_type.png"), fig_uni_C)

save(joinpath(@__DIR__, "data/UNI", "figs", "HM_DeltaArea_C_sigma_$grid_type.png"), fig_HM_Cσ)

# ---------------------------
# EXAMPLE SWEEP B: (R95, σ) → control=:R95, C free
# ---------------------------
R95s   = Int.(range(4.0, 40.0; length=10))
sweep_Rσ = NamedTuple[(; R95=r, σ=s) for r in R95s for s in Sigmas]
fixed_Rσ = (; Cgrid=Cgrid, align=default_align, σ=default_σ, R95=default_R95,
            C=default_C, S=S, basal_frac=basal_frac, τA=τA, kreq=kreq)

res_Rσ = replicate_sweep(rng, sweep_Rσ; fixed=fixed_Rσ, replicates=replicates)
@assert all(res_Rσ.control .== :R95)

# Heatmap ΔA(R95, σ)
fig_HM_Rσ = heatmap_from_results(res_Rσ, :R95, :σ;
    title="ΔArea (AM–BAM) — R95 vs σ $(grid_type)")
save(joinpath(@__DIR__, "data/UNI", "figs", "HM_DeltaArea_R_sigma_$grid_type.png"), fig_HM_Rσ)

# Univariate: ΔA vs R95, colored by realized Creal (free)
# res_midσ2 = filter(row -> isapprox(row.σ, median(Sigmas); atol=1e-8), res_Rσ)
fig_uni_R = plot_univariate_with_free(res_Rσ, :R95)
save(joinpath(@__DIR__, "data/UNI", "figs", "UNI_Delta_vs_R95_col_by_C_$grid_type.png"), fig_uni_R)

###############################################################################
# ---------------------------
# SWEEP C: (C, align) → control=:C, R95 free
# ---------------------------
Aligns = range(0.0, 1.0; length=16)   # local range for this sweep

sweep_C_align = NamedTuple[(; C=c, align=a) for c in Cs for a in Aligns]

fixed_C_align = (;
    Cgrid      = Cgrid,
    align      = default_align,
    σ          = default_σ,
    R95        = default_R95,   # (ignored; R95 is free when control=:C)
    C          = default_C,
    S          = S,
    basal_frac = basal_frac,
    τA         = τA,
    kreq       = kreq
)

res_C_align = replicate_sweep(rng, sweep_C_align; fixed=fixed_C_align, replicates=replicates)
@assert all(res_C_align.control .== :C)

# Heatmap ΔArea(C, align)
fig_HM_C_align = heatmap_from_results(res_C_align, :C, :align;
    title="ΔArea (AM–BAM) — C vs align $(grid_type)")
save(joinpath(@__DIR__, "data/UNI/figs", "HM_DeltaArea_C_align_$(grid_type).png"), fig_HM_C_align)

