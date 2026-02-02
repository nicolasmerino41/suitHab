#!/usr/bin/env julia
# Julia 1.11
# Run:  julia --threads auto fig1_sweep_only.jl
#
# PURPOSE
# - This script reuses the FUNCTIONS from your previous “FINAL” script
#   (traffic_extinction_figures_FINAL.jl) and ONLY produces “Figure 1” sweeps.
#
# OUTPUT
# - 5 figures, each a 2×3 panel:
#     rows = 2 divergence regimes (ordered from highest → lowest divergence)
#     cols = 3 habitat-loss geometries (random / cluster / front)
# - Total regimes = 10 (two per figure).
# - Each row title contains the exact parametrisation used.
#
# REQUIREMENT
# - Your FINAL script must NOT auto-run when included. Apply the patch above.

using Random
using Statistics
using CairoMakie
using Printf
using Dates

# ---------------------------------------------
# 1) Include the function library (your FINAL)
# ---------------------------------------------
# Adjust filename if yours differs.
# include("traffic_extinction_figures_FINAL.jl")

# ---------------------------------------------
# 2) User-facing knobs for this sweep
# ---------------------------------------------
const WHICH_FIG1 = :consumers     # change to :all if you want the all-species version
const NREP = NREP_FIG1            # reuse constant from included script (or set your own)
const DO_PRINT_DIAG_FIRST = true  # prints diagnostics only for case 1, random geometry, rep 1

ts = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
OUTDIR_SWEEP = joinpath(pwd(), "fig1_sweep_" * ts)
isdir(OUTDIR_SWEEP) || mkpath(OUTDIR_SWEEP)

println("Saving sweep figures to: ", OUTDIR_SWEEP)
println("WHICH_FIG1 = $(WHICH_FIG1) ; NREP = $(NREP)")

# ---------------------------------------------
# 3) Define the 10 regimes (highest → lowest divergence)
# ---------------------------------------------
# “Divergence ladder”:
# - Connectance increases (more redundancy → less divergence)
# - Correlation increases (consumer niches align with prey → less divergence)
# - Niche scenario broadens (SkewNarrow → VarOcc → SkewBroad)
# - prey_sigma_factor increases to 1 (prey less bottlenecked)
#
# NOTE: these are intentionally “article illustrative”, not a statistical fit.
#
struct Fig1Regime
    C::Float64
    r::Float64
    scen::NicheScenario
    prey_sigma_factor::Float64
end

function scen_name(s::NicheScenario)
    s isa SkewNarrow      && return "SkewNarrow"
    s isa VarOccVarSigma  && return "VarOcc+VarSigma"
    s isa VarOccMildSigma && return "VarOcc+MildSigma"
    s isa SkewBroad       && return "SkewBroad"
    return "Unknown"
end

const REGIMES = Fig1Regime[
    Fig1Regime(0.02, -0.50, SkewNarrow(),      0.65),
    Fig1Regime(0.03, -0.35, SkewNarrow(),      0.70),
    Fig1Regime(0.04, -0.20, VarOccVarSigma(),  0.75),
    Fig1Regime(0.05, -0.05, VarOccVarSigma(),  0.80),
    Fig1Regime(0.07,  0.10, VarOccVarSigma(),  0.85),
    Fig1Regime(0.09,  0.25, VarOccMildSigma(), 0.90),
    Fig1Regime(0.12,  0.40, VarOccMildSigma(), 0.95),
    Fig1Regime(0.15,  0.55, SkewBroad(),       1.00),
    Fig1Regime(0.18,  0.75, SkewBroad(),       1.00),
    Fig1Regime(0.20,  0.90, SkewBroad(),       1.00),
]

# ---------------------------------------------
# 4) Fig1 replicate simulation, parameterised
# ---------------------------------------------
function simulate_replicate_fig1_param!(
    rng::AbstractRNG,
    ws::CCWorkspace,
    geometry::Symbol,
    regime::Fig1Regime,
    which::Symbol;
    do_print_diag::Bool=false
)
    nb = round(Int, BASAL_FRAC * S)
    basal_mask = BitVector(falses(S))
    basal_mask[1:nb] .= true

    env1, env2 = make_climate(rng)
    prey = make_metaweb(rng, basal_mask, regime.C)
    mu1, mu2, rgot = assign_centroids_with_target_corr(rng, prey, basal_mask, regime.r)

    # sigmas (occupancy scenario + viability), prey narrower when prey_sigma_factor < 1
    s1, s2, clim, frac_extA0_all = assign_sigmas_with_viability(
        rng, ws, env1, env2, mu1, mu2, basal_mask, regime.scen;
        prey_sigma_factor=regime.prey_sigma_factor
    )

    # soften globally if baseline extant fraction under A-only is too low
    if frac_extA0_all < BASELINE_EXTANT_TARGET
        factor = (BASELINE_EXTANT_TARGET / max(frac_extA0_all, 1e-6))^(0.35)
        for sp in 1:S
            s1[sp] *= factor
            s2[sp] *= factor
            clim[sp] = suitability_mask(env1, env2, mu1[sp], mu2[sp], s1[sp], s2[sp], SUIT_THRESH)
        end
    end

    order = geometry == :random  ? order_random(rng) :
            geometry == :cluster ? order_clustered(rng) :
            geometry == :front   ? order_front(rng) :
            error("Unknown geometry")

    habitat0 = BitVector(trues(NCELLS))
    A0 = [clim[sp] .& habitat0 for sp in 1:S]

    # AB fixed-point baseline for SAR_eff (non-circular)
    AB0_fp = fixed_point_AB(A0, prey, basal_mask)

    # baseline richness for anchoring SAR curves
    extA0 = extant_mask!(ws, A0, which, basal_mask)
    presAB0_casc, extAB0_all = extinction_cascade_AB!(rng, ws, A0, prey, basal_mask)

    extAB0 = if which == :all
        extAB0_all
    else
        tmp = copy(extAB0_all)
        @inbounds for sp in 1:S
            if basal_mask[sp]; tmp[sp] = false; end
        end
        tmp
    end

    S0A  = float(gamma_richness_from_ext(extA0, which, basal_mask))
    S0AB = float(gamma_richness_from_ext(extAB0, which, basal_mask))

    # SAR fits (area-only)
    Apts,  Spts  = sar_points_gamma(A0, which, basal_mask)
    z_sar = fit_z_only(Apts, Spts)

    Apts2, Spts2 = sar_points_gamma(AB0_fp, which, basal_mask)
    z_eff = fit_z_only(Apts2, Spts2)

    # effective supported area baseline (geometry sensitive but non-circular)
    supp0 = supported_union_mask(AB0_fp, which, basal_mask)
    Aeff0 = effective_area(ws, supp0)
    A0_total = float(NCELLS)

    # curves
    nH = length(HL)
    gamma_A  = zeros(Float64, nH)
    gamma_AB = zeros(Float64, nH)
    sar_baseline  = zeros(Float64, nH)
    sar_effective = zeros(Float64, nH)

    for (k,h) in enumerate(HL)
        kdestroy = round(Int, h * NCELLS)
        habitat = habitat_mask_from_order(order, kdestroy)
        A_pres = [clim[sp] .& habitat for sp in 1:S]

        extA  = extant_mask!(ws, A_pres, which, basal_mask)
        gamma_A[k] = gamma_richness_from_ext(extA, which, basal_mask)

        _, extAB_all = extinction_cascade_AB!(rng, ws, A_pres, prey, basal_mask)
        extAB = if which == :all
            extAB_all
        else
            tmp = copy(extAB_all)
            @inbounds for sp in 1:S
                if basal_mask[sp]; tmp[sp] = false; end
            end
            tmp
        end
        gamma_AB[k] = gamma_richness_from_ext(extAB, which, basal_mask)

        # SAR baseline on remaining area
        Arem = (1.0 - h) * A0_total
        sar_baseline[k] = sar_predict(S0A, A0_total, z_sar, Arem)

        # SAR effective on effective supported area (AB fixed-point, no extinction)
        AB_fp = fixed_point_AB(A_pres, prey, basal_mask)
        supp  = supported_union_mask(AB_fp, which, basal_mask)
        Aeff  = effective_area(ws, supp)
        sar_effective[k] = sar_predict(S0AB, max(1.0, Aeff0), z_eff, max(1.0, Aeff))
    end

    if ENFORCE_MONOTONE
        enforce_monotone_nonincreasing!(gamma_A)
        enforce_monotone_nonincreasing!(gamma_AB)
        enforce_monotone_nonincreasing!(sar_baseline)
        enforce_monotone_nonincreasing!(sar_effective)
    end

    if do_print_diag
        md, sd, z0 = prey_stats(prey, basal_mask)
        println("---- DIAGNOSTICS (Sweep replicate) ----")
        println("Geometry=$(geometry)  WHICH=$(which)")
        println("C=$(regime.C)  r=$(regime.r)  scen=$(scen_name(regime.scen))  prey_sigma_factor=$(regime.prey_sigma_factor)")
        println("Achieved centroid corr (approx): $(round(rgot,digits=3))")
        println("Consumers prey-degree mean=$(round(md,digits=2)) sd=$(round(sd,digits=2)), zero-prey consumers=$(z0)")
        println("---------------------------------------")
    end

    return (gamma_A=gamma_A, gamma_AB=gamma_AB, sar_baseline=sar_baseline, sar_effective=sar_effective)
end

# ---------------------------------------------
# 5) Mean curves over replicates (threaded)
# ---------------------------------------------
function mean_curves_param(curve_list::Vector{NamedTuple})
    out = Dict{Symbol,Vector{Float64}}()
    for k in (:gamma_A, :gamma_AB, :sar_baseline, :sar_effective)
        mats = reduce(hcat, [getfield(c,k) for c in curve_list])
        out[k] = vec(mean(mats; dims=2))
        ENFORCE_MONOTONE && enforce_monotone_nonincreasing!(out[k])
    end
    return out
end

function compute_case_curves(case_idx::Int, geom::Symbol, which::Symbol)
    reps = Vector{NamedTuple}(undef, NREP)
    regime = REGIMES[case_idx]

    Threads.@threads for r in 1:NREP
        tid = Threads.threadid()
        rng = rngs[tid]
        ws  = wss[tid]

        seed = BASE_SEED +
               1_000_000 * case_idx +
               10_000    * (findfirst(==(geom), [:random,:cluster,:front])) +
               100       * (which == :consumers ? 2 : 1) +
               r
        Random.seed!(rng, seed)

        reps[r] = simulate_replicate_fig1_param!(
            rng, ws, geom, regime, which;
            do_print_diag = (DO_PRINT_DIAG_FIRST && case_idx==1 && geom==:random && r==1)
        )
    end

    return mean_curves_param(reps)
end

# ---------------------------------------------
# 6) Plotting: 2×3 for two regimes
# ---------------------------------------------
function regime_label(k::Int, reg::Fig1Regime)
    return @sprintf(
        "Case %02d: C=%.3f, r=%.2f, scen=%s, preyσ×=%.2f",
        k, reg.C, reg.r, scen_name(reg.scen), reg.prey_sigma_factor
    )
end

function plot_two_case_panel(case_a::Int, case_b::Int; which::Symbol)
    f = Figure(size=(1650, 900))
    geoms = [:random, :cluster, :front]
    geom_titles = Dict(:random=>"Random loss", :cluster=>"Clustered loss", :front=>"Front loss")

    # Big title
    Label(f[0, :],
        "Figure 1 sweep — $(which == :consumers ? "consumers only" : "all species") — Cases $(case_a)–$(case_b)",
        fontsize=20
    )

    # Row labels
    Label(f[1, 0], regime_label(case_a, REGIMES[case_a]), rotation=pi/2, fontsize=14)
    Label(f[2, 0], regime_label(case_b, REGIMES[case_b]), rotation=pi/2, fontsize=14)

    # Make room for row labels column
    colsize!(f.layout, 0, Relative(0.06))

    for (row, case_idx) in enumerate([case_a, case_b])
        for (col, g) in enumerate(geoms)
            ax = Axis(f[row, col],
                title = geom_titles[g],
                xlabel = "Habitat loss (proportion destroyed)",
                ylabel = "Gamma richness ($(which == :consumers ? "consumers" : "all"))"
            )
            curves = compute_case_curves(case_idx, g, which)
            add_lines!(ax, HL, curves[:gamma_A], curves[:gamma_AB], curves[:sar_baseline], curves[:sar_effective])
        end
    end

    return f
end

# ---------------------------------------------
# 7) Run: produce 5 figures (10 regimes)
# ---------------------------------------------
case_pairs = [(1,2), (3,4), (5,6), (7,8), (9,10)]

for (a,b) in case_pairs
    println("Computing sweep panel for cases $(a)-$(b) ...")
    fig = plot_two_case_panel(a, b; which=WHICH_FIG1)
    fn = joinpath(OUTDIR_SWEEP, @sprintf("fig1_sweep_cases_%02d_%02d_%s.png", a, b, Symbol(WHICH_FIG1)))
    save(fn, fig)
    display(fig)
end

println("\nDone. Sweep output directory:")
println(OUTDIR_SWEEP)
