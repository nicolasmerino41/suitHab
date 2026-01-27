###############################################################################
# patchonly_final_figures.jl
#
# Goal (copy-paste + run):
#   Produce ARTICLE-READY figures from your PATCH-only outputs:
#     Fig1 (main): 6 line plots (2 example cases × 3 HL geometries) under Emin_patch
#                  + (optional) same figure under viability for SI.
#     Fig2 (main): parameter sweep heatmaps (corr × connectance-or-k_prey),
#                  shown for multiple niche scenarios, using ΔAUC (raw + relative).
#     Fig3 (test): mechanism panel "support amount vs support connectivity"
#
# Also:
#   - Reports AUC (raw) and relative ΔAUC
#   - Enforces a monotonic sanity check on richness-vs-HL curves (optional correction)
#
# Expected inputs (in RESULTS_DIR):
#   - AUC summary CSV: auc_summary.csv
#       with at least: scenario_name, k_prey, corr, geometry, auc_emin, auc_viability
#       (or similarly named columns)
#   - Curve CSV (needed for Fig1 and to recompute AUC from curves):
#       One CSV containing richness curves vs habitat-loss f, for both models A-only and AB
#       and for metrics Emin_patch and viability.
#
# This script is defensive: it will try to autodetect the curve file + column names.
###############################################################################

using CSV, DataFrames, Statistics, Printf
using CairoMakie
using Glob

# ----------------------------- USER CONFIG -----------------------------------
const RESULTS_DIR = joinpath(pwd(), "Figures", "patchonly_36")  # change if needed
const OUT_DIR     = joinpath(pwd(), "Figures", "paper_final")
const SAVE_FIGS   = true
const MAKE_SI_VIABILITY = true

# Fig1: pick two example cases (one low divergence, one high divergence)
# If your data differs, just edit these selectors.
const EXAMPLE_LOW  = (scenario_name = "Skewed to broader niches (higher mean occupancy), still variable",
                      k_prey = 9, corr = 1.0)
const EXAMPLE_HIGH = (scenario_name = "Skewed to narrower niches (lower mean occupancy), still variable",
                      k_prey = 1, corr = 0.0)

# Fig2: which geometry to show as “main” in heatmaps (cluster is usually strongest)
const HEATMAP_GEOMETRY = "cluster"   # "random", "cluster", "front"
# Which scenarios to show in Fig2 (up to 4 recommended for a clean main figure)
# If empty, script will auto-pick up to 4 most frequent scenario_names.
const FIG2_SCENARIOS = String[]  # e.g. ["Very variable occupancy + very variable sigma", "Same occupancy range, less extreme sigma variability", ...]

# Monotonic sanity check options (applies to curves if present)
const ENFORCE_MONOTONE = true
const MONOTONE_FIX     = true   # if true: isotonic-like fix via cumulative minimum

# ----------------------------- STYLE -----------------------------------------
set_theme!(Theme(
    fontsize = 15,
    Axis = (xlabelsize = 15, ylabelsize = 15, titlesize = 16,
            xticklabelsize = 13, yticklabelsize = 13),
    Legend = (labelsize = 13,),
))

mkpath(OUT_DIR)

# ----------------------------- UTILITIES -------------------------------------
# robust float compare for selecting corr
isapprox01(a, b; atol=1e-9) = abs(a - b) ≤ atol

function try_read_csv(path::String)
    @info "Reading CSV: $path"
    return CSV.read(path, DataFrame)
end

function find_first_csv_with_columns(dir::String, required::Vector{Symbol})
    cands = sort(glob("*.csv", dir))
    for f in cands
        df = try
            CSV.read(f, DataFrame; limit=5)
        catch
            continue
        end
        cols = Set(Symbol.(names(df)))
        if all(r -> r in cols, required)
            return f
        end
    end
    return nothing
end

# trapezoid AUC assuming x sorted
function auc_trapz(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    n == length(y) || error("auc_trapz: x and y length mismatch")
    n < 2 && return 0.0
    s = 0.0
    @inbounds for i in 1:(n-1)
        dx = float(x[i+1] - x[i])
        s += 0.5 * dx * float(y[i] + y[i+1])
    end
    return s
end

# Monotone decreasing sanity check + optional fix
function monotone_report_and_fix!(df::DataFrame; group_cols::Vector{Symbol}, xcol::Symbol, ycol::Symbol)
    # returns: (n_groups, n_viol_groups, viol_fraction, max_increase)
    gdf = groupby(df, group_cols)
    n_groups = length(gdf)
    n_viol = 0
    max_inc = 0.0
    for sub in gdf
        sort!(sub, xcol)
        y = Vector{Float64}(sub[!, ycol])
        inc = maximum(diff(y); init=0.0)
        if inc > 1e-12
            n_viol += 1
            max_inc = max(max_inc, inc)
            if ENFORCE_MONOTONE && MONOTONE_FIX
                # cumulative minimum from left to right => non-increasing (since HL increases)
                yfix = similar(y)
                yfix[1] = y[1]
                for i in 2:length(y)
                    yfix[i] = min(yfix[i-1], y[i])
                end
                sub[!, ycol] .= yfix
            end
        end
    end
    viol_fraction = n_groups == 0 ? 0.0 : n_viol / n_groups
    return (n_groups, n_viol, viol_fraction, max_inc)
end

# ----------------------------- LOAD AUC --------------------------------------
function load_auc_summary(results_dir::String)
    # prefer exact name
    p = joinpath(results_dir, "auc_summary.csv")
    if isfile(p)
        df = try_read_csv(p)
    else
        # try to find something that looks like it
        p2 = find_first_csv_with_columns(results_dir, [:scenario_name, :geometry, :corr])
        p2 === nothing && error("Could not find auc_summary.csv (or any CSV with scenario_name/geometry/corr) in $results_dir")
        df = try_read_csv(p2)
    end

    # normalize types
    for c in [:scenario_name, :geometry]
        if c in names(df)
            df[!, c] = string.(df[!, c])
        end
    end
    if :corr in names(df)
        df[!, :corr] = Float64.(df[!, :corr])
    end
    if :k_prey in names(df)
        df[!, :k_prey] = Int.(round.(Float64.(df[!, :k_prey])))
    end

    # detect AUC columns
    # common: auc_emin, auc_viability
    if !(:auc_emin in names(df)) && (:aucE in names(df))
        rename!(df, :aucE => :auc_emin)
    end
    if !(:auc_viability in names(df)) && (:aucV in names(df))
        rename!(df, :aucV => :auc_viability)
    end
    for c in [:auc_emin, :auc_viability]
        if c in names(df)
            df[!, c] = Float64.(df[!, c])
        end
    end
    return df
end

# ----------------------------- LOAD CURVES -----------------------------------
"""
Attempt to load a curve CSV with richness vs habitat loss f.

Supported “long” formats:
  columns include: f, scenario_name, geometry, corr, k_prey, metric, model, richness
where:
  metric ∈ {"Emin_patch","viability"} (or variants),
  model ∈ {"A-only","AB"} (or variants).

Supported “wide” formats:
  columns include: f, scenario_name, geometry, corr, k_prey, metric, richness_Aonly, richness_AB
"""
function load_curve_data(results_dir::String)
    # Try common filenames first
    candidates = [
        joinpath(results_dir, "curve_summary.csv"),
        joinpath(results_dir, "curves_summary.csv"),
        joinpath(results_dir, "curves.csv"),
        joinpath(results_dir, "richness_curves.csv"),
        joinpath(results_dir, "patch_richness_curves.csv"),
        joinpath(results_dir, "timeseries.csv"),
    ]
    curve_path = nothing
    for p in candidates
        if isfile(p)
            curve_path = p
            break
        end
    end
    # fallback: pick first CSV that has an f column and at least one richness column
    if curve_path === nothing
        for p in sort(glob("*.csv", results_dir))
            dfh = try
                CSV.read(p, DataFrame; limit=5)
            catch
                continue
            end
            cols = Set(Symbol.(names(dfh)))
            if (:f in cols) && ( (:richness in cols) || (:patch_richness in cols) ||
                                 (:richness_Aonly in cols) || (:richness_AB in cols) )
                curve_path = p
                break
            end
        end
    end
    curve_path === nothing && return nothing

    df = try_read_csv(curve_path)

    # normalize
    if :scenario_name in names(df); df[!, :scenario_name] = string.(df[!, :scenario_name]) end
    if :geometry in names(df);      df[!, :geometry]      = string.(df[!, :geometry])      end
    if :metric in names(df);        df[!, :metric]        = string.(df[!, :metric])        end
    if :model in names(df);         df[!, :model]         = string.(df[!, :model])         end
    if :corr in names(df);          df[!, :corr]          = Float64.(df[!, :corr])         end
    if :k_prey in names(df);        df[!, :k_prey]        = Int.(round.(Float64.(df[!, :k_prey]))) end
    if :f in names(df);             df[!, :f]             = Float64.(df[!, :f])            end

    # unify richness column name
    if :patch_richness in names(df) && !(:richness in names(df))
        rename!(df, :patch_richness => :richness)
    end

    # If wide: convert to long
    if (:richness_Aonly in names(df)) && (:richness_AB in names(df)) && !(:model in names(df))
        dfa = select(df, Not([:richness_AB]))
        rename!(dfa, :richness_Aonly => :richness)
        dfa[!, :model] .= "A-only"
        dfb = select(df, Not([:richness_Aonly]))
        rename!(dfb, :richness_AB => :richness)
        dfb[!, :model] .= "AB"
        df = vcat(dfa, dfb; cols=:union)
    end

    # If model exists but different labels, normalize to "A-only"/"AB"
    if :model in names(df)
        df[!, :model] = replace.(df[!, :model],
            "Aonly"=>"A-only", "A_only"=>"A-only", "A-only"=>"A-only",
            "AB"=>"AB", "A+B"=>"AB"
        )
    end

    # Normalize metric naming (best-effort)
    if :metric in names(df)
        df[!, :metric] = replace.(df[!, :metric],
            "emin"=>"Emin_patch", "Emin"=>"Emin_patch", "Emin_patch"=>"Emin_patch",
            "viab"=>"viability", "Viability"=>"viability", "viability"=>"viability"
        )
    end

    # Validate minimal columns for plotting
    required = [:scenario_name, :geometry, :corr, :k_prey, :f, :model, :richness]
    missing = filter(c -> !(c in names(df)), required)
    !isempty(missing) && error("Curve file found ($curve_path) but missing columns: $(missing)\nNames: $(names(df))")

    df[!, :richness] = Float64.(df[!, :richness])

    return df
end

# --------------------- COMPUTE AUCs FROM CURVES ------------------------------
function compute_auc_table_from_curves(curves::DataFrame)
    # returns DataFrame with AUC_Aonly, AUC_AB, ΔAUC_raw, ΔAUC_rel for each metric
    hasmetric = (:metric in names(curves))
    if !hasmetric
        curves[!, :metric] .= "Emin_patch"
    end

    group_cols = [:scenario_name, :geometry, :corr, :k_prey, :metric, :model]
    # monotone check per (scenario,geometry,corr,k,metric,model)
    nG, nV, fracV, maxInc = monotone_report_and_fix!(curves;
        group_cols=group_cols, xcol=:f, ycol=:richness)

    @info @sprintf("Monotonic check: %d groups, %d with violations (%.1f%%). Max increase=%.4g. Fix applied=%s",
        nG, nV, 100*fracV, maxInc, (ENFORCE_MONOTONE && MONOTONE_FIX))

    g = groupby(curves, group_cols)
    tmp = DataFrame()
    tmp.scenario_name = String[]
    tmp.geometry = String[]
    tmp.corr = Float64[]
    tmp.k_prey = Int[]
    tmp.metric = String[]
    tmp.model = String[]
    tmp.auc = Float64[]
    for sub in g
        sort!(sub, :f)
        a = auc_trapz(sub.f, sub.richness)
        push!(tmp, (sub.scenario_name[1], sub.geometry[1], sub.corr[1], sub.k_prey[1], sub.metric[1], sub.model[1], a))
    end

    # pivot to A-only vs AB
    keycols = [:scenario_name, :geometry, :corr, :k_prey, :metric]
    aonly = filter(:model => ==("A-only"), tmp)
    ab    = filter(:model => ==("AB"),     tmp)
    rename!(aonly, :auc => :auc_Aonly)
    rename!(ab,    :auc => :auc_AB)
    select!(aonly, Not(:model))
    select!(ab,    Not(:model))

    df = outerjoin(aonly, ab; on=keycols, makeunique=true)

    # ΔAUC and relative ΔAUC (bounded 0..1-ish if AB ≤ A-only)
    df[!, :dAUC_raw] = df.auc_Aonly .- df.auc_AB
    denom = max.(df.auc_Aonly, 1e-9)
    df[!, :dAUC_rel] = df.dAUC_raw ./ denom

    return df
end

# ----------------------------- FIGURE 1 --------------------------------------
function make_fig1(curves::DataFrame; metric="Emin_patch",
                   ex_low=EXAMPLE_LOW, ex_high=EXAMPLE_HIGH,
                   outbase="Fig1_examples_Emin")
    geoms = ["random", "cluster", "front"]
    fig = Figure(size=(1400, 700))

    function plot_case!(row::Int, titleprefix::String, sel)
        for (j, g) in enumerate(geoms)
            ax = Axis(fig[row, j],
                title = "$(titleprefix) — $(g)",
                xlabel = (row == 2 ? "habitat loss f" : ""),
                ylabel = (j == 1 ? "patch richness" : "")
            )
            sub = filter(r ->
                r.metric == metric &&
                r.geometry == g &&
                r.scenario_name == sel.scenario_name &&
                r.k_prey == sel.k_prey &&
                isapprox01(r.corr, sel.corr),
                curves
            )
            if nrow(sub) == 0
                text!(ax, 0.5, 0.5, text="(no data)", align=(:center,:center))
                continue
            end
            for (lab, mod) in [("A-only","A-only"), ("AB","AB")]
                s2 = filter(:model => ==(mod), sub)
                sort!(s2, :f)
                lines!(ax, s2.f, s2.richness, label=lab, linewidth=3)
            end
            if row == 1 && j == 1
                axislegend(ax; position=:lb, framevisible=true)
            end
            xlims!(ax, 0, 1)
            ylims!(ax, 0, maximum(sub.richness) * 1.05)
        end
    end

    plot_case!(1, "Low divergence",  ex_low)
    plot_case!(2, "High divergence", ex_high)

    Label(fig[0, :], "PATCH-only — examples under $(metric)", fontsize=18)

    display(fig)
    if SAVE_FIGS
        save(joinpath(OUT_DIR, outbase * ".png"), fig, px_per_unit=2)
        save(joinpath(OUT_DIR, outbase * ".pdf"), fig)
        @info "Saved Fig1 to $(OUT_DIR)"
    end
    return fig
end

# ----------------------------- FIGURE 2 --------------------------------------
function pick_fig2_scenarios(df_auc::DataFrame)
    if !isempty(FIG2_SCENARIOS)
        return FIG2_SCENARIOS
    end
    # auto-pick up to 4 most frequent scenarios
    if :scenario_name ∉ names(df_auc)
        return String[]
    end
    counts = combine(groupby(df_auc, :scenario_name), nrow => :n)
    sort!(counts, :n, rev=true)
    return counts.scenario_name[1:min(4, nrow(counts))]
end

function make_heatmap_panel!(ax::Axis, d::DataFrame; xcol::Symbol, ycol::Symbol, zcol::Symbol,
                            title::String, xlabel::String, ylabel::String,
                            clim=nothing)
    # pivot to grid
    xs = sort(unique(d[!, xcol]))
    ys = sort(unique(d[!, ycol]))
    x_to_i = Dict(x => i for (i,x) in enumerate(xs))
    y_to_i = Dict(y => i for (i,y) in enumerate(ys))
    Z = fill(NaN, length(ys), length(xs))
    for r in eachrow(d)
        Z[y_to_i[r[ycol]], x_to_i[r[xcol]]] = r[zcol]
    end
    hm = heatmap!(ax, xs, ys, Z; colormap=:viridis, nan_color=:transparent,
                  colorrange = clim === nothing ? Makie.automatic : clim)
    ax.title = title
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    return hm
end

function make_fig2(aucdf::DataFrame; metric=:auc_emin, outbase="Fig2_sweep_Emin")
    # We want ΔAUC across A-only vs AB, but aucdf currently stores (A-only - AB) already? No:
    # Your auc_summary.csv is AUC(A-only - AB) under Emin and viability (i.e., ΔAUC already).
    # So here:
    #   raw ΔAUC = auc_emin (or auc_viability)
    #   relative = raw / (raw + something) is not possible without A-only baseline.
    #
    # Therefore: if your auc_summary.csv is already ΔAUC, we plot it as "raw ΔAUC".
    # If you compute from curves (preferred), pass that table instead (it has dAUC_raw and dAUC_rel).
    #
    has_rel = (:dAUC_rel in names(aucdf)) && (:dAUC_raw in names(aucdf))
    zraw = has_rel ? :dAUC_raw : metric
    zrel = has_rel ? :dAUC_rel : nothing

    scenarios = pick_fig2_scenarios(aucdf)
    isempty(scenarios) && error("No scenarios found for Fig2. Check AUC table columns.")

    # choose x-axis: realised connectance if present, else k_prey
    xcol = if (:connectance_real in names(aucdf))
        :connectance_real
    elseif (:connectance in names(aucdf))
        :connectance
    else
        :k_prey
    end

    # We'll show corr on x (horizontal) and connectance/k on y (vertical) to match your earlier plots
    # (swap if you prefer).
    # Here: x=corr, y=xcol (connectance or k_prey)
    xaxis = :corr
    yaxis = xcol

    # Main figure: 2x2 (or 1xN) panels across scenarios, fixed geometry
    n = length(scenarios)
    nrow = n <= 2 ? 1 : 2
    ncol = ceil(Int, n / nrow)

    # Raw ΔAUC figure
    fig_raw = Figure(size=(1500, 700))
    hms = Any[]
    for (i, sc) in enumerate(scenarios)
        r = (i-1) ÷ ncol + 1
        c = (i-1) % ncol + 1
        ax = Axis(fig_raw[r, c])
        sub = filter(row -> row.scenario_name == sc && row.geometry == HEATMAP_GEOMETRY, aucdf)
        if nrow(sub) == 0
            text!(ax, 0.5, 0.5, text="(no data)", align=(:center,:center))
            continue
        end
        hm = make_heatmap_panel!(ax, sub; xcol=xaxis, ycol=yaxis, zcol=zraw,
            title = "$(sc) — $(HEATMAP_GEOMETRY)",
            xlabel = "corr (consumer–prey niche correlation)",
            ylabel = (yaxis == :k_prey ? "k_prey (proxy for connectance)" : string(yaxis)),
        )
        push!(hms, hm)
    end
    if !isempty(hms)
        Colorbar(fig_raw[:, end+1], hms[1], label = has_rel ? "raw ΔAUC (A-only − AB)" : "ΔAUC (from auc_summary)")
    end
    Label(fig_raw[0, :], has_rel ? "Raw ΔAUC (A-only − AB) — $(HEATMAP_GEOMETRY)" : "ΔAUC (from auc_summary) — $(HEATMAP_GEOMETRY)", fontsize=18)

    display(fig_raw)
    if SAVE_FIGS
        save(joinpath(OUT_DIR, outbase * "_raw.png"), fig_raw, px_per_unit=2)
        save(joinpath(OUT_DIR, outbase * "_raw.pdf"), fig_raw)
    end

    # Relative ΔAUC figure (only if available)
    fig_rel = nothing
    if has_rel
        fig_rel = Figure(size=(1500, 700))
        hms2 = Any[]
        for (i, sc) in enumerate(scenarios)
            r = (i-1) ÷ ncol + 1
            c = (i-1) % ncol + 1
            ax = Axis(fig_rel[r, c])
            sub = filter(row -> row.scenario_name == sc && row.geometry == HEATMAP_GEOMETRY, aucdf)
            if nrow(sub) == 0
                text!(ax, 0.5, 0.5, text="(no data)", align=(:center,:center))
                continue
            end
            hm = make_heatmap_panel!(ax, sub; xcol=xaxis, ycol=yaxis, zcol=:dAUC_rel,
                title = "$(sc) — $(HEATMAP_GEOMETRY)",
                xlabel = "corr (consumer–prey niche correlation)",
                ylabel = (yaxis == :k_prey ? "k_prey (proxy for connectance)" : string(yaxis)),
                clim = (0.0, 1.0)
            )
            push!(hms2, hm)
        end
        if !isempty(hms2)
            Colorbar(fig_rel[:, end+1], hms2[1], label="relative ΔAUC = (A-only − AB)/A-only")
        end
        Label(fig_rel[0, :], "Relative ΔAUC — $(HEATMAP_GEOMETRY)", fontsize=18)

        display(fig_rel)
        if SAVE_FIGS
            save(joinpath(OUT_DIR, outbase * "_rel.png"), fig_rel, px_per_unit=2)
            save(joinpath(OUT_DIR, outbase * "_rel.pdf"), fig_rel)
        end
    end

    return fig_raw, fig_rel
end

# ----------------------------- FIGURE 3 (MECH) -------------------------------
"""
Mechanism panel: expects a CSV with columns:
  f, mean_phi, supported_lcc, frag_fail
If not found, it will skip gracefully.

Place it as: RESULTS_DIR/mechanism_summary.csv (or any CSV in RESULTS_DIR matching these columns).
"""
function make_fig3_mechanism(results_dir::String; outbase="Fig3_mechanism")
    mech_path = nothing
    # try exact first
    p = joinpath(results_dir, "mechanism_summary.csv")
    if isfile(p)
        mech_path = p
    else
        # find any csv with the needed columns
        mech_path = find_first_csv_with_columns(results_dir, [:f, :mean_phi])
        if mech_path === nothing
            # allow alternate column spellings
            for cand in sort(glob("*.csv", results_dir))
                dfh = try CSV.read(cand, DataFrame; limit=5) catch; continue end
                cols = Set(Symbol.(names(dfh)))
                if (:f in cols) && ( (:meanφ in cols) || (:mean_phi in cols) )
                    mech_path = cand
                    break
                end
            end
        end
    end

    mech_path === nothing && begin
        @warn "No mechanism CSV found (mechanism_summary.csv). Skipping Fig3."
        return nothing
    end

    df = try_read_csv(mech_path)
    df[!, :f] = Float64.(df[!, :f])

    # normalize names
    if :meanφ in names(df) && !(:mean_phi in names(df)); rename!(df, :meanφ => :mean_phi) end
    if :supportedLCC in names(df) && !(:supported_lcc in names(df)); rename!(df, :supportedLCC => :supported_lcc) end
    if :fragFail in names(df) && !(:frag_fail in names(df)); rename!(df, :fragFail => :frag_fail) end

    needed = [:mean_phi, :supported_lcc, :frag_fail]
    miss = filter(c -> !(c in names(df)), needed)
    !isempty(miss) && begin
        @warn "Mechanism CSV found ($mech_path) but missing columns $(miss). Skipping Fig3."
        return nothing
    end

    for c in needed
        df[!, c] = Float64.(df[!, c])
    end
    sort!(df, :f)

    fig = Figure(size=(1400, 350))
    geoms = ["random", "cluster", "front"]

    # If your mechanism file includes geometry, split; otherwise plot single curve repeated
    has_geom = :geometry in names(df)
    if has_geom
        df[!, :geometry] = string.(df[!, :geometry])
    end

    for (j, g) in enumerate(geoms)
        ax = Axis(fig[1, j], title="Mechanism: support amount vs support connectivity — $(g)",
                  xlabel="habitat loss f", ylabel="mean φ (consumers)")
        axr = Axis(fig[1, j], yaxisposition = :right, ylabel="supported-LCC / frag-fail")
        hidespines!(axr); hidexdecorations!(axr)

        sub = has_geom ? filter(:geometry => ==(g), df) : df
        if nrow(sub) == 0
            text!(ax, 0.5, 0.5, text="(no data)", align=(:center,:center))
            continue
        end
        lines!(ax, sub.f, sub.mean_phi, linewidth=3, label="mean φ")
        lines!(axr, sub.f, sub.supported_lcc, linestyle=:dash, linewidth=3)
        lines!(axr, sub.f, sub.frag_fail, linestyle=:dot, linewidth=3)
        if j == 1
            axislegend(ax; position=:lb, framevisible=true)
        end
        xlims!(ax, 0, 1)
        ylims!(ax, 0, 1)
        ylims!(axr, 0, 1)
    end

    display(fig)
    if SAVE_FIGS
        save(joinpath(OUT_DIR, outbase * ".png"), fig, px_per_unit=2)
        save(joinpath(OUT_DIR, outbase * ".pdf"), fig)
    end
    return fig
end

# ----------------------------- REPORTING -------------------------------------
function report_auc_summary(aucdf::DataFrame)
    @info "AUC table columns: $(names(aucdf))"
    if (:dAUC_raw in names(aucdf)) && (:dAUC_rel in names(aucdf))
        # quick ranked summary by scenario
        g = combine(groupby(aucdf, [:metric, :scenario_name, :geometry]),
            :dAUC_raw => mean => :mean_dAUC_raw,
            :dAUC_rel => mean => :mean_dAUC_rel
        )
        sort!(g, [:metric, :mean_dAUC_rel], rev=true)
        @info "Top mean relative ΔAUC (by metric/scenario/geometry):"
        show(first(g, min(12, nrow(g))); allcols=true, truncatelines=true); println()
    else
        # assume aucdf already is ΔAUC in auc_emin / auc_viability
        for col in (:auc_emin, :auc_viability)
            if col in names(aucdf)
                g = combine(groupby(aucdf, [:scenario_name, :geometry]),
                    col => mean => Symbol("mean_" * String(col))
                )
                sort!(g, Symbol("mean_" * String(col)), rev=true)
                @info "Top mean $(col) (interpreted as ΔAUC):"
                show(first(g, min(12, nrow(g))); allcols=true, truncatelines=true); println()
            end
        end
    end
end

# ----------------------------- MAIN ------------------------------------------
function main()
    # 1) Load AUC summary (needed at minimum for Fig2)
    auc_summary = load_auc_summary(RESULTS_DIR)

    # 2) Try load curves; if present, recompute AUC and use that for Fig2 (better, gives relative ΔAUC)
    curves = load_curve_data(RESULTS_DIR)

    if curves !== nothing
        # Ensure metric column exists (if absent, treat as Emin_patch only)
        if :metric ∉ names(curves)
            curves[!, :metric] .= "Emin_patch"
        end

        # 2a) Compute AUC table from curves (includes raw + relative ΔAUC)
        auc_from_curves = compute_auc_table_from_curves(curves)
        report_auc_summary(auc_from_curves)

        # 3) Fig1 (Emin main)
        if "Emin_patch" in unique(curves.metric)
            make_fig1(curves; metric="Emin_patch", outbase="Fig1_examples_Emin")
        else
            @warn "No Emin_patch curves detected; skipping Fig1 Emin."
        end

        # SI Fig1 (viability)
        if MAKE_SI_VIABILITY && ("viability" in unique(curves.metric))
            make_fig1(curves; metric="viability", outbase="FigS1_examples_viability")
        end

        # 4) Fig2 (use Emin and optionally viability)
        if any(auc_from_curves.metric .== "Emin_patch")
            dfE = filter(:metric => ==("Emin_patch"), auc_from_curves)
            make_fig2(dfE; metric=:dAUC_raw, outbase="Fig2_sweep_Emin")
        end
        if MAKE_SI_VIABILITY && any(auc_from_curves.metric .== "viability")
            dfV = filter(:metric => ==("viability"), auc_from_curves)
            make_fig2(dfV; metric=:dAUC_raw, outbase="FigS2_sweep_viability")
        end

    else
        @warn "No curve CSV found -> cannot build Fig1 or relative ΔAUC from curves."
        @warn "Proceeding with Fig2 using auc_summary.csv as ΔAUC (raw only)."

        report_auc_summary(auc_summary)

        # Fig2 (raw only) from existing ΔAUC columns
        if :auc_emin in names(auc_summary)
            make_fig2(auc_summary; metric=:auc_emin, outbase="Fig2_sweep_Emin_from_aucsummary")
        end
        if MAKE_SI_VIABILITY && (:auc_viability in names(auc_summary))
            make_fig2(auc_summary; metric=:auc_viability, outbase="FigS2_sweep_viability_from_aucsummary")
        end
    end

    # 5) Fig3 mechanism (optional CSV)
    make_fig3_mechanism(RESULTS_DIR)
end

main()
###############################################################################
# Notes (important for “definitive” quality):
# - Fig1 depends on having the richness curves CSV. If you don’t yet export it,
#   export one “long” CSV with richness vs f for each (scenario, geometry, corr, k_prey, metric, model).
# - Fig2 looks best if you run a denser parameter grid (corr and connectance/k)
#   so heatmaps show continuous transitions. The plotting code handles arbitrary grids.
# - If you add realised connectance to your outputs (column :connectance_real),
#   the heatmap y-axis will automatically switch from k_prey to that connectance.
###############################################################################
