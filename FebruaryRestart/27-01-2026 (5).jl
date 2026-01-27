using CSV, DataFrames, Statistics
using CairoMakie

# -------------------------
# CONFIG
# -------------------------
function default_auc_path()
    candidates = String[
        joinpath(pwd(), "Figures", "patchonly_36", "auc_summary.csv"),
        joinpath(pwd(), "auc_summary.csv"),
        joinpath(pwd(), "Figures", "patchonly", "auc_summary.csv"),
    ]
    for p in candidates
        if isfile(p)
            return p
        end
    end
    error("Could not find auc_summary.csv. Pass path as ARGS[1].")
end

AUC_PATH = length(ARGS) >= 1 ? ARGS[1] : default_auc_path()

ENFORCE_SANITY = true
TOL_NEG        = 1e-9
TOL_REL        = 1e-3
VERBOSE_REPORT = true

# -------------------------
# HELPERS
# -------------------------
function pickcol(df::DataFrame, candidates::Vector{Symbol})
    for c in candidates
        c in names(df) && return c
    end
    return nothing
end

to_str(x) = string(x)
to_int(x) = Int(round(Float64(x)))
to_flt(x) = Float64(x)

nice_key(niche, geom) = string(niche, " | ", geom)

# -------------------------
# LOAD
# -------------------------
df = CSV.read(AUC_PATH, DataFrame)

# Explicit column names (fixed schema)
col_niche = :scenario_name
col_geom  = :geometry
col_k     = :k_prey
col_corr  = :corr

if any(isnothing, (col_niche, col_geom, col_k, col_corr))
    error("Missing required ID columns. Found: $(names(df))")
end

col_aucE = :auc_emin
col_aucV = :auc_viability

if col_aucE === nothing || col_aucV === nothing
    error("Missing AUC columns. Found: $(names(df))")
end

col_S0A  = pickcol(df, [:S0_A, :S0_Aonly, :S0_A_only, :S0])
col_fmin = pickcol(df, [:fmin, :hl_min])
col_fmax = pickcol(df, [:fmax, :hl_max])

# -------------------------
# COERCE TYPES
# -------------------------
df[!, col_niche] = to_str.(df[!, col_niche])
df[!, col_geom]  = to_str.(df[!, col_geom])
df[!, col_k]     = to_int.(df[!, col_k])
df[!, col_corr]  = to_flt.(df[!, col_corr])
df[!, col_aucE]  = to_flt.(df[!, col_aucE])
df[!, col_aucV]  = to_flt.(df[!, col_aucV])

if col_S0A !== nothing
    df[!, col_S0A] = to_flt.(df[!, col_S0A])
end
if col_fmin !== nothing && col_fmax !== nothing
    Δf = median(df[!, col_fmax] .- df[!, col_fmin])
else
    Δf = 1.0
end

# -------------------------
# LONG FORMAT
# -------------------------
function build_long(df::DataFrame; metric_name::String, auc_col::Symbol)
    out = DataFrame(
        niche    = df[!, col_niche],
        geometry = df[!, col_geom],
        k_prey   = df[!, col_k],
        corr     = df[!, col_corr],
        metric   = fill(metric_name, nrow(df)),
        auc_raw  = df[!, auc_col],
    )

    if col_S0A !== nothing
        max_auc = df[!, col_S0A] .* Δf
        out.auc_rel = out.auc_raw ./ max_auc
    else
        out.auc_rel = fill(NaN, nrow(df))
    end
    return out
end

d = vcat(
    build_long(df; metric_name="Emin_patch", auc_col=col_aucE),
    build_long(df; metric_name="viability",  auc_col=col_aucV)
)

if col_S0A === nothing
    for m in unique(d.metric)
        idxx = findall(d.metric .== m)
        mx  = maximum(d.auc_raw[idxx])
        d.auc_rel[idxx] .= d.auc_raw[idxx] ./ (mx == 0 ? 1.0 : mx)
    end
end

# -------------------------
# SANITY CHECKS
# -------------------------
d.flag_neg = d.auc_raw .< -TOL_NEG
d.flag_rel = (d.auc_rel .< -TOL_REL) .| (d.auc_rel .> 1 + TOL_REL)

if VERBOSE_REPORT
    println("\n--- AUC POST-PASS REPORT ---")
    println("File: ", AUC_PATH)
    println("Rows total: ", nrow(d))
    println("Negative ΔAUC rows: ", count(d.flag_neg))
    println("Relative outside [0,1]: ", count(d.flag_rel))
end

if ENFORCE_SANITY
    d = d[.!d.flag_neg .& .!d.flag_rel, :]
end

d.key = nice_key.(d.niche, d.geometry)

# -------------------------
# PLOTS (DISPLAY ONLY)
# -------------------------
function plot_bar_means(d::DataFrame; value_col::Symbol, title::String, ylabel::String)
    g = combine(groupby(d, [:key, :metric]), value_col => mean => :mean_auc)
    keys = unique(g.key)
    metrics = ["Emin_patch", "viability"]

    meanE = Float64[]
    meanV = Float64[]
    for k in keys
        rE = g[(g.key .== k) .& (g.metric .== metrics[1]), :mean_auc]
        rV = g[(g.key .== k) .& (g.metric .== metrics[2]), :mean_auc]

        push!(meanE, isempty(rE) ? NaN : rE[1])
        push!(meanV, isempty(rV) ? NaN : rV[1])
    end


    fig = Figure(size=(1600, 550))
    ax = Axis(fig[1, 1], title=title, ylabel=ylabel)

    x = 1:length(keys)
    w = 0.35
    barplot!(ax, x .- w/2, meanE; width=w, label=metrics[1])
    barplot!(ax, x .+ w/2, meanV; width=w, label=metrics[2])

    ax.xticks = (x, keys)
    ax.xticklabelrotation = pi/4
    ax.xticklabelalign = (:right, :center)
    axislegend(ax, position=:rt)

    display(fig)
end

plot_bar_means(d; value_col=:auc_raw,
    title="Mean raw ΔAUC (A-only − AB)",
    ylabel="mean ΔAUC"
)

plot_bar_means(d; value_col=:auc_rel,
    title="Mean relative ΔAUC (A-only − AB)",
    ylabel="relative ΔAUC"
)

# -------------------------
# HEATMAP FACETS
# -------------------------
function plot_heatmap_facets(d::DataFrame; metric::String, value_col::Symbol,
                             title::String, cbar_label::String)

    dd = d[d.metric .== metric, :]
    niches = unique(dd.niche)
    geoms  = unique(dd.geometry)
    ks     = sort(unique(dd.k_prey))
    cors   = sort(unique(dd.corr))

    fig = Figure(size=(1600, 1200))
    Label(fig[0, 1:length(geoms)], title, fontsize=18)

    vals = skipmissing(dd[!, value_col])
    vmin, vmax = extrema(vals)
    vmin == vmax && (vmin -= 1; vmax += 1)

    hm_last = nothing
    for (i, niche) in enumerate(niches), (j, geom) in enumerate(geoms)
        ax = Axis(fig[i, j],
            title = "$niche — $geom",
            xlabel = "corr",
            ylabel = "k_prey"
        )

        sub = dd[(dd.niche .== niche) .& (dd.geometry .== geom), :]
        M = fill(NaN, length(ks), length(cors))

        for r in eachrow(sub)
            M[findfirst(==(r.k_prey), ks),
              findfirst(==(r.corr), cors)] = r[value_col]
        end

        hm_last = heatmap!(ax, cors, ks, M; colorrange=(vmin, vmax))
    end

    Colorbar(fig[:, end+1], hm_last, label=cbar_label)
    display(fig)
end

plot_heatmap_facets(d; metric="Emin_patch", value_col=:auc_raw,
    title="Raw ΔAUC — Emin_patch", cbar_label="ΔAUC")

plot_heatmap_facets(d; metric="viability", value_col=:auc_raw,
    title="Raw ΔAUC — viability", cbar_label="ΔAUC")

plot_heatmap_facets(d; metric="Emin_patch", value_col=:auc_rel,
    title="Relative ΔAUC — Emin_patch", cbar_label="relative ΔAUC")

plot_heatmap_facets(d; metric="viability", value_col=:auc_rel,
    title="Relative ΔAUC — viability", cbar_label="relative ΔAUC")

println("\nDone. (Plots displayed; nothing saved.)")
