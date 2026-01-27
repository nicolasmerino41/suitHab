using CSV, DataFrames, Statistics
using CairoMakie

# -------------------------
# 1) Load CSV
# -------------------------
csvpath = joinpath(pwd(), "Figures", "patchonly_36", "auc_summary.csv")
df = CSV.read(csvpath, DataFrame)

# Explicit column names (fixed schema)
col_niche = :scenario_name
col_geom  = :geometry
col_k     = :k_prey
col_corr  = :corr

col_aucE = :auc_emin
col_aucV = :auc_viability

# Coerce to clean types
df[!, col_niche] = string.(df[!, col_niche])
df[!, col_geom]  = string.(df[!, col_geom])
df[!, col_k]     = Int.(round.(Float64.(df[!, col_k])))
df[!, col_corr]  = Float64.(df[!, col_corr])
df[!, col_aucE]  = Float64.(df[!, col_aucE])
df[!, col_aucV]  = Float64.(df[!, col_aucV])

# Unique ordered axes
niches = sort(unique(df[!, col_niche]))
geoms  = sort(unique(df[!, col_geom]))
ks     = sort(unique(df[!, col_k]))
cors   = sort(unique(df[!, col_corr]))

# -------------------------
# 2) Helper: mean heatmap matrix for a given niche × geom
# -------------------------
function mean_matrix(sub::DataFrame, col_metric::Symbol, ks::Vector{Int}, cors::Vector{Float64})
    M = fill(NaN, length(ks), length(cors))  # rows=k, cols=corr
    for (i, k) in pairs(ks), (j, c) in pairs(cors)
        s = sub[(sub[!, col_k] .== k) .& (sub[!, col_corr] .== c), :]
        if nrow(s) > 0
            M[i, j] = mean(skipmissing(s[!, col_metric]))
        end
    end
    return M
end

# -------------------------
# 3) Plot grid of heatmaps: niches (rows) × geometries (cols)
# -------------------------
function heatmap_grid(df::DataFrame; col_metric::Symbol, title::String)
    fig = Figure(size = (360 * length(geoms), 260 * length(niches)))

    # global color limits for comparability
    vals = df[!, col_metric]
    vmin, vmax = quantile(vals, 0.02), quantile(vals, 0.98)
    if isapprox(vmin, vmax; atol=1e-12)
        vmin -= 1e-6
        vmax += 1e-6
    end

    for (ri, niche) in pairs(niches)
        for (ci, geom) in pairs(geoms)
            ax = Axis(fig[ri, ci],
                title = "$(niche) — $(geom)",
                xlabel = "corr (consumer–prey niche correlation)",
                ylabel = "k_prey",
                titlesize = 9
            )

            sub = df[(df[!, col_niche] .== niche) .& (df[!, col_geom] .== geom), :]
            M = mean_matrix(sub, col_metric, ks, cors)

            heatmap!(ax, cors, ks, M; colorrange = (vmin, vmax))
            ax.xticks = (cors, string.(cors))
            ax.yticks = (ks, string.(ks))
        end
    end

    Colorbar(fig[:, end+1], limits = (vmin, vmax), label = title)
    Label(fig[0, :], title, fontsize = 18)

    display(fig)
end

# -------------------------
# 4) Heatmaps (shown immediately)
# -------------------------
figE = heatmap_grid(df;
    col_metric = col_aucE,
    title = "AUC(A-only − AB) under Emin_patch"
)

figV = heatmap_grid(df;
    col_metric = col_aucV,
    title = "AUC(A-only − AB) under viability"
)

# Difference grid: viability − Emin
df_diff = deepcopy(df)
df_diff[!, :auc_diff] = df_diff[!, col_aucV] .- df_diff[!, col_aucE]

figD = heatmap_grid(df_diff;
    col_metric = :auc_diff,
    title = "ΔAUC = viability − Emin (A-only − AB)"
)

# -------------------------
# 5) Bar summary: average over k_prey × corr
# -------------------------
g = combine(groupby(df, [col_niche, col_geom]),
    col_aucE => mean => :mean_auc_emin,
    col_aucV => mean => :mean_auc_viab
)

begin
    

    figB = Figure(size = (1100, 700))
    ax1 = Axis(figB[1, 1],
        title = "Mean AUC(A-only − AB) by niche × geometry",
        ylabel = "mean AUC"
    )

    cats = ["$(g[i, col_niche]) | $(g[i, col_geom])" for i in 1:nrow(g)]
    x = 1:length(cats)

    barplot!(ax1, x .- 0.2, g.mean_auc_emin; width = 0.35, label = "Emin_patch")
    barplot!(ax1, x .+ 0.2, g.mean_auc_viab; width = 0.35, label = "viability")

    ax1.xticks = (x, cats)
    ax1.xticklabelrotation = pi/4
    ax1.xticklabelalign = (:right, :center)

    axislegend(ax1; position = :rt)

    display(figB)
end