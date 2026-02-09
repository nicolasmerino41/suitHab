#!/usr/bin/env julia
# ============================================================
# Olden (CTmax) × GloBI metaweb
# FULL METRIC PIPELINE (single metric: CTmax)
#
# Inputs:
#   - FebruaryRestart/Comte_Olden_Data_Imputed.csv
#   - outputs_imputed_globi_edges.csv   (pred, prey)
#
# Outputs (folder):
#   outputs_empirical_olden_globi_ctmax/CTmax/
#     - edge_table_CTmax.csv
#     - predator_level_CTmax.csv
#     - summary_CTmax.csv
#     - edge_pred_vs_prey_CTmax.png
#     - node_pred_vs_meanprey_CTmax.png
#     - null_mean_absdiff_CTmax.png
# ============================================================

using CSV, DataFrames, Statistics, Random
using CairoMakie

# ============================================================
# 0) PATHS
# ============================================================
script_dir   = @__DIR__
project_root = script_dir

imputed_csv = joinpath(project_root, "FebruaryRestart", "Comte_Olden_Data_Imputed.csv")
edges_csv   = joinpath(project_root, "outputs_imputed_globi_edges.csv")

OUTROOT = joinpath(project_root, "outputs_empirical_olden_globi_ctmax", "CTmax")
isdir(OUTROOT) || mkpath(OUTROOT)

# ============================================================
# 1) CANONICAL TAXON
# ============================================================
function canon_taxon(s::AbstractString)::String
    t = replace(String(s), r"[_(),\[\]]" => " ")
    t = replace(t, r"\s+" => " ")
    t = strip(t)
    parts = split(t)
    length(parts) < 2 && return ""

    genus   = String(filter(isletter, parts[1]))
    species = String(filter(isletter, parts[2]))

    # FIX: correct boolean logic (avoid precedence bugs)
    (isempty(genus) || isempty(species)) && return ""
    lowercase(species) in ("sp","spp") && return ""

    return uppercase(genus[1]) * lowercase(genus[2:end]) * " " * lowercase(species)
end

# ============================================================
# 2) HELPERS
# ============================================================
parse_comma_float(x) =
    x === missing ? NaN :
    x isa Real ? Float64(x) :
    x isa AbstractString ? (y = tryparse(Float64, replace(strip(x), "," => ".")); y === nothing ? NaN : y) :
    NaN

pearson_r(x,y) = (length(x) ≥ 3 && std(x) > 0 && std(y) > 0) ? cor(x,y) : NaN

function fit_line(x,y)
    b = cov(x,y) / var(x)
    a = mean(y) - b * mean(x)
    return a,b
end

function scatter_with_fit(x,y,title,xlab,ylab,file)
    good = isfinite.(x) .& isfinite.(y)
    sum(good) < 2 && return

    x = x[good]; y = y[good]

    r = pearson_r(x,y)
    a,b = fit_line(x,y)

    fig = Figure(size=(900,570))
    ax = Axis(fig[1,1], title=title, xlabel=xlab, ylabel=ylab)
    scatter!(ax, x, y)

    xs = range(minimum(x), maximum(x), length=200)
    lines!(ax, xs, a .+ b .* xs)

    Label(fig[2,1:2], "r=$(round(r,digits=2)), n=$(length(x))", fontsize=13)
    save(file, fig)
    display(fig)
end

function edge_swap(preds, rng; nsweeps=10)
    preds = copy(preds)
    E = length(preds)
    for _ in 1:(nsweeps * E)
        i,j = rand(rng, 1:E, 2)
        i == j && continue
        preds[i], preds[j] = preds[j], preds[i]
    end
    return preds
end

function null_mean_absdiff(base_preds, base_preys, lookup; nperm=1000, nsweeps=10, seed=1)
    rng = MersenneTwister(seed)
    E = length(base_preds)

    vals = Float64[]
    sizehint!(vals, nperm)

    for _ in 1:nperm
        p2 = edge_swap(base_preds, rng; nsweeps=nsweeps)
        s = 0.0
        for i in 1:E
            s += abs(lookup[p2[i]] - lookup[base_preys[i]])
        end
        push!(vals, s / E)
    end
    return vals
end

# ============================================================
# 3) LOAD OLDEN (CTmax) + BUILD LOOKUP
# ============================================================
imputed_df = CSV.read(imputed_csv, DataFrame; missingstring="", ntasks=1, pool=false)

# Canonize Species and parse CTmax
imputed_df = transform(imputed_df,
    :Species => (col -> canon_taxon.(string.(coalesce.(col, "")))) => :canon,
    Symbol("CTmax") => (col -> parse_comma_float.(col)) => :CTmax
)

# Filter valid
imputed_df = imputed_df[
    (.!isempty.(imputed_df.canon)) .&
    isfinite.(imputed_df.CTmax),
    :
]

# If multiple CTmax per canon, average them
ct_sum = combine(groupby(imputed_df, :canon), :CTmax => mean => :CTmax)
lookup = Dict(r.canon => r.CTmax for r in eachrow(ct_sum))

@info "Olden CTmax species (canon, mean aggregated)" length(lookup)

# ============================================================
# 4) LOAD EDGE LIST AND MAP METRIC
# ============================================================
mw = CSV.read(edges_csv, DataFrame; missingstring="", ntasks=1, pool=false)

# Be robust to column names (your file likely has :pred, :prey)
pred_col = ("pred" in names(mw)) ? :pred : (("predator" in names(mw)) ? :predator : error("No pred/predator column found"))
prey_col = ("prey" in names(mw)) ? :prey : error("No prey column found")

mw2 = DataFrame(
    pred = canon_taxon.(string.(mw[!, pred_col])),
    prey = canon_taxon.(string.(mw[!, prey_col]))
)

edges0 = mw2[
    haskey.(Ref(lookup), mw2.pred) .&
    haskey.(Ref(lookup), mw2.prey),
    :
]

@info "Edges with CTmax on both ends" nrow(edges0)

E = nrow(edges0)
E < 5 && @warn "Very few edges after CTmax mapping" E

pred_val = [lookup[p] for p in edges0.pred]
prey_val = [lookup[p] for p in edges0.prey]
absdiff  = abs.(pred_val .- prey_val)

edges = DataFrame(
    pred = edges0.pred,
    prey = edges0.prey,
    pred_val = pred_val,
    prey_val = prey_val,
    absdiff = absdiff
)

CSV.write(joinpath(OUTROOT, "edge_table_CTmax.csv"), edges)

# ============================================================
# 5) NODE-LEVEL TABLE (predator vs mean prey)
# ============================================================
edges_clean = edges[isfinite.(edges.pred_val) .& isfinite.(edges.prey_val), :]

pred_grp = groupby(edges_clean, :pred)
pred_level = combine(pred_grp,
    :pred_val => mean => :pred_val,
    :prey_val => mean => :mean_prey_val,
    nrow => :nprey
)

CSV.write(joinpath(OUTROOT, "predator_level_CTmax.csv"), pred_level)

# ============================================================
# 6) NULL MODEL (swap predator identities)
# ============================================================
base_preds = collect(edges_clean.pred)
base_preys = collect(edges_clean.prey)

null_vals = null_mean_absdiff(base_preds, base_preys, lookup; nperm=1000, nsweeps=10, seed=1)

obs = mean(edges_clean.absdiff)
μ0  = mean(filter(isfinite, null_vals))
σ0  = std(filter(isfinite, null_vals))
z   = (obs - μ0) / (σ0 > 0 ? σ0 : NaN)
p   = mean(null_vals .<= obs)  # one-sided: unusually small mean |Δ|?

# ============================================================
# 7) PLOTS
# ============================================================
scatter_with_fit(
    edges_clean.pred_val, edges_clean.prey_val,
    "Edge-level: predator vs prey (CTmax)",
    "Predator CTmax",
    "Prey CTmax",
    joinpath(OUTROOT, "edge_pred_vs_prey_CTmax.png")
)

scatter_with_fit(
    pred_level.pred_val, pred_level.mean_prey_val,
    "Node-level: predator vs mean(prey) (CTmax)",
    "Predator CTmax",
    "Mean prey CTmax",
    joinpath(OUTROOT, "node_pred_vs_meanprey_CTmax.png")
)

fig = Figure(size=(900,570))
ax = Axis(fig[1,1],
    title="Null: mean |ΔCTmax| (pred-swap)",
    xlabel="Mean |ΔCTmax|",
    ylabel="Count"
)

clean = filter(isfinite, null_vals)
isempty(clean) || hist!(ax, clean, bins=30)
vlines!(ax, [obs])
save(joinpath(OUTROOT, "null_mean_absdiff_CTmax.png"), fig)
display(fig)

# ============================================================
# 8) SUMMARY
# ============================================================
summary = DataFrame(
    metric = "CTmax",
    obs = obs,
    null_mean = μ0,
    null_sd = σ0,
    z = z,
    p_one_sided = p,
    edges = nrow(edges_clean),
    predators = length(unique(edges_clean.pred)),
    prey = length(unique(edges_clean.prey))
)

CSV.write(joinpath(OUTROOT, "summary_CTmax.csv"), summary)

@info "DONE" summary
