#!/usr/bin/env julia
# ============================================================
# MERGED CTmax × (EU metaweb + GloBI) metaweb pipeline
#
# Thermal sources (CTmax):
#   1) ThermoFresh/ThermoTol (thermtol_comb_final.csv; metric == "CTmax")
#   2) Olden imputed (Comte_Olden_Data_Imputed.csv; column "CTmax")
#   3) GlobTherm (GlobalTherm_upload_02_11_17.csv; tries to find CTmax column)
#
# Interaction sources:
#   - EU metaweb edges (CSV with predator/prey columns)
#   - GloBI edges built from:
#       * ThermoFresh-filtered GloBI fish metaweb
#       * Olden-filtered GloBI edges
#   Dedup rule inside GloBI: prefer ThermoFresh-supported edges over Olden-supported.
#
# Outputs:
#   outputs_merged_all_ctmax/
#     - edge_table.csv
#     - predator_level.csv
#     - summary.csv
#     - edge_pred_vs_prey_color_by_interactions.png
#     - edge_pred_vs_prey_color_by_ctmax_source.png
#     - node_pred_vs_meanprey_color_by_ctmax_source.png
#     - null_mean_absdiff.png
# ============================================================
using CSV, DataFrames, Statistics, Random
using CairoMakie

# ============================================================
# 0) PATHS (EDIT THESE)
# ============================================================
script_dir   = @__DIR__
project_root = script_dir

# Thermal data
thermtol_csv = joinpath(project_root, "thermtol_comb_final.csv")  # ThermoFresh/ThermoTol long table
olden_csv    = joinpath(project_root, "Comte_Olden_Data_Imputed.csv")
globtherm_csv = joinpath(project_root, "GlobalTherm_upload_02_11_17.csv")

# Interaction data
eu_metaweb_csv = joinpath(project_root, "../data/TetraEU_pairwise_interactions.csv")
globi_tf_csv   = joinpath(project_root, "thermofresh_globi_metaweb_fish_predators.csv")
globi_old_csv  = joinpath(project_root, "../outputs_imputed_globi_edges.csv")

# Output
OUTROOT = joinpath(project_root, "outputs_merged_all_ctmax")
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

as_float(x) =
    x === missing ? NaN :
    x isa Real ? Float64(x) :
    x isa AbstractString ? (y=tryparse(Float64, strip(x)); y===nothing ? NaN : y) :
    NaN

pearson_r(x,y) = (length(x) ≥ 3 && std(x)>0 && std(y)>0) ? cor(x,y) : NaN

function fit_line(x,y)
    b = cov(x,y)/var(x)
    a = mean(y) - b*mean(x)
    return a,b
end

# Categorical-color scatter with fit + legend
function scatter_with_fit_color(x, y, group; title="", xlab="", ylab="", file="")
    good = isfinite.(x) .& isfinite.(y) .& .!ismissing.(group)
    sum(good) < 2 && return

    x = x[good]; y = y[good]
    g = String.(group[good])

    r = pearson_r(x, y)
    a, b = fit_line(x, y)

    levs = sort(unique(g))
    K = length(levs)

    # A discrete categorical palette (Set1_9 has 9 colors; recycle if >9)
    basepal = Makie.ColorSchemes.Set1_9.colors
    palette = [basepal[mod1(i, length(basepal))] for i in 1:K]

    idx = Dict(l => i for (i, l) in enumerate(levs))
    point_colors = [palette[idx[v]] for v in g]

    fig = Figure(size=(980, 620))
    ax  = Axis(fig[1, 1], title=title, xlabel=xlab, ylabel=ylab)

    scatter!(ax, x, y; color=point_colors)  # <-- explicit colors (no colormap)

    xs = range(minimum(x), maximum(x), length=200)
    lines!(ax, xs, a .+ b .* xs)

    # Legend using the same palette
    handles = [MarkerElement(color=palette[i], marker=:circle) for i in 1:K]
    Legend(fig[1, 2], handles, levs, "Group"; tellwidth=false)

    Label(fig[2, 1:2], "r=$(round(r, digits=2)), n=$(length(x))", fontsize=13)

    file != "" && save(file, fig)
    display(fig)
end

function edge_swap(preds, rng; nsweeps=10)
    preds = copy(preds)
    E = length(preds)
    for _ in 1:(nsweeps * E)
        i,j = rand(rng, 1:E, 2)
        i==j && continue
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

# Find a CTmax-ish column in GlobTherm (best-effort)
function find_ctmax_col(df::DataFrame)
    nms = names(df)
    # prioritize exact-ish matches first
    candidates = [
        r"^ctmax$"i,
        r"ctmax"i,
        r"critical.*max"i,
        r"tmax"i
    ]
    for rx in candidates
        for nm in nms
            occursin(rx, String(nm)) && return nm
        end
    end
    return nothing
end

# ============================================================
# 3) BUILD CTmax LOOKUP WITH PRIORITY
#     ThermoFresh > Olden > GlobTherm
# ============================================================
# --- 3a) ThermoFresh/ThermoTol (long format) ---
tf_df = CSV.read(thermtol_csv, DataFrame; missingstring="", ntasks=1, pool=false)

# select CTmax rows robustly (case-insensitive substring match)
metric_str = lowercase.(string.(coalesce.(tf_df.metric, "")))
tf_ct = tf_df[occursin.("ctmax", metric_str), :]

tf_ct = transform(tf_ct,
    :species => (col -> canon_taxon.(string.(coalesce.(col, "")))) => :canon,
    :tol     => (col -> as_float.(col)) => :ctmax
)

tf_ct = tf_ct[.!isempty.(tf_ct.canon) .& isfinite.(tf_ct.ctmax), :]

# If multiple CTmax measures per species, average them
tf_ct_sum = combine(groupby(tf_ct, :canon), :ctmax => mean => :ctmax)
tf_lookup = Dict(r.canon => r.ctmax for r in eachrow(tf_ct_sum))

@info "ThermoFresh/ThermoTol CTmax species" length(tf_lookup)

# --- 3b) Olden imputed ---
ol_df = CSV.read(olden_csv, DataFrame; missingstring="", ntasks=1, pool=false)
ol_df = transform(ol_df,
    :Species => (col -> canon_taxon.(string.(coalesce.(col, "")))) => :canon,
    Symbol("CTmax") => (col -> parse_comma_float.(col)) => :ctmax
)
ol_df = ol_df[.!isempty.(ol_df.canon) .& isfinite.(ol_df.ctmax), :]
ol_sum = combine(groupby(ol_df, :canon), :ctmax => mean => :ctmax)
ol_lookup = Dict(r.canon => r.ctmax for r in eachrow(ol_sum))
@info "Olden CTmax species" length(ol_lookup)

# --- 3c) GlobTherm ---
gt_df = CSV.read(globtherm_csv, DataFrame; missingstring="", ntasks=1, pool=false)

# Defensive string extraction: handles missing AND #undef
function safe_string_vec(col)
    out = Vector{String}(undef, length(col))
    for i in eachindex(col)
        v = ""
        try
            x = col[i]
            v = (x === missing) ? "" : string(x)
        catch e
            if !(e isa UndefRefError)
                rethrow()
            end
            v = ""
        end
        out[i] = v
    end
    return out
end

gt_genus   = safe_string_vec(gt_df[!, :Genus])
gt_species = safe_string_vec(gt_df[!, :Species])

gt_canon_raw = strip.(gt_genus .* " " .* gt_species)
gt_canon = canon_taxon.(gt_canon_raw)

t1 = as_float.(gt_df[!, :Tmax])
t2 = (:Tmax_2 in names(gt_df)) ? as_float.(gt_df[!, :Tmax_2]) : fill(NaN, nrow(gt_df))

# mean of available finite values
gt_ctmax = similar(t1)
for i in eachindex(t1)
    a, b = t1[i], t2[i]
    if isfinite(a) && isfinite(b)
        gt_ctmax[i] = (a + b) / 2
    elseif isfinite(a)
        gt_ctmax[i] = a
    elseif isfinite(b)
        gt_ctmax[i] = b
    else
        gt_ctmax[i] = NaN
    end
end

tmp = DataFrame(canon = gt_canon, ctmax = gt_ctmax)
tmp = tmp[.!isempty.(tmp.canon) .& isfinite.(tmp.ctmax), :]

tmp_sum = combine(groupby(tmp, :canon), :ctmax => mean => :ctmax)
gt_lookup = Dict(r.canon => r.ctmax for r in eachrow(tmp_sum))

@info "GlobTherm CTmax species" length(gt_lookup)

# --- 3d) Merge lookups with priority ---
# Final lookup and CTmax source label
ctmax_lookup = Dict{String,Float64}()
ctmax_source = Dict{String,String}()

# lowest priority first, then overwrite
for (sp,v) in gt_lookup
    ctmax_lookup[sp] = v
    ctmax_source[sp] = "GlobTherm"
end
for (sp,v) in ol_lookup
    ctmax_lookup[sp] = v
    ctmax_source[sp] = "Olden"
end
for (sp,v) in tf_lookup
    ctmax_lookup[sp] = v
    ctmax_source[sp] = "ThermoFresh"
end

@info "Final merged CTmax species" length(ctmax_lookup)

# ============================================================
# 4) LOAD & MERGE INTERACTIONS (EU + GloBI)
# ============================================================
# --- 4a) EU metaweb ---
eu_raw = CSV.read(eu_metaweb_csv, DataFrame; missingstring="", ntasks=1, pool=false)

eu = DataFrame(
    pred = canon_taxon.(string.(eu_raw.sourceTaxonName)),
    prey = canon_taxon.(string.(eu_raw.targetTaxonName))
)

filter!(r -> !isempty(r.pred) && !isempty(r.prey), eu)

# Tag EU edges
eu.interaction_source .= "EU"
eu.evidence_source    .= "EU"

# --- 4b) GloBI edges from ThermoFresh-filtered build ---
g_tf_raw = CSV.read(globi_tf_csv, DataFrame; missingstring="", ntasks=1, pool=false)

# Adjust column names if needed (your TF GloBI file may be predator/prey or pred/prey)
tf_pred_col = ("predator" in names(g_tf_raw)) ? :predator : (("pred" in names(g_tf_raw)) ? :pred : error("TF GloBI: no pred/predator column"))
tf_prey_col = ("prey" in names(g_tf_raw)) ? :prey : error("TF GloBI: no prey column")

g_tf = DataFrame(
    pred = canon_taxon.(string.(g_tf_raw[!, tf_pred_col])),
    prey = canon_taxon.(string.(g_tf_raw[!, tf_prey_col]))
)

filter!(r -> !isempty(r.pred) && !isempty(r.prey), g_tf)

g_tf.interaction_source .= "GloBI"
g_tf.evidence_source    .= "ThermoFresh"  # this is what wins ties inside GloBI

# --- 4c) GloBI edges from Olden-filtered build ---
g_ol_raw = CSV.read(globi_old_csv, DataFrame; missingstring="", ntasks=1, pool=false)

ol_pred_col = ("pred" in names(g_ol_raw)) ? :pred : (("predator" in names(g_ol_raw)) ? :predator : error("Olden GloBI: no pred/predator column"))
ol_prey_col = ("prey" in names(g_ol_raw)) ? :prey : error("Olden GloBI: no prey column")

g_ol = DataFrame(
    pred = canon_taxon.(string.(g_ol_raw[!, ol_pred_col])),
    prey = canon_taxon.(string.(g_ol_raw[!, ol_prey_col]))
)

filter!(r -> !isempty(r.pred) && !isempty(r.prey), g_ol)

g_ol.interaction_source .= "GloBI"
g_ol.evidence_source    .= "Olden"

# --- 4d) Apply preference within GloBI: ThermoFresh > Olden ---
globi_all = vcat(g_tf, g_ol)

evidence_rank(s::AbstractString) = (s == "ThermoFresh") ? 3 :
                                  (s == "Olden")       ? 2 : 1

globi_best = combine(groupby(globi_all, [:pred, :prey])) do sdf
    i = argmax(evidence_rank.(String.(sdf.evidence_source)))
    sdf[i:i, :]
end

@info "GloBI edges (all, before preference)" nrow(globi_all)
@info "GloBI edges (after preference)" nrow(globi_best)

# --- 4e) Merge EU + best GloBI, then collapse duplicates across DBs ---
all_edges = vcat(eu, globi_best)

merged_edges = combine(groupby(all_edges, [:pred, :prey])) do sdf
    inter = join(sort(unique(String.(sdf.interaction_source))), ";")
    evid  = join(sort(unique(String.(sdf.evidence_source))), ";")
    DataFrame(
        interaction_sources = [inter],
        evidence_sources    = [evid]
    )
end

@info "EU edges" nrow(eu)
@info "Merged unique interactions (EU + GloBI-best)" nrow(merged_edges)

# ============================================================
# 5) FILTER TO EDGES WITH CTmax ON BOTH ENDS + BUILD EDGE TABLE
# ============================================================
edges0 = merged_edges[
    haskey.(Ref(ctmax_lookup), merged_edges.pred) .&
    haskey.(Ref(ctmax_lookup), merged_edges.prey),
    :
]

@info "Edges with CTmax on both ends (after merge)" nrow(edges0)

pred_val = [ctmax_lookup[p] for p in edges0.pred]
prey_val = [ctmax_lookup[p] for p in edges0.prey]
absdiff  = abs.(pred_val .- prey_val)

ctsrc_pred = [ctmax_source[p] for p in edges0.pred]
ctsrc_prey = [ctmax_source[p] for p in edges0.prey]

edges = DataFrame(
    pred = edges0.pred,
    prey = edges0.prey,
    interaction_sources = edges0.interaction_sources,
    evidence_sources = edges0.evidence_sources,
    pred_ctmax = pred_val,
    prey_ctmax = prey_val,
    absdiff = absdiff,
    ctmax_source_pred = ctsrc_pred,
    ctmax_source_prey = ctsrc_prey
)

CSV.write(joinpath(OUTROOT, "edge_table.csv"), edges)

# ============================================================
# 6) NODE-LEVEL TABLE (predator vs mean prey)
# ============================================================
edges_clean = edges[isfinite.(edges.pred_ctmax) .& isfinite.(edges.prey_ctmax), :]

pred_grp = groupby(edges_clean, :pred)
pred_level = combine(pred_grp,
    :pred_ctmax => mean => :pred_ctmax,
    :prey_ctmax => mean => :mean_prey_ctmax,
    :ctmax_source_pred => (x -> first(x)) => :ctmax_source_pred,
    nrow => :nprey
)

CSV.write(joinpath(OUTROOT, "predator_level.csv"), pred_level)

# ============================================================
# 7) NULL MODEL (swap predator identities)
# ============================================================
base_preds = collect(edges_clean.pred)
base_preys = collect(edges_clean.prey)

null_vals = null_mean_absdiff(base_preds, base_preys, ctmax_lookup; nperm=1000, nsweeps=10, seed=1)

obs = mean(edges_clean.absdiff)
clean_null = filter(isfinite, null_vals)
μ0  = isempty(clean_null) ? NaN : mean(clean_null)
σ0  = isempty(clean_null) ? NaN : std(clean_null)
z   = (σ0 > 0) ? (obs - μ0) / σ0 : NaN
p   = mean(null_vals .<= obs)  # one-sided: unusually small mean |Δ|?

# ============================================================
# 8) PLOTS (colored)
# ============================================================
# Edge-level colored by interaction DB(s)
scatter_with_fit_color(
    edges_clean.pred_ctmax, edges_clean.prey_ctmax, edges_clean.interaction_sources;
    title="Edge-level: predator vs prey (CTmax)",
    xlab="Predator CTmax",
    ylab="Prey CTmax",
    file=joinpath(OUTROOT, "edge_pred_vs_prey_color_by_interactions.png")
)

# Edge-level colored by CTmax source of predator (ThermoFresh/Olden/GlobTherm)
scatter_with_fit_color(
    edges_clean.pred_ctmax, edges_clean.prey_ctmax, edges_clean.ctmax_source_pred;
    title="Edge-level: predator vs prey (CTmax)",
    xlab="Predator CTmax",
    ylab="Prey CTmax",
    file=joinpath(OUTROOT, "edge_pred_vs_prey_color_by_ctmax_source.png")
)

# Node-level colored by CTmax source of predator
scatter_with_fit_color(
    pred_level.pred_ctmax, pred_level.mean_prey_ctmax, pred_level.ctmax_source_pred;
    title="Node-level: predator vs mean(prey) (CTmax)",
    xlab="Predator CTmax",
    ylab="Mean prey CTmax",
    file=joinpath(OUTROOT, "node_pred_vs_meanprey_color_by_ctmax_source.png")
)

# Null histogram
begin
    fig = Figure(size=(900,570))
    ax = Axis(fig[1,1],
        title="Null: mean |ΔCTmax| (pred-swap)",
        xlabel="Mean |ΔCTmax|",
        ylabel="Count"
    )
    isempty(clean_null) || hist!(ax, clean_null, bins=30)
    vlines!(ax, [obs])
    save(joinpath(OUTROOT, "null_mean_absdiff.png"), fig)
    display(fig)
end

# ============================================================
# 9) SUMMARY
# ============================================================
summary = DataFrame(
    metric = "CTmax",
    edges = nrow(edges_clean),
    predators = length(unique(edges_clean.pred)),
    prey = length(unique(edges_clean.prey)),
    obs_mean_absdiff = obs,
    null_mean = μ0,
    null_sd = σ0,
    z = z,
    p_one_sided = p
)
CSV.write(joinpath(OUTROOT, "summary.csv"), summary)

@info "DONE" summary
