#!/usr/bin/env julia
# ============================================================
# END-TO-END: CTmax pipeline + Order (phylo) confounding check
#
# What this script does (all in one run):
#   1) Rebuilds the *exact* CTmax analysis edge set:
#        EU + (GloBI TF/Olden with TF preferred), dedup,
#        then filter to edges where BOTH endpoints have CTmax.
#   2) Builds species->order from:
#        - ThermoFresh/ThermoTol (expects "order" column; case-insensitive)
#        - GlobTherm (expects "Order" column; case-insensitive)
#        - Olden: no order column (ignored)
#   3) Finds species missing order within the CTmax analysis set.
#   4) Auto-fills missing orders via GBIF (internet required).
#   5) Creates TWO datasets:
#        A) "all_edges": original CTmax edges
#        B) "filtered_no_same_order": drops edges where pred & prey share order
#           (only when both orders are known)
#   6) Re-runs your CTmax outputs for BOTH A and B:
#        - edge_table.csv
#        - predator_level.csv
#        - summary.csv
#        - edge scatter (colored by interaction sources)
#        - edge scatter (colored by CTmax source)
#        - node scatter (colored by CTmax source)
#        - null histogram
#   7) Writes a comparison summary CSV.
#
# Outputs:
#   outputs_ctmax_order_confounding/
#     - orders/
#         species_order_lookup_before_gbif.csv
#         missing_order_species_before_gbif.csv
#         gbif_order_patch.csv
#         gbif_unresolved.csv
#         species_order_lookup_after_gbif.csv
#     - all_edges/
#         (same files as original CTmax pipeline)
#     - filtered_no_same_order/
#         (same files as original CTmax pipeline)
#     - comparison_summary.csv
#
# Notes / choices:
#   - If one endpoint is missing order, we KEEP the edge in filtered set.
#     (Only drop when both orders are known AND identical.)
#   - CTmax priority: ThermoFresh > Olden > GlobTherm (same as your script)
#   - GloBI edge preference: ThermoFresh-supported > Olden-supported (same)
#   - GBIF lookup uses /species/match with kingdom=Animalia then fallback.
# ============================================================

using CSV, DataFrames, Statistics, Random
using CairoMakie
using HTTP, JSON3, URIs
using Dates

# ============================================================
# 0) PATHS (EDIT THESE)
# ============================================================
script_dir   = @__DIR__
project_root = script_dir

# Thermal data
thermtol_csv  = joinpath(project_root, "thermtol_comb_final.csv")
olden_csv     = joinpath(project_root, "Comte_Olden_Data_Imputed.csv")
globtherm_csv = joinpath(project_root, "GlobalTherm_upload_02_11_17.csv")

# Interaction data
eu_metaweb_csv = joinpath(project_root, "../data/TetraEU_pairwise_interactions.csv")
globi_tf_csv   = joinpath(project_root, "thermofresh_globi_metaweb_fish_predators.csv")
globi_old_csv  = joinpath(project_root, "../outputs_imputed_globi_edges.csv")

OUTROOT = joinpath(project_root, "outputs_ctmax_order_confounding")
isdir(OUTROOT) || mkpath(OUTROOT)

ORDDIR = joinpath(OUTROOT, "orders")
isdir(ORDDIR) || mkpath(ORDDIR)

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

function scatter_with_fit_color(x, y, group; title="", xlab="", ylab="", file="")
    good = isfinite.(x) .& isfinite.(y) .& .!ismissing.(group)
    sum(good) < 2 && return

    x = x[good]; y = y[good]
    g = String.(group[good])

    r = pearson_r(x, y)
    a, b = fit_line(x, y)

    levs = sort(unique(g))
    K = length(levs)

    basepal = Makie.ColorSchemes.Set1_9.colors
    palette = [basepal[mod1(i, length(basepal))] for i in 1:K]

    idx = Dict(l => i for (i, l) in enumerate(levs))
    point_colors = [palette[idx[v]] for v in g]

    fig = Figure(size=(980, 620))
    ax  = Axis(fig[1, 1], title=title, xlabel=xlab, ylabel=ylab)

    scatter!(ax, x, y; color=point_colors)
    xs = range(minimum(x), maximum(x), length=200)
    lines!(ax, xs, a .+ b .* xs)

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

function find_col_ci(df::DataFrame, candidates::Vector{Symbol})
    nms = names(df)
    nms_l = lowercase.(String.(nms))
    for c in candidates
        i = findfirst(==(lowercase(String(c))), nms_l)
        i === nothing || return nms[i]
    end
    return nothing
end

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

function norm_order(s::AbstractString)::String
    t = strip(String(s))
    t == "" && return ""
    t = replace(t, r"\s+" => " ")
    return t
end

# ============================================================
# 3) BUILD CTmax LOOKUP + CTmax SOURCE (priority TF > Olden > GlobTherm)
# ============================================================
tf_df = CSV.read(thermtol_csv, DataFrame; missingstring="", ntasks=1, pool=false)
metric_str = lowercase.(string.(coalesce.(tf_df.metric, "")))
tf_ct = tf_df[occursin.("ctmax", metric_str), :]

tf_ct = transform(tf_ct,
    :species => (col -> canon_taxon.(string.(coalesce.(col, "")))) => :canon,
    :tol     => (col -> as_float.(col)) => :ctmax
)
tf_ct = tf_ct[.!isempty.(tf_ct.canon) .& isfinite.(tf_ct.ctmax), :]
tf_ct_sum = combine(groupby(tf_ct, :canon), :ctmax => mean => :ctmax)
tf_lookup = Dict(r.canon => r.ctmax for r in eachrow(tf_ct_sum))

ol_df = CSV.read(olden_csv, DataFrame; missingstring="", ntasks=1, pool=false)
ol_tmp = transform(ol_df,
    :Species => (col -> canon_taxon.(string.(coalesce.(col, "")))) => :canon,
    Symbol("CTmax") => (col -> parse_comma_float.(col)) => :ctmax
)
ol_tmp = ol_tmp[.!isempty.(ol_tmp.canon) .& isfinite.(ol_tmp.ctmax), :]
ol_sum = combine(groupby(ol_tmp, :canon), :ctmax => mean => :ctmax)
ol_lookup = Dict(r.canon => r.ctmax for r in eachrow(ol_sum))

gt_df = CSV.read(globtherm_csv, DataFrame; missingstring="", ntasks=1, pool=false)
gt_genus_col   = find_col_ci(gt_df, [:Genus, :genus])
gt_species_col = find_col_ci(gt_df, [:Species, :species])

gt_genus   = gt_genus_col === nothing ? fill("", nrow(gt_df)) : safe_string_vec(gt_df[!, gt_genus_col])
gt_species = gt_species_col === nothing ? fill("", nrow(gt_df)) : safe_string_vec(gt_df[!, gt_species_col])

gt_canon_raw = strip.(gt_genus .* " " .* gt_species)
gt_canon = canon_taxon.(gt_canon_raw)

tmax_col  = find_col_ci(gt_df, [:Tmax, :tmax])
tmax2_col = find_col_ci(gt_df, [:Tmax_2, :tmax_2])

t1 = tmax_col === nothing ? fill(NaN, nrow(gt_df)) : as_float.(gt_df[!, tmax_col])
t2 = tmax2_col === nothing ? fill(NaN, nrow(gt_df)) : as_float.(gt_df[!, tmax2_col])

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

gt_tmp = DataFrame(canon=gt_canon, ctmax=gt_ctmax)
gt_tmp = gt_tmp[.!isempty.(gt_tmp.canon) .& isfinite.(gt_tmp.ctmax), :]
gt_sum = combine(groupby(gt_tmp, :canon), :ctmax => mean => :ctmax)
gt_lookup = Dict(r.canon => r.ctmax for r in eachrow(gt_sum))

ctmax_lookup = Dict{String,Float64}()
ctmax_source = Dict{String,String}()

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

@info "CTmax species: TF=$(length(tf_lookup)) Olden=$(length(ol_lookup)) GlobTherm=$(length(gt_lookup)) merged=$(length(ctmax_lookup))"

# ============================================================
# 4) LOAD & MERGE INTERACTIONS (EU + GloBI-best)
# ============================================================
eu_raw = CSV.read(eu_metaweb_csv, DataFrame; missingstring="", ntasks=1, pool=false)
eu = DataFrame(
    pred = canon_taxon.(string.(eu_raw.sourceTaxonName)),
    prey = canon_taxon.(string.(eu_raw.targetTaxonName))
)
filter!(r -> !isempty(r.pred) && !isempty(r.prey), eu)
eu.interaction_source .= "EU"
eu.evidence_source    .= "EU"

g_tf_raw = CSV.read(globi_tf_csv, DataFrame; missingstring="", ntasks=1, pool=false)
tf_pred_col = ("predator" in names(g_tf_raw)) ? :predator : (("pred" in names(g_tf_raw)) ? :pred : error("TF GloBI: no pred/predator column"))
tf_prey_col = ("prey" in names(g_tf_raw)) ? :prey : error("TF GloBI: no prey column")

g_tf = DataFrame(
    pred = canon_taxon.(string.(g_tf_raw[!, tf_pred_col])),
    prey = canon_taxon.(string.(g_tf_raw[!, tf_prey_col]))
)
filter!(r -> !isempty(r.pred) && !isempty(r.prey), g_tf)
g_tf.interaction_source .= "GloBI"
g_tf.evidence_source    .= "ThermoFresh"

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

globi_all = vcat(g_tf, g_ol)
evidence_rank(s::AbstractString) = (s == "ThermoFresh") ? 3 :
                                  (s == "Olden")       ? 2 : 1

globi_best = combine(groupby(globi_all, [:pred, :prey])) do sdf
    i = argmax(evidence_rank.(String.(sdf.evidence_source)))
    sdf[i:i, :]
end

all_edges = vcat(eu, globi_best)

merged_edges = combine(groupby(all_edges, [:pred, :prey])) do sdf
    inter = join(sort(unique(String.(sdf.interaction_source))), ";")
    evid  = join(sort(unique(String.(sdf.evidence_source))), ";")
    DataFrame(interaction_sources=[inter], evidence_sources=[evid])
end

@info "Merged unique interactions (EU + GloBI-best)" nrow(merged_edges)

edges_ct = merged_edges[
    haskey.(Ref(ctmax_lookup), merged_edges.pred) .&
    haskey.(Ref(ctmax_lookup), merged_edges.prey),
    :
]
@info "CTmax analysis edges (both endpoints have CTmax)" nrow(edges_ct)

species_ct = sort(unique(vcat(edges_ct.pred, edges_ct.prey)))

# ============================================================
# 5) SPECIES->ORDER LOOKUP FROM TF + GLOBTHERM (Olden has none)
#     Priority: TF > GlobTherm
# ============================================================
tf_order_col = find_col_ci(tf_df, [:order, :Order, :ORDER])
tf_ord_lookup = Dict{String,String}()
if tf_order_col !== nothing
    tf_canon = canon_taxon.(string.(coalesce.(tf_df[!, :species], "")))
    tf_ord   = norm_order.(string.(coalesce.(tf_df[!, tf_order_col], "")))
    for i in eachindex(tf_canon)
        sp = tf_canon[i]; ord = tf_ord[i]
        (!isempty(sp) && !isempty(ord)) || continue
        tf_ord_lookup[sp] = ord
    end
end
@info "ThermoFresh order species" length(tf_ord_lookup) "order_col" (tf_order_col === nothing ? "NONE" : String(tf_order_col))

gt_order_col = find_col_ci(gt_df, [:Order, :order, :ORDER])
gt_ord_lookup = Dict{String,String}()
if gt_order_col !== nothing
    gt_ord = norm_order.(string.(coalesce.(gt_df[!, gt_order_col], "")))
    for i in eachindex(gt_canon)
        sp = gt_canon[i]; ord = gt_ord[i]
        (!isempty(sp) && !isempty(ord)) || continue
        gt_ord_lookup[sp] = ord
    end
end
@info "GlobTherm order species" length(gt_ord_lookup) "order_col" (gt_order_col === nothing ? "NONE" : String(gt_order_col))

order_lookup = Dict{String,String}()
order_source = Dict{String,String}()

for (sp,o) in gt_ord_lookup
    order_lookup[sp] = o
    order_source[sp] = "GlobTherm"
end
for (sp,o) in tf_ord_lookup
    order_lookup[sp] = o
    order_source[sp] = "ThermoFresh"
end

lookup_before = DataFrame(
    species = species_ct,
    order = [get(order_lookup, sp, "") for sp in species_ct],
    order_source = [get(order_source, sp, "") for sp in species_ct]
)
CSV.write(joinpath(ORDDIR, "species_order_lookup_before_gbif.csv"), lookup_before)

missing_before = lookup_before[lookup_before.order .== "", :]
CSV.write(joinpath(ORDDIR, "missing_order_species_before_gbif.csv"), missing_before)
@info "Order coverage BEFORE GBIF: missing species" nrow(missing_before)

# ============================================================
# 6) GBIF AUTO-FILL FOR MISSING ORDERS
# ============================================================
GBIF_MATCH_URL = "https://api.gbif.org/v1/species/match"
RATE_LIMIT_SECONDS = 0.25
TIMEOUT_SECONDS = 20
MAX_RETRIES = 3
DEFAULT_KINGDOM = "Animalia"

function clean_species_name(s::AbstractString)::String
    t = strip(String(s))
    t = replace(t, r"\s+" => " ")
    return t
end

function gbif_match(name::String; kingdom::Union{Nothing,String}=DEFAULT_KINGDOM)
    q = Dict("name" => name)
    if kingdom !== nothing && kingdom != ""
        q["kingdom"] = kingdom
    end
    url = string(URI(GBIF_MATCH_URL; query=q))

    last_err = nothing
    for _ in 1:MAX_RETRIES
        try
            resp = HTTP.get(url; connect_timeout=TIMEOUT_SECONDS, read_timeout=TIMEOUT_SECONDS)
            if resp.status == 200
                return JSON3.read(String(resp.body))
            else
                last_err = "HTTP $(resp.status)"
            end
        catch e
            last_err = e
        end
        sleep(RATE_LIMIT_SECONDS)
    end
    return last_err
end

function extract_gbif_order(obj)
    if !(obj isa JSON3.Object)
        return ("", "ERROR", missing, missing, missing, string(obj))
    end
    mt = string(get(obj, :matchType, ""))
    conf = try Int(get(obj, :confidence, missing)) catch; missing end
    ord = string(get(obj, :order, ""))   # <--- GBIF field for order
    usageKey = get(obj, :usageKey, missing)
    familyKey = get(obj, :familyKey, missing) # kept for metadata, though not needed
    return (ord, mt, conf, usageKey, familyKey, "")
end

function gbif_order_for_species(sp::String)
    obj1 = gbif_match(sp; kingdom=DEFAULT_KINGDOM)
    ord1, mt1, conf1, uk1, fk1, note1 = extract_gbif_order(obj1)
    if ord1 != "" && mt1 != "NONE"
        return (ord1, "kingdom=Animalia", mt1, conf1, uk1, fk1, note1)
    end
    obj2 = gbif_match(sp; kingdom=nothing)
    ord2, mt2, conf2, uk2, fk2, note2 = extract_gbif_order(obj2)
    if ord2 != "" && mt2 != "NONE"
        return (ord2, "no_kingdom", mt2, conf2, uk2, fk2, note2)
    end
    return ("", "no_match", mt2, conf2, uk2, fk2, note2)
end

miss_species = unique(clean_species_name.(string.(coalesce.(missing_before.species, ""))))
miss_species = filter(!isempty, miss_species)

gbif_rows = DataFrame(
    species = String[],
    order = String[],
    gbif_mode = String[],
    matchType = String[],
    confidence = Union{Missing,Int}[],
    usageKey = Union{Missing,Int}[],
    familyKey = Union{Missing,Int}[],
    note = String[]
)

println("GBIF order lookup for $(length(miss_species)) missing-order species…")
println("Started: ", Dates.format(now(), dateformat"yyyy-mm-dd HH:MM:SS"))

for (i, sp) in enumerate(miss_species)
    ord, mode, mt, conf, uk, fk, note = gbif_order_for_species(sp)
    push!(gbif_rows, (sp, ord, mode, mt, conf, uk, fk, note))
    if i % 25 == 0 || i == length(miss_species)
        println("  $(i)/$(length(miss_species))  resolved=$(sum(gbif_rows.order .!= ""))  last='$(sp)'")
    end
    sleep(RATE_LIMIT_SECONDS)
end

CSV.write(joinpath(ORDDIR, "gbif_order_patch.csv"), gbif_rows)

gbif_resolved = gbif_rows[gbif_rows.order .!= "", :]
gbif_unresolved = gbif_rows[gbif_rows.order .== "", :]
CSV.write(joinpath(ORDDIR, "gbif_unresolved.csv"), gbif_unresolved)

# Patch order_lookup with GBIF where still missing
for r in eachrow(gbif_resolved)
    sp = r.species
    if get(order_lookup, sp, "") == ""
        order_lookup[sp] = r.order
        order_source[sp] = "GBIF"
    end
end

lookup_after = DataFrame(
    species = species_ct,
    order = [get(order_lookup, sp, "") for sp in species_ct],
    order_source = [get(order_source, sp, "") for sp in species_ct]
)
CSV.write(joinpath(ORDDIR, "species_order_lookup_after_gbif.csv"), lookup_after)

missing_after = lookup_after[lookup_after.order .== "", :]
CSV.write(joinpath(ORDDIR, "missing_order_species_after_gbif.csv"), missing_after)
@info "Order coverage AFTER GBIF: missing species" nrow(missing_after)

# ============================================================
# 7) BUILD EDGE TABLES (ALL vs FILTERED_NO_SAME_ORDER)
# ============================================================
pred_ord = [get(order_lookup, sp, "") for sp in edges_ct.pred]
prey_ord = [get(order_lookup, sp, "") for sp in edges_ct.prey]

both_known = (pred_ord .!= "") .& (prey_ord .!= "")
same_order = both_known .& (pred_ord .== prey_ord)

edges_ct_ord = DataFrame(
    pred = edges_ct.pred,
    prey = edges_ct.prey,
    interaction_sources = edges_ct.interaction_sources,
    evidence_sources = edges_ct.evidence_sources,
    pred_order = pred_ord,
    prey_order = prey_ord,
    both_orders_known = both_known,
    same_order = same_order
)

CSV.write(joinpath(ORDDIR, "ctmax_edges_with_order_after_gbif.csv"), edges_ct_ord)

edges_all = edges_ct_ord
edges_filtered = edges_ct_ord[.!edges_ct_ord.same_order, :]   # keep if missing order on either end OR different order

@info "Edges total (CTmax analysis)" nrow(edges_all)
@info "Edges removed (same-order, both known)" sum(edges_ct_ord.same_order)
@info "Edges remaining after filter" nrow(edges_filtered)

# ============================================================
# 8) RUN YOUR ORIGINAL CTmax PIPELINE OUTPUTS FOR A GIVEN EDGE SET
# ============================================================
function run_ctmax_outputs(edges_in::DataFrame, outdir::String; label::String)
    isdir(outdir) || mkpath(outdir)

    pred_val = [ctmax_lookup[p] for p in edges_in.pred]
    prey_val = [ctmax_lookup[p] for p in edges_in.prey]
    absdiff  = abs.(pred_val .- prey_val)

    ctsrc_pred = [ctmax_source[p] for p in edges_in.pred]
    ctsrc_prey = [ctmax_source[p] for p in edges_in.prey]

    edges = DataFrame(
        pred = edges_in.pred,
        prey = edges_in.prey,
        interaction_sources = edges_in.interaction_sources,
        evidence_sources = edges_in.evidence_sources,
        pred_ctmax = pred_val,
        prey_ctmax = prey_val,
        absdiff = absdiff,
        ctmax_source_pred = ctsrc_pred,
        ctmax_source_prey = ctsrc_prey
    )
    CSV.write(joinpath(outdir, "edge_table.csv"), edges)

    edges_clean = edges[isfinite.(edges.pred_ctmax) .& isfinite.(edges.prey_ctmax), :]

    pred_level = if isempty(edges_clean)
        DataFrame(pred=String[], pred_ctmax=Float64[], mean_prey_ctmax=Float64[], ctmax_source_pred=String[], nprey=Int[])
    else
        pred_grp = groupby(edges_clean, :pred)
        combine(pred_grp,
            :pred_ctmax => mean => :pred_ctmax,
            :prey_ctmax => mean => :mean_prey_ctmax,
            :ctmax_source_pred => (x -> first(x)) => :ctmax_source_pred,
            nrow => :nprey
        )
    end
    CSV.write(joinpath(outdir, "predator_level.csv"), pred_level)

    obs = isempty(edges_clean) ? NaN : mean(edges_clean.absdiff)

    null_vals = Float64[]
    clean_null = Float64[]
    μ0 = NaN; σ0 = NaN; z = NaN; p = NaN

    if nrow(edges_clean) >= 2
        base_preds = collect(edges_clean.pred)
        base_preys = collect(edges_clean.prey)
        null_vals = null_mean_absdiff(base_preds, base_preys, ctmax_lookup; nperm=1000, nsweeps=10, seed=1)
        clean_null = filter(isfinite, null_vals)
        μ0 = isempty(clean_null) ? NaN : mean(clean_null)
        σ0 = isempty(clean_null) ? NaN : std(clean_null)
        z  = (σ0 > 0 && isfinite(obs)) ? (obs - μ0) / σ0 : NaN
        p  = isempty(null_vals) || !isfinite(obs) ? NaN : mean(null_vals .<= obs)
    end

    if nrow(edges_clean) >= 2
        scatter_with_fit_color(
            edges_clean.pred_ctmax, edges_clean.prey_ctmax, edges_clean.interaction_sources;
            title="Edge-level: predator vs prey (CTmax) — $(label)",
            xlab="Predator CTmax",
            ylab="Prey CTmax",
            file=joinpath(outdir, "edge_pred_vs_prey_color_by_interactions.png")
        )

        scatter_with_fit_color(
            edges_clean.pred_ctmax, edges_clean.prey_ctmax, edges_clean.ctmax_source_pred;
            title="Edge-level: predator vs prey (CTmax) — $(label)",
            xlab="Predator CTmax",
            ylab="Prey CTmax",
            file=joinpath(outdir, "edge_pred_vs_prey_color_by_ctmax_source.png")
        )
    end

    if nrow(pred_level) >= 2
        scatter_with_fit_color(
            pred_level.pred_ctmax, pred_level.mean_prey_ctmax, pred_level.ctmax_source_pred;
            title="Node-level: predator vs mean(prey) (CTmax) — $(label)",
            xlab="Predator CTmax",
            ylab="Mean prey CTmax",
            file=joinpath(outdir, "node_pred_vs_meanprey_color_by_ctmax_source.png")
        )
    end

    begin
        fig = Figure(size=(900,570))
        ax = Axis(fig[1,1],
            title="Null: mean |ΔCTmax| (pred-swap) — $(label)",
            xlabel="Mean |ΔCTmax|",
            ylabel="Count"
        )
        isempty(clean_null) || hist!(ax, clean_null, bins=30)
        isfinite(obs) && vlines!(ax, [obs])
        save(joinpath(outdir, "null_mean_absdiff.png"), fig)
        display(fig)
    end

    summary = DataFrame(
        metric = "CTmax",
        dataset = label,
        edges = nrow(edges_clean),
        predators = length(unique(edges_clean.pred)),
        prey = length(unique(edges_clean.prey)),
        obs_mean_absdiff = obs,
        null_mean = μ0,
        null_sd = σ0,
        z = z,
        p_one_sided = p
    )
    CSV.write(joinpath(outdir, "summary.csv"), summary)
    return summary
end

# ============================================================
# 9) RUN BOTH + WRITE COMPARISON
# ============================================================
out_all = joinpath(OUTROOT, "all_edges")
out_fil = joinpath(OUTROOT, "filtered_no_same_order")

summ_all = run_ctmax_outputs(edges_all, out_all; label="ALL_EDGES")
summ_fil = run_ctmax_outputs(edges_filtered, out_fil; label="NO_SAME_ORDER")

comparison = vcat(summ_all, summ_fil)
CSV.write(joinpath(OUTROOT, "comparison_summary.csv"), comparison)

@info "DONE" OUTROOT
println("\n=== comparison_summary ===")
show(comparison, allrows=true, allcols=true)
println()
