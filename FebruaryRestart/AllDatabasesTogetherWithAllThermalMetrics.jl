#!/usr/bin/env julia
# ============================================================
# MERGED THERMAL METRICS × (EU metaweb + GloBI) metaweb pipeline
#
# Runs the CTmax-style pipeline for multiple thermal traits at once:
#   - lt50, ltmin, ltmax, ctmin   (ctmax optional)
#
# Thermal sources (trait-dependent; auto-detected):
#   1) ThermoFresh/ThermoTol long table (thermtol_comb_final.csv; uses metric field)
#   2) Olden imputed (Comte_Olden_Data_Imputed.csv; scans columns)
#   3) GlobTherm (GlobalTherm_upload_02_11_17.csv; scans columns)
#
# Priority when merging trait values:
#   ThermoFresh > Olden > GlobTherm
#
# Interactions:
#   - EU metaweb edges
#   - GloBI edges built from:
#       * ThermoFresh-filtered GloBI fish metaweb
#       * Olden-filtered GloBI edges
#   Dedup inside GloBI: prefer ThermoFresh-supported edges over Olden-supported.
#
# Outputs:
#   outputs_merged_all_metrics/<trait>/
#     - edge_table.csv
#     - predator_level.csv
#     - summary.csv
#     - edge_pred_vs_prey_color_by_interactions.png
#     - edge_pred_vs_prey_color_by_trait_source.png
#     - node_pred_vs_meanprey_color_by_trait_source.png
#     - null_mean_absdiff.png
#   outputs_merged_all_metrics/summary_all_metrics.csv
# ============================================================

using CSV, DataFrames, Statistics, Random
using CairoMakie

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

# Output root
OUTROOT = joinpath(project_root, "outputs_merged_all_metrics")
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

# Try to get a column by name (case-insensitive) from a DataFrame
function find_col_ci(df::DataFrame, candidates::Vector{Symbol})
    nms = names(df)
    nms_l = lowercase.(String.(nms))
    for c in candidates
        i = findfirst(==(lowercase(String(c))), nms_l)
        i === nothing || return nms[i]
    end
    return nothing
end

# Defensive string extraction (handles missing and UndefRefError)
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

# Normalize trait names from metric strings / column names
function normalize_trait_name(raw::AbstractString)::String
    s = lowercase(strip(String(raw)))
    s2 = replace(s, r"[\s_\-]" => "")

    if occursin("ctmax", s2) || occursin("criticalthermalmaximum", s2)
        return "ctmax"
    elseif occursin("ctmin", s2) || occursin("criticalthermalminimum", s2)
        return "ctmin"
    end

    if occursin("lt50", s2)
        return "lt50"
    elseif occursin("ltmin", s2) || (occursin("lowerlethal", s2) && occursin("temp", s2))
        return "ltmin"
    elseif occursin("ltmax", s2) || (occursin("upperlethal", s2) && occursin("temp", s2))
        return "ltmax"
    end

    if s2 == "tmin"
        return "tmin"
    elseif s2 == "tmax"
        return "tmax"
    end

    return ""
end

# ============================================================
# 3) LOAD INTERACTIONS ONCE (EU + GloBI-best)
# ============================================================
# --- 3a) EU metaweb ---
eu_raw = CSV.read(eu_metaweb_csv, DataFrame; missingstring="", ntasks=1, pool=false)

eu = DataFrame(
    pred = canon_taxon.(string.(eu_raw.sourceTaxonName)),
    prey = canon_taxon.(string.(eu_raw.targetTaxonName))
)
filter!(r -> !isempty(r.pred) && !isempty(r.prey), eu)
eu.interaction_source .= "EU"
eu.evidence_source    .= "EU"

# --- 3b) GloBI edges from ThermoFresh-filtered build ---
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

# --- 3c) GloBI edges from Olden-filtered build ---
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

# --- 3d) Apply preference within GloBI: ThermoFresh > Olden ---
globi_all = vcat(g_tf, g_ol)
evidence_rank(s::AbstractString) = (s == "ThermoFresh") ? 3 : (s == "Olden") ? 2 : 1

globi_best = combine(groupby(globi_all, [:pred, :prey])) do sdf
    i = argmax(evidence_rank.(String.(sdf.evidence_source)))
    sdf[i:i, :]
end

# --- 3e) Merge EU + best GloBI, then collapse duplicates across DBs ---
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
@info "GloBI edges (after preference)" nrow(globi_best)
@info "Merged unique interactions (EU + GloBI-best)" nrow(merged_edges)

# ============================================================
# 4) LOAD THERMAL SOURCES ONCE
# ============================================================
tf_df = CSV.read(thermtol_csv, DataFrame; missingstring="", ntasks=1, pool=false)
ol_df = CSV.read(olden_csv, DataFrame; missingstring="", ntasks=1, pool=false)
gt_df = CSV.read(globtherm_csv, DataFrame; missingstring="", ntasks=1, pool=false)

# Precompute canonical species fields
tf_species_canon = canon_taxon.(string.(coalesce.(tf_df.species, "")))
tf_metric_raw    = string.(coalesce.(tf_df.metric, ""))
tf_metric_norm   = normalize_trait_name.(tf_metric_raw)
tf_vals_all       = as_float.(tf_df.tol)

ol_species_canon = canon_taxon.(string.(coalesce.(ol_df[!, :Species], "")))

# GlobTherm canonical species (try common genus/species column variants)
gt_genus_col   = find_col_ci(gt_df, [:Genus, :genus])
gt_species_col = find_col_ci(gt_df, [:Species, :species])

gt_genus_vec   = gt_genus_col === nothing ? fill("", nrow(gt_df)) : safe_string_vec(gt_df[!, gt_genus_col])
gt_species_vec = gt_species_col === nothing ? fill("", nrow(gt_df)) : safe_string_vec(gt_df[!, gt_species_col])
gt_canon_raw   = strip.(gt_genus_vec .* " " .* gt_species_vec)
gt_species_canon = canon_taxon.(gt_canon_raw)

# Candidate columns in Olden/GlobTherm by normalized trait name
function candidate_columns_by_trait(df::DataFrame)
    out = Dict{String, Vector{Symbol}}()
    for nm in names(df)
        tr = normalize_trait_name(String(nm))
        tr == "" && continue
        push!(get!(out, tr, Symbol[]), Symbol(nm))  # <-- key fix
    end
    return out
end

ol_cols_by_trait = candidate_columns_by_trait(ol_df)
gt_cols_by_trait = candidate_columns_by_trait(gt_df)

# ============================================================
# 5) BUILD LOOKUP FOR A GIVEN TRAIT (priority: TF > Olden > GlobTherm)
#     Returns (lookup::Dict{String,Float64}, source::Dict{String,String})
# ============================================================
function build_trait_lookup(trait::String;
    tf_species=tf_species_canon, tf_trait=tf_metric_norm, tf_vals=tf_vals_all,
    ol_df=ol_df, ol_species=ol_species_canon, ol_cols=ol_cols_by_trait,
    gt_df=gt_df, gt_species=gt_species_canon, gt_cols=gt_cols_by_trait
)
    # --- ThermoFresh long table ---
    tf_idx = tf_trait .== trait
    tf_tmp = DataFrame(canon=tf_species[tf_idx], val=tf_vals[tf_idx])
    tf_tmp = tf_tmp[.!isempty.(tf_tmp.canon) .& isfinite.(tf_tmp.val), :]
    tf_sum = isempty(tf_tmp) ? DataFrame(canon=String[], val=Float64[]) :
             combine(groupby(tf_tmp, :canon), :val => mean => :val)
    tf_lookup = Dict(r.canon => r.val for r in eachrow(tf_sum))

    # --- Olden (columns) ---
    ol_lookup = Dict{String,Float64}()
    if haskey(ol_cols, trait)
        cols = ol_cols[trait]
        if !isempty(cols)
            vals = fill(NaN, nrow(ol_df))
            for c in cols
                v = parse_comma_float.(ol_df[!, c])
                for i in eachindex(vals)
                    if !isfinite(vals[i]) && isfinite(v[i])
                        vals[i] = v[i]
                    end
                end
            end
            tmp = DataFrame(canon=ol_species, val=vals)
            tmp = tmp[.!isempty.(tmp.canon) .& isfinite.(tmp.val), :]
            sumdf = isempty(tmp) ? DataFrame(canon=String[], val=Float64[]) :
                    combine(groupby(tmp, :canon), :val => mean => :val)
            ol_lookup = Dict(r.canon => r.val for r in eachrow(sumdf))
        end
    end

    # --- GlobTherm (columns) ---
    gt_lookup = Dict{String,Float64}()
    if haskey(gt_cols, trait)
        cols = gt_cols[trait]
        if !isempty(cols)
            vals = fill(NaN, nrow(gt_df))
            for c in cols
                v = as_float.(gt_df[!, c])
                for i in eachindex(vals)
                    if !isfinite(vals[i]) && isfinite(v[i])
                        vals[i] = v[i]
                    end
                end
            end
            tmp = DataFrame(canon=gt_species, val=vals)
            tmp = tmp[.!isempty.(tmp.canon) .& isfinite.(tmp.val), :]
            sumdf = isempty(tmp) ? DataFrame(canon=String[], val=Float64[]) :
                    combine(groupby(tmp, :canon), :val => mean => :val)
            gt_lookup = Dict(r.canon => r.val for r in eachrow(sumdf))
        end
    end

    # Merge with priority (lowest first)
    lookup = Dict{String,Float64}()
    source = Dict{String,String}()

    for (sp,v) in gt_lookup
        lookup[sp] = v
        source[sp] = "GlobTherm"
    end
    for (sp,v) in ol_lookup
        lookup[sp] = v
        source[sp] = "Olden"
    end
    for (sp,v) in tf_lookup
        lookup[sp] = v
        source[sp] = "ThermoFresh"
    end

    # For reporting (helps confirm which DBs contributed)
    @info "Trait=$(trait): species w/ values" length(lookup) "TF" length(tf_lookup) "Olden" length(ol_lookup) "GlobTherm" length(gt_lookup)

    return lookup, source
end

# ============================================================
# 6) RUN PIPELINE FOR ALL TRAITS WITH DATA
# ============================================================
TRAITS = ["lt50", "ltmin", "ltmax", "ctmin"]

summary_all = DataFrame(
    metric=String[],
    edges=Int[],
    predators=Int[],
    prey=Int[],
    obs_mean_absdiff=Float64[],
    null_mean=Float64[],
    null_sd=Float64[],
    z=Float64[],
    p_one_sided=Float64[]
)

for trait in TRAITS
    outdir = joinpath(OUTROOT, trait)
    isdir(outdir) || mkpath(outdir)

    lookup, src = build_trait_lookup(trait)

    # Filter edges to those with trait on both ends
    edges0 = merged_edges[
        haskey.(Ref(lookup), merged_edges.pred) .&
        haskey.(Ref(lookup), merged_edges.prey),
        :
    ]

    @info "Trait=$(trait): edges with trait on both ends" nrow(edges0)

    if nrow(edges0) == 0
        # still write a summary stub so you can see it ran
        summ = DataFrame(
            metric = trait,
            edges = 0,
            predators = 0,
            prey = 0,
            obs_mean_absdiff = NaN,
            null_mean = NaN,
            null_sd = NaN,
            z = NaN,
            p_one_sided = NaN
        )
        CSV.write(joinpath(outdir, "summary.csv"), summ)
        push!(summary_all, summ[1, :])
        continue
    end

    pred_val = [lookup[p] for p in edges0.pred]
    prey_val = [lookup[p] for p in edges0.prey]
    absdiff  = abs.(pred_val .- prey_val)

    src_pred = [src[p] for p in edges0.pred]
    src_prey = [src[p] for p in edges0.prey]

    edges = DataFrame(
        pred = edges0.pred,
        prey = edges0.prey,
        interaction_sources = edges0.interaction_sources,
        evidence_sources = edges0.evidence_sources,
        pred_trait = pred_val,
        prey_trait = prey_val,
        absdiff = absdiff,
        trait_source_pred = src_pred,
        trait_source_prey = src_prey
    )
    CSV.write(joinpath(outdir, "edge_table.csv"), edges)

    # Node-level predator vs mean prey
    edges_clean = edges[isfinite.(edges.pred_trait) .& isfinite.(edges.prey_trait), :]

    pred_level = if isempty(edges_clean)
        DataFrame(pred=String[], pred_trait=Float64[], mean_prey_trait=Float64[], trait_source_pred=String[], nprey=Int[])
    else
        pred_grp = groupby(edges_clean, :pred)
        combine(pred_grp,
            :pred_trait => mean => :pred_trait,
            :prey_trait => mean => :mean_prey_trait,
            :trait_source_pred => (x -> first(x)) => :trait_source_pred,
            nrow => :nprey
        )
    end
    CSV.write(joinpath(outdir, "predator_level.csv"), pred_level)

    # Null model
    obs = isempty(edges_clean) ? NaN : mean(edges_clean.absdiff)

    μ0 = NaN; σ0 = NaN; z = NaN; p = NaN
    null_vals = Float64[]
    clean_null = Float64[]

    if nrow(edges_clean) >= 2
        base_preds = collect(edges_clean.pred)
        base_preys = collect(edges_clean.prey)
        null_vals = null_mean_absdiff(base_preds, base_preys, lookup; nperm=1000, nsweeps=10, seed=1)
        clean_null = filter(isfinite, null_vals)
        μ0 = isempty(clean_null) ? NaN : mean(clean_null)
        σ0 = isempty(clean_null) ? NaN : std(clean_null)
        z  = (σ0 > 0 && isfinite(obs)) ? (obs - μ0) / σ0 : NaN
        p  = isempty(null_vals) || !isfinite(obs) ? NaN : mean(null_vals .<= obs)  # unusually small mean |Δ|?
    end

    # Plots (same as CTmax-style; just relabeled)
    if nrow(edges_clean) >= 2
        scatter_with_fit_color(
            edges_clean.pred_trait, edges_clean.prey_trait, edges_clean.interaction_sources;
            title="Edge-level: predator vs prey ($(trait))",
            xlab="Predator $(trait)",
            ylab="Prey $(trait)",
            file=joinpath(outdir, "edge_pred_vs_prey_color_by_interactions.png")
        )

        scatter_with_fit_color(
            edges_clean.pred_trait, edges_clean.prey_trait, edges_clean.trait_source_pred;
            title="Edge-level: predator vs prey ($(trait))",
            xlab="Predator $(trait)",
            ylab="Prey $(trait)",
            file=joinpath(outdir, "edge_pred_vs_prey_color_by_trait_source.png")
        )
    end

    if nrow(pred_level) >= 2
        scatter_with_fit_color(
            pred_level.pred_trait, pred_level.mean_prey_trait, pred_level.trait_source_pred;
            title="Node-level: predator vs mean(prey) ($(trait))",
            xlab="Predator $(trait)",
            ylab="Mean prey $(trait)",
            file=joinpath(outdir, "node_pred_vs_meanprey_color_by_trait_source.png")
        )
    end

    # Null histogram
    begin
        fig = Figure(size=(900,570))
        ax = Axis(fig[1,1],
            title="Null: mean |Δ$(trait)| (pred-swap)",
            xlabel="Mean |Δ$(trait)|",
            ylabel="Count"
        )
        isempty(clean_null) || hist!(ax, clean_null, bins=30)
        isfinite(obs) && vlines!(ax, [obs])
        save(joinpath(outdir, "null_mean_absdiff.png"), fig)
        display(fig)
    end

    # Summary
    summ = DataFrame(
        metric = trait,
        edges = nrow(edges_clean),
        predators = length(unique(edges_clean.pred)),
        prey = length(unique(edges_clean.prey)),
        obs_mean_absdiff = obs,
        null_mean = μ0,
        null_sd = σ0,
        z = z,
        p_one_sided = p
    )
    CSV.write(joinpath(outdir, "summary.csv"), summ)
    push!(summary_all, summ[1, :])
end

CSV.write(joinpath(OUTROOT, "summary_all_metrics.csv"), summary_all)
@info "DONE" OUTROOT
println("\n=== summary_all_metrics ===")
show(summary_all, allrows=true, allcols=true)
println()
