#!/usr/bin/env julia
# ============================================================
# Thermal trait inventory + counts across sources
#
# Sources:
#   1) ThermoFresh/ThermoTol long table (thermtol_comb_final.csv)
#      - expects columns like: species, metric, tol
#   2) Olden imputed (Comte_Olden_Data_Imputed.csv)
#      - scans column names for thermal traits
#   3) GlobTherm (GlobalTherm_upload_02_11_17.csv)
#      - scans column names for thermal traits
#
# Outputs (in outputs_thermal_trait_inventory/):
#   - tf_unique_metrics.csv
#   - olden_candidate_columns.csv
#   - globtherm_candidate_columns.csv
#   - trait_counts_by_source.csv
#   - trait_counts_overall.csv
#   - trait_species_overlap.csv
# ============================================================
using CSV, DataFrames, Statistics

# =========================
# 0) PATHS (EDIT THESE)
# =========================
script_dir   = @__DIR__
project_root = script_dir

thermtol_csv  = joinpath(project_root, "thermtol_comb_final.csv")
olden_csv     = joinpath(project_root, "Comte_Olden_Data_Imputed.csv")
globtherm_csv = joinpath(project_root, "GlobalTherm_upload_02_11_17.csv")

OUTROOT = joinpath(project_root, "outputs_thermal_trait_inventory")
isdir(OUTROOT) || mkpath(OUTROOT)

# =========================
# 1) CANONICAL TAXON
# =========================
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

# =========================
# 2) VALUE PARSERS
# =========================
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

# =========================
# 3) TRAIT NAME NORMALIZATION
# =========================
"""
Map a raw metric/column name to a standard trait key.
Returns "" if not recognized.
"""
function normalize_trait_name(raw::AbstractString)::String
    s = lowercase(strip(String(raw)))

    # remove common separators
    s2 = replace(s, r"[\s_\-]" => "")

    # CTmax / CTmin
    if occursin("ctmax", s2) || occursin("criticalthermalmaximum", s2)
        return "ctmax"
    elseif occursin("ctmin", s2) || occursin("criticalthermalminimum", s2)
        return "ctmin"
    end

    # lethal temps (LT50, LTmin, LTmax) — allow variants like "lt_50", "lt50c", etc.
    if occursin("lt50", s2)
        return "lt50"
    elseif occursin("ltmin", s2)
        return "ltmin"
    elseif occursin("ltmax", s2)
        return "ltmax"
    end

    # generic Tmin/Tmax (many datasets use these)
    # NOTE: we keep these distinct from CTmin/CTmax
    if s2 == "tmin" || occursin(r"^tmin[^a-z0-9]*"i, s)
        return "tmin"
    elseif s2 == "tmax" || occursin(r"^tmax[^a-z0-9]*"i, s)
        return "tmax"
    end

    # sometimes "lowerlethaltemperature" / "upperlethaltemperature"
    if occursin("lowerlethal", s2) && occursin("temp", s2)
        return "ltmin"
    elseif occursin("upperlethal", s2) && occursin("temp", s2)
        return "ltmax"
    end

    return ""
end

# Standard traits you explicitly care about (we’ll still record others if they map)
STANDARD_TRAITS = ["tmin", "lt50", "ltmin", "ltmax", "ctmin", "ctmax", "tmax"]

# =========================
# 4) COUNTING HELPERS
# =========================
"""
Given vectors of canon species and values, return (n_records, n_species).
"""
function count_records_species(canon::Vector{String}, vals::Vector{Float64})
    good = .!isempty.(canon) .& isfinite.(vals)
    n_records = sum(good)
    n_species = length(unique(canon[good]))
    return n_records, n_species
end

"""
Append a count row to a results DataFrame.
"""
function push_count!(df::DataFrame; source::String, trait::String, n_records::Int, n_species::Int)
    push!(df, (source=source, trait=trait, n_records=n_records, n_species=n_species))
end

# ============================================================
# 5) THERMOFRESH / THERMOTOL (LONG TABLE): scan unique metrics + counts
# ============================================================
tf_df = CSV.read(thermtol_csv, DataFrame; missingstring="", ntasks=1, pool=false)

# Save discovered metrics (raw)
tf_metric_raw = string.(coalesce.(tf_df.metric, ""))
tf_unique = sort(unique(tf_metric_raw))
CSV.write(joinpath(OUTROOT, "tf_unique_metrics.csv"), DataFrame(metric=tf_unique))

# Build counts for mapped traits
tf_species = canon_taxon.(string.(coalesce.(tf_df.species, "")))
tf_vals    = as_float.(tf_df.tol)

tf_trait_norm = normalize_trait_name.(tf_metric_raw)

tf_counts = DataFrame(source=String[], trait=String[], n_records=Int[], n_species=Int[])
tf_trait_species_sets = Dict{Tuple{String,String}, Set{String}}()  # (source,trait) => Set(species)

for tr in STANDARD_TRAITS
    idx = tf_trait_norm .== tr
    canon = tf_species[idx]
    vals  = tf_vals[idx]
    nrec, nsp = count_records_species(canon, vals)
    push_count!(tf_counts; source="ThermoFresh/ThermoTol", trait=tr, n_records=nrec, n_species=nsp)

    good = .!isempty.(canon) .& isfinite.(vals)
    tf_trait_species_sets[("ThermoFresh/ThermoTol", tr)] = Set(canon[good])
end

# ============================================================
# 6) OLDEN: scan column names, pick candidate thermal columns, counts
# ============================================================
ol_df = CSV.read(olden_csv, DataFrame; missingstring="", ntasks=1, pool=false)

# Canon species column (assumes :Species like your CTmax pipeline)
ol_canon = canon_taxon.(string.(coalesce.(ol_df[!, :Species], "")))

# Identify candidate columns by attempting normalization
ol_candidates = DataFrame(column=String[], normalized_trait=String[])
for nm in names(ol_df)
    tr = normalize_trait_name(String(nm))
    if tr != ""
        push!(ol_candidates, (column=String(nm), normalized_trait=tr))
    end
end
CSV.write(joinpath(OUTROOT, "olden_candidate_columns.csv"), ol_candidates)

ol_counts = DataFrame(source=String[], trait=String[], n_records=Int[], n_species=Int[])
ol_trait_species_sets = Dict{Tuple{String,String}, Set{String}}()

for tr in STANDARD_TRAITS
    # all columns that map to this trait (sometimes there are multiple)
    cols = ol_candidates.column[ol_candidates.normalized_trait .== tr]
    if isempty(cols)
        push_count!(ol_counts; source="Olden", trait=tr, n_records=0, n_species=0)
        ol_trait_species_sets[("Olden", tr)] = Set{String}()
        continue
    end

    # combine multiple columns by taking first finite among them per row
    vals = fill(NaN, nrow(ol_df))
    for c in cols
        v = parse_comma_float.(ol_df[!, Symbol(c)])
        for i in eachindex(vals)
            if !isfinite(vals[i]) && isfinite(v[i])
                vals[i] = v[i]
            end
        end
    end

    nrec, nsp = count_records_species(ol_canon, vals)
    push_count!(ol_counts; source="Olden", trait=tr, n_records=nrec, n_species=nsp)

    good = .!isempty.(ol_canon) .& isfinite.(vals)
    ol_trait_species_sets[("Olden", tr)] = Set(ol_canon[good])
end

# ============================================================
# 7) GLOBTHERM: scan columns, counts (flexible genus/species + trait columns)
# ============================================================
gt_df = CSV.read(globtherm_csv, DataFrame; missingstring="", ntasks=1, pool=false)

# Robust genus/species extraction: try common names
function getcol(df::DataFrame, sym::Symbol)
    sym in names(df) ? df[!, sym] : fill(missing, nrow(df))
end

gt_genus   = string.(coalesce.(getcol(gt_df, :Genus), ""))
gt_species = string.(coalesce.(getcol(gt_df, :Species), ""))

gt_canon_raw = strip.(gt_genus .* " " .* gt_species)
gt_canon = canon_taxon.(gt_canon_raw)

# Candidate columns
gt_candidates = DataFrame(column=String[], normalized_trait=String[])
for nm in names(gt_df)
    tr = normalize_trait_name(String(nm))
    if tr != ""
        push!(gt_candidates, (column=String(nm), normalized_trait=tr))
    end
end
CSV.write(joinpath(OUTROOT, "globtherm_candidate_columns.csv"), gt_candidates)

gt_counts = DataFrame(source=String[], trait=String[], n_records=Int[], n_species=Int[])
gt_trait_species_sets = Dict{Tuple{String,String}, Set{String}}()

for tr in STANDARD_TRAITS
    cols = gt_candidates.column[gt_candidates.normalized_trait .== tr]
    if isempty(cols)
        push_count!(gt_counts; source="GlobTherm", trait=tr, n_records=0, n_species=0)
        gt_trait_species_sets[("GlobTherm", tr)] = Set{String}()
        continue
    end

    vals = fill(NaN, nrow(gt_df))
    for c in cols
        v = as_float.(gt_df[!, Symbol(c)])
        for i in eachindex(vals)
            if !isfinite(vals[i]) && isfinite(v[i])
                vals[i] = v[i]
            end
        end
    end

    nrec, nsp = count_records_species(gt_canon, vals)
    push_count!(gt_counts; source="GlobTherm", trait=tr, n_records=nrec, n_species=nsp)

    good = .!isempty.(gt_canon) .& isfinite.(vals)
    gt_trait_species_sets[("GlobTherm", tr)] = Set(gt_canon[good])
end

# ============================================================
# 8) MERGE COUNTS + OVERALL TOTALS + SPECIES OVERLAPS
# ============================================================
counts_by_source = vcat(tf_counts, ol_counts, gt_counts)
CSV.write(joinpath(OUTROOT, "trait_counts_by_source.csv"), counts_by_source)

# Overall totals (sum records, union species across sources)
overall = DataFrame(trait=String[], total_records=Int[], total_species=Int[])
overlap_rows = DataFrame(trait=String[],
    species_TF=Int[], species_Olden=Int[], species_GlobTherm=Int[],
    union_species=Int[], intersection_all3=Int[])

for tr in STANDARD_TRAITS
    s_tf = get(tf_trait_species_sets, ("ThermoFresh/ThermoTol", tr), Set{String}())
    s_ol = get(ol_trait_species_sets, ("Olden", tr), Set{String}())
    s_gt = get(gt_trait_species_sets, ("GlobTherm", tr), Set{String}())

    union_s = union(s_tf, s_ol, s_gt)
    inter3  = intersect(s_tf, s_ol, s_gt)

    tot_records = sum(counts_by_source.n_records[counts_by_source.trait .== tr])
    push!(overall, (trait=tr, total_records=tot_records, total_species=length(union_s)))

    push!(overlap_rows, (
        trait=tr,
        species_TF=length(s_tf),
        species_Olden=length(s_ol),
        species_GlobTherm=length(s_gt),
        union_species=length(union_s),
        intersection_all3=length(inter3)
    ))
end

CSV.write(joinpath(OUTROOT, "trait_counts_overall.csv"), overall)
CSV.write(joinpath(OUTROOT, "trait_species_overlap.csv"), overlap_rows)

@info "DONE: wrote outputs to" OUTROOT
println("\n=== QUICK VIEW: overall trait totals ===")
show(overall, allrows=true, allcols=true)
println("\n\n=== QUICK VIEW: species overlap by trait ===")
show(overlap_rows, allrows=true, allcols=true)
println()
