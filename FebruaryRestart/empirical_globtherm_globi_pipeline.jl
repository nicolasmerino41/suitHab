#!/usr/bin/env julia
# ============================================================
# GloBI × GlobTherm linkage diagnostic
# PURPOSE:
#   • Pull trophic interactions from GloBI
#   • Match BOTH predator and prey against GlobTherm
#   • Report how many links survive
#
# OUTPUT:
#   • counts only (no plots)
# ============================================================

using CSV, DataFrames, Statistics
using HTTP, JSON3
using ProgressMeter

# ============================================================
# 0) PATHS
# ============================================================
script_dir   = @__DIR__
project_root = script_dir

globtherm_csv = joinpath(project_root, "GlobalTherm_upload_02_11_17.csv")  # <-- adjust if needed
out_csv       = joinpath(project_root, "globi_globtherm_links.csv")

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

    # FIX: correct boolean logic
    (isempty(genus) || isempty(species)) && return ""
    lowercase(species) in ("sp","spp") && return ""

    return uppercase(genus[1]) * lowercase(genus[2:end]) * " " * lowercase(species)
end

# ============================================================
# 2) LOAD GLOBTHERM
# ============================================================
gt = CSV.read(globtherm_csv, DataFrame; missingstring="", ntasks=1, pool=false)

gt = transform(gt,
    :Species => (col -> canon_taxon.(string.(coalesce.(col, "")))) => :canon
)

gt = gt[.!isempty.(gt.canon), :]

globtherm_species = Set(gt.canon)
@info "GlobTherm species loaded" length(globtherm_species)

# ============================================================
# 3) GLOBI QUERY SETUP
# ============================================================
GLOBI_URL = "https://api.globalbioticinteractions.org/interaction"
TROPHIC = Set(["eats","preysOn","feedsOn","eatsParasiteOf"])

function globi_for_species(sp::String)::DataFrame
    q = "?sourceTaxon=" * HTTP.escapeuri(strip(sp))

    r = HTTP.get(GLOBI_URL * q; status_exception=false)
    r.status == 200 || return DataFrame()

    obj = JSON3.read(r.body)
    (haskey(obj, :columns) && haskey(obj, :data)) || return DataFrame()

    cols = Symbol.(obj.columns)
    rows = obj.data
    isempty(rows) && return DataFrame()

    # SAME FIX AS BEFORE: rows -> NamedTuples -> DataFrame
    nt_rows = NamedTuple[]
    for row in rows
        length(row) == length(cols) || continue
        push!(nt_rows, NamedTuple{Tuple(cols)}(Tuple(row)))
    end

    isempty(nt_rows) && return DataFrame()
    return DataFrame(nt_rows)
end

# ============================================================
# 4) QUERY GLOBI (PREDATORS = GlobTherm species)
# ============================================================
edges = Vector{Tuple{String,String}}()

@showprogress for sp in globtherm_species
    df = globi_for_species(sp)
    isempty(df) && continue

    required = (:source_taxon_name, :target_taxon_name, :interaction_type)
    all(c -> c in propertynames(df), required) || continue

    df = df[in.(df.interaction_type, Ref(TROPHIC)), :]
    isempty(df) && continue

    for r in eachrow(df)
        (ismissing(r.source_taxon_name) || ismissing(r.target_taxon_name)) && continue

        pred = canon_taxon(String(r.source_taxon_name))
        prey = canon_taxon(String(r.target_taxon_name))

        isempty(pred) && continue
        isempty(prey) && continue

        pred ∈ globtherm_species || continue
        prey ∈ globtherm_species || continue

        push!(edges, (pred, prey))
    end
end

edges_df = unique(DataFrame(edges, [:predator, :prey]))
CSV.write(out_csv, edges_df)

# ============================================================
# 5) SUMMARY
# ============================================================
@info "RESULT SUMMARY" (
    globtherm_species = length(globtherm_species),
    matched_edges     = nrow(edges_df),
    predators         = length(unique(edges_df.predator)),
    prey              = length(unique(edges_df.prey))
)

println("\n=== DONE ===")
