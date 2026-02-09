#!/usr/bin/env julia
using CSV, DataFrames
using HTTP, JSON3
using ProgressMeter

# =========================
# 0) CANONICAL TAXON
# =========================
function canon_taxon(s::AbstractString)::String
    t = replace(String(s), r"[_(),\[\]]" => " ")
    t = replace(t, r"\s+" => " ")
    t = strip(t)

    parts = split(t)
    length(parts) < 2 && return ""

    genus_raw   = parts[1]
    species_raw = parts[2]

    genus   = String(filter(isletter, genus_raw))
    species = String(filter(isletter, species_raw))

    # FIX: must return if either is empty
    (isempty(genus) || isempty(species)) && return ""

    lowercase(species) in ("sp","spp") && return ""

    return uppercase(genus[1]) * lowercase(genus[2:end]) * " " * lowercase(species)
end

# =========================
# 1) LOAD IMPUTED DATA
# =========================
imputed_df = CSV.read(
    "FebruaryRestart/Comte_Olden_Data_Imputed.csv",
    DataFrame
)

parse_comma_float(x) =
    x === missing ? NaN :
    x isa Real ? Float64(x) :
    x isa AbstractString ? parse(Float64, replace(x, "," => ".")) :
    NaN

imputed_df = transform(
    imputed_df,
    :Species => (s -> canon_taxon.(coalesce.(String.(s), ""))) => :canon,
    Symbol("CTmax") => (v -> parse_comma_float.(v)) => :CTmax
)

imputed_df = imputed_df[
    (.!isempty.(imputed_df.canon)) .&
    isfinite.(imputed_df.CTmax),
    :
]

ct_lookup = Dict(r.canon => r.CTmax for r in eachrow(imputed_df))
@info "Imputed CTmax species" length(ct_lookup)

# =========================
# 2) GLOBI HELPERS
# =========================
GLOBI_URL = "https://api.globalbioticinteractions.org/interaction"

TROPHIC = Set([
    "eats",
    "preysOn",
    "feedsOn",
    "eatsParasiteOf"
])

function globi_for_species(sp::String)::DataFrame
    # Safer than manual space replacement, and fixes your url_encode/urlencode mismatch
    q = "?sourceTaxon=" * HTTP.escapeuri(strip(sp))

    r = HTTP.get(GLOBI_URL * q; status_exception=false)
    r.status == 200 || return DataFrame()

    obj = JSON3.read(r.body)

    (haskey(obj, :columns) && haskey(obj, :data)) || return DataFrame()

    cols = Symbol.(obj.columns)
    rows = obj.data
    isempty(rows) && return DataFrame()

    # Build row-wise (same fix as before)
    nt_rows = NamedTuple[]
    for row in rows
        length(row) == length(cols) || continue
        push!(nt_rows, NamedTuple{Tuple(cols)}(Tuple(row)))
    end

    isempty(nt_rows) && return DataFrame()
    return DataFrame(nt_rows)
end

# =========================
# 3) QUERY GLOBI
# =========================
edges = Vector{Tuple{String,String}}()

species = collect(keys(ct_lookup))

@showprogress for sp in species
    df = globi_for_species(sp)
    isempty(df) && continue

    required = (:source_taxon_name, :target_taxon_name, :interaction_type)
    all(c -> c in propertynames(df), required) || continue

    df = df[in.(df.interaction_type, Ref(TROPHIC)), :]
    isempty(df) && continue

    for r in eachrow(df)
        # Guard against missing values
        (ismissing(r.source_taxon_name) || ismissing(r.target_taxon_name)) && continue

        pred = canon_taxon(String(r.source_taxon_name))
        prey = canon_taxon(String(r.target_taxon_name))

        isempty(pred) && continue
        isempty(prey) && continue

        haskey(ct_lookup, pred) || continue
        haskey(ct_lookup, prey) || continue

        push!(edges, (pred, prey))
    end
end

# =========================
# 4) RESULTS
# =========================
edges_df = unique(DataFrame(edges, [:pred, :prey]))

@info "Matched GloBI edges (CTmax on both ends)" nrow(edges_df)
@info "Predators with edges" length(unique(edges_df.pred))
@info "Prey with edges" length(unique(edges_df.prey))

CSV.write(
    "outputs_imputed_globi_edges.csv",
    edges_df
)

@info "Saved matched edges to outputs_imputed_globi_edges.csv"
