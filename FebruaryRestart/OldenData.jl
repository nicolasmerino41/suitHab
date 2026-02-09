using XLSX, DataFrames, CSV

xlsx_path = "FebruaryRestart/Comte_Olden_Data_Imputed.csv"
imputed_df = CSV.read(xlsx_path, DataFrame)
original_df = CSV.read("FebruaryRestart/Comte_Olden_Data_Original.csv", DataFrame)
# helper: convert "38,36" -> 38.36
parse_comma_float(x) =
    x === missing ? NaN :
    x isa AbstractString ? parse(Float64, replace(x, "," => ".")) :
    NaN

# canonicalize species names
imputed_df = transform(imputed_df,
    :Species => (s -> canon_taxon.(String.(s))) => :canon,
    Symbol("CTmax") => (v -> parse_comma_float.(v)) => :CTmax
)
# canonicalize species
orig = transform(
    original_df,
    :Species => (s -> canon_taxon.(String.(s))) => :canon
)
# drop failures
orig = orig[
    (.!isempty.(orig.canon)),
    :
]

@info "Rows after canonicalization" nrow(orig)
@info "Unique species (raw)" length(unique(original_df.Species))
@info "Unique species (canonical)" length(unique(orig.canon))

# keep only valid rows
imputed_df = imputed_df[
    (.!isempty.(imputed_df.canon)) .&
    isfinite.(imputed_df.CTmax),
    :
]
@info "Comte & Olden imputed CTmax"
@info "Species with CTmax" nrow(imputed_df)

ctmax_lookup = Dict(
    r.canon => r.CTmax for r in eachrow(imputed_df)
)

web_species = union(web.pred, web.prey)
web_species = Set(web_species)

@info "Total species in web" length(web_species)

n_species_ctmax = count(sp -> haskey(ctmax_lookup, sp), web_species)

@info "Species with CTmax" n_species_ctmax
@info "Species coverage (%)" 100 * n_species_ctmax / length(web_species)

