using CSV, DataFrames

gt = CSV.read(
    "FebruaryRestart/GlobalTherm_upload_02_11_17.csv",
    DataFrame;
    missingstring = "",
    ntasks = 1,
    pool = false
)
# Force fully-assigned Vector{Union{Missing,String}} storage (kills #undef, SubString, pooled weirdness)
function force_string_union!(df::DataFrame, cols::Vector{Symbol})
    for col in cols
        col ∈ propertynames(df) || continue
        v = df[!, col]
        n = length(v)

        newv = Vector{Union{Missing,String}}(undef, n)
        for i in 1:n
            if isassigned(v, i)
                x = v[i]
                if x === missing
                    newv[i] = missing
                elseif x isa AbstractString
                    newv[i] = String(x)  # materialize SubString etc.
                elseif x isa AbstractVector{<:AbstractString} && !isempty(x)
                    newv[i] = String(x[1])  # scalarize vector-of-strings cells
                else
                    newv[i] = missing
                end
            else
                newv[i] = missing  # <- converts #undef to missing
            end
        end

        df[!, col] = newv
    end
    return df
end

force_string_union!(gt, [:Genus, :Species, :Class, :Order, :Family])

# ==================================================
# Paths (anchored to FebruaryRestart/)
# ==================================================

# Directory of THIS script: suitHab/FebruaryRestart
script_dir = @__DIR__

# Project root: suitHab/
project_root = dirname(script_dir)

# Data directory
data_dir = joinpath(project_root, "data")

pairwise_file = joinpath(data_dir, "TetraEU_pairwise_interactions.csv")
diet_file     = joinpath(data_dir, "TetraEU_generic_diet.csv")

# ==================================================
# Load pairwise predator–prey interactions
# ==================================================

@info "Loading pairwise interactions" pairwise_file

raw_web = CSV.read(pairwise_file, DataFrame)

required_cols = ("sourceTaxonName", "targetTaxonName")
missingg = setdiff(required_cols, names(raw_web))
isempty(missingg) || error("Missing required columns in pairwise file: $missing")

web = DataFrame(
    predator = string.(raw_web.sourceTaxonName),
    prey     = string.(raw_web.targetTaxonName)
)

# --------------------------------------------------
# Species summaries
# --------------------------------------------------
unique_predators = unique(web.predator)
unique_preys     = unique(web.prey)
unique_species   = unique(vcat(unique_predators, unique_preys))

@info "Metaweb summary" (
    interactions = nrow(web),
    predators    = length(unique_predators),
    preys        = length(unique_preys),
    species      = length(unique_species)
)

# --------------------------------------------------
# Utility: check if a species is a predator
# --------------------------------------------------

function is_predator(species::AbstractString, web::DataFrame)
    any(web.predator .== species)
end

# Example:
# is_predator("Pelobates cultripes", web)

# ==================================================
# Load generic diet data
# ==================================================

@info "Loading generic diet data" diet_file

raw_diets = CSV.read(diet_file, DataFrame)

required_cols = ("sourceTaxonName", "targetGenericItemName")
missingg = setdiff(required_cols, names(raw_diets))
isempty(missingg) || error("Missing required columns in diet file: $missing")

diets = DataFrame(
    species = string.(raw_diets.sourceTaxonName),
    item    = string.(raw_diets.targetGenericItemName)
)

@info "Diet summary" (
    rows     = nrow(diets),
    species  = length(unique(diets.species))
)

# ==================================================
# Return objects for downstream scripts
# ==================================================
(; web, diets, unique_species)

