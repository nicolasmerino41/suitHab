using CSV, DataFrames
using HTTP, JSON3, DataFrames

globi = CSV.read(
    "FebruaryRestart/interactions.tsv",
    DataFrame;
    delim = '\t',
    missingstring = ""
)

# Keep only fish–fish trophic interactions
globi_fish = globi[
    (globi.interactionType .∈ ["eats", "preysOn"]) .&
    (globi.taxonClass .== "Actinopterygii") .&
    (globi.targetTaxonClass .== "Actinopterygii"),
    [:taxonName, :targetTaxonName]
]

rename!(globi_fish,
    :taxonName => :consumer,
    :targetTaxonName => :prey
)

# Drop duplicates
unique!(globi_fish)

function globi_interactions_for_species(sp::AbstractString)
    base = "https://api.globalbioticinteractions.org/interaction"
    q = HTTP.escapeuri(sp)
    r = HTTP.get("$base?sourceTaxon=$q&limit=1000")

    r.status != 200 && return DataFrame()

    obj = JSON3.read(r.body)

    cols = Vector{String}(obj["columns"])
    rows = obj["data"]

    # Convert rows → DataFrame
    df = DataFrame([ [row[i] for row in rows] for i in eachindex(cols) ],
                    Symbol.(cols))

    return df
end

TROPHIC = Set([
    "eats",
    "preysOn",
    "eatsParasiteOf",
    "feedsOn"
])

df_trophic = df[in.(df.interaction_type, Ref(TROPHIC)), :]
