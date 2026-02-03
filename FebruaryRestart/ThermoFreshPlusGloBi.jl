############################################################
# GloBI → ThermoFresh trophic metaweb builder
# Predators: fish only (Actinopteri)
# Prey: any ThermoFresh species
############################################################
using HTTP, JSON3
using DataFrames, CSV
using ProgressMeter

# =========================
# 1) USER INPUT
# =========================

# ThermoFresh DataFrame (already loaded)
tf = CSV.read(
    "FebruaryRestart/thermtol_taxonomy_final.csv",
    DataFrame;
    missingstring = "",
    ntasks = 1,
    pool = false
)

# Column names in ThermoFresh
TF_SPECIES_COL = :species
TF_CLASS_COL   = :class

# Output file
OUTFILE = "thermofresh_globi_metaweb_fish_predators.csv"

# =========================
# 2) CONSTANTS
# =========================

GLOBI_URL = "https://api.globalbioticinteractions.org/interaction"

TROPHIC_INTERACTIONS = Set([
    "eats",
    "preysOn",
    "feedsOn",
    "eatsParasiteOf"
])

# =========================
# 3) HELPERS
# =========================

function url_encode_species(sp::String)
    replace(strip(sp), " " => "%20")
end

function globi_interactions_for_species(sp::String)
    query = "?sourceTaxon=$(url_encode_species(sp))"
    r = HTTP.get(GLOBI_URL * query)
    r.status != 200 && return DataFrame()
    return DataFrame(JSON3.read(r.body))
end

# =========================
# 4) SPECIES SETS (🔒 CORE LOGIC)
# =========================

# --- ALL ThermoFresh species (allowed prey) ---
tf_species_all = Set(
    strip.(String.(skipmissing(tf[!, TF_SPECIES_COL])))
)

@info "Total ThermoFresh species" length(tf_species_all)

# --- ONLY fish species (allowed predators) ---
tf_nomiss = tf[.!ismissing.(tf.class), :]
tf_fish = tf_nomiss[tf_nomiss.class .== "Actinopteri", :]

fish_predators = Set(
    strip.(String.(skipmissing(tf_fish[!, TF_SPECIES_COL])))
)

@info "Fish predators (Actinopteri)" length(fish_predators)

# =========================
# 5) MAIN LOOP (🔒 PREDATORS = FISH ONLY)
# =========================
edges = Vector{Tuple{String,String}}()

@showprogress for sp in fish_predators   # 🔒 iterate ONLY over fish

    df = globi_interactions_for_species(sp)
    isempty(df) && continue

    # Required columns safety check
    required = (:source_taxon_name, :target_taxon_name, :interaction_type)
    all(c -> c in propertynames(df), required) || continue

    # Keep trophic interactions only
    df = df[in.(df.interaction_type, Ref(TROPHIC_INTERACTIONS)), :]
    isempty(df) && continue

    for r in eachrow(df)
        predator = strip(r.source_taxon_name)
        prey     = strip(r.target_taxon_name)

        # 🔒 HARD CONSTRAINTS
        predator ∈ fish_predators    || continue   # predator must be fish
        prey     ∈ tf_species_all    || continue   # prey must be in ThermoFresh

        push!(edges, (predator, prey))
    end
end

@info "Raw trophic edges collected" length(edges)

# =========================
# 6) BUILD META-WEB
# =========================
metaweb = unique(DataFrame(edges, [:predator, :prey]))

@info "Final metaweb summary" (
    edges      = nrow(metaweb),
    predators  = length(unique(metaweb.predator)),
    prey       = length(unique(metaweb.prey))
)

# =========================
# 7) SANITY CHECKS (DO NOT SKIP)
# =========================
@info "Non-fish predators (should be empty)" setdiff(
    unique(metaweb.predator),
    fish_predators
)

# =========================
# 8) SAVE OUTPUT
# =========================
CSV.write(OUTFILE, metaweb)
@info "Saved metaweb to" OUTFILE

############################################################
# END
############################################################
