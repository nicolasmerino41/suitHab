using CSV, DataFrames, Statistics

# -----------------------------
# SETTINGS
# -----------------------------
path = "Metawebs/BajaCaliforniaMetaweb.txt"

use_projected_metaweb = true   # true = all links; false = only PresentAtAnyEstuary == 1
use_node_ids = true            # true = ConsumerNodeID/ResourceNodeID; false = ConsumerSpeciesID/ResourceSpeciesID
allow_self = false             # include i->i links?

# Optional: filter specific LinkTypeIDs (leave empty to keep all)
# Example: keep only LinkTypeID in [1,2,3]
keep_linktype_ids = Int[]      # e.g. [1,2,3]

# -----------------------------
# Helpers
# -----------------------------
stripstr(x) = strip(string(x))

function is_one(x)
    # robustly interpret 1/0, true/false, "1"/"0", "yes"/"no"
    if x === missing
        return false
    end
    s = lowercase(stripstr(x))
    return s in ("1", "true", "t", "yes", "y")
end

# -----------------------------
# Load (tab-delimited)
# -----------------------------
df = CSV.read(path, DataFrame; delim='\t', ignorerepeated=true, missingstring="")

# -----------------------------
# Choose which links to include
# -----------------------------
if !use_projected_metaweb
    @assert :PresentAtAnyEstuary in propertynames(df) "Missing column PresentAtAnyEstuary"
    df = df[[is_one(v) for v in df.PresentAtAnyEstuary], :]
end

if !isempty(keep_linktype_ids)
    @assert :LinkTypeID in propertynames(df) "Missing column LinkTypeID"
    # LinkTypeID might be read as Int, Float, or String → normalize to Int when possible
    lt = Vector{Int}(undef, nrow(df))
    for i in 1:nrow(df)
        v = df.LinkTypeID[i]
        lt[i] = v === missing ? -1 : parse(Int, stripstr(v))
    end
    mask = [lt[i] in keep_linktype_ids for i in eachindex(lt)]
    df = df[mask, :]
end

# -----------------------------
# Extract endpoints (nodes)
# -----------------------------
if use_node_ids
    @assert (:ConsumerNodeID in propertynames(df)) && (:ResourceNodeID in propertynames(df)) "Missing ConsumerNodeID/ResourceNodeID"
    consumers = stripstr.(df.ConsumerNodeID)
    resources = stripstr.(df.ResourceNodeID)
else
    @assert (:ConsumerSpeciesID in propertynames(df)) && (:ResourceSpeciesID in propertynames(df)) "Missing ConsumerSpeciesID/ResourceSpeciesID"
    consumers = stripstr.(df.ConsumerSpeciesID)
    resources = stripstr.(df.ResourceSpeciesID)
end

# Drop empty endpoints (just in case)
keep = .!(isempty.(consumers) .| isempty.(resources))
consumers = consumers[keep]
resources = resources[keep]

# Drop self-links if desired
if !allow_self
    mask = consumers .!= resources
    consumers = consumers[mask]
    resources = resources[mask]
end

# Deduplicate realized links
links = DataFrame(resource = resources, consumer = consumers)
unique!(links)

# -----------------------------
# Core metrics
# -----------------------------
nodes = unique(vcat(links.resource, links.consumer))
S = length(nodes)
L = nrow(links)

Lpossible = allow_self ? S^2 : S * (S - 1)        # directed possible links
Connectance = Lpossible == 0 ? NaN : L / Lpossible

# "preys per node" = in-degree across ALL nodes (include zeros)
indeg_map = Dict{String, Int}()
tmp = combine(groupby(links, :consumer), :resource => (x -> length(unique(x))) => :k)
for i in 1:nrow(tmp)
    indeg_map[tmp.consumer[i]] = tmp.k[i]
end
indeg_all = [get(indeg_map, node, 0) for node in nodes]

Mean_prey_per_node = mean(Float64.(indeg_all))  # equals L/S if all nodes counted consistently
SD_prey_per_node   = std(Float64.(indeg_all))

# -----------------------------
# Print summary + compilation row
# -----------------------------
println("=== Baja California estuaries metaweb summary ===")
println("Using projected metaweb?       ", use_projected_metaweb)
println("Using node IDs (stages)?       ", use_node_ids)
println("Nodes (S):                     ", S)
println("Links (L, unique):             ", L)
println("Connectance L/(S*(S-1)):        ", Connectance)
println("Mean prey per node (in-degree): ", Mean_prey_per_node)
println("SD prey per node:               ", SD_prey_per_node)
println()

# Ready-to-paste TSV row (blanks where unknown)
Name  = "BajaCaliforniaMetaweb"
Ecosystem = "California/Baja California estuaries"
Taxa = "parasite-inclusive food web"
TypeOf = "directed trophic network"
Levels = ""              # unknown -> blank
Region = ""              # unknown -> blank
Temporal_range = ""      # unknown -> blank
HighestTaxRes = ""       # unknown -> blank
Reference = "Hechinger et al. 2011 (Ecological Archives E092-066), doi:10.1890/10-1383.1"

println("Name\tEcosystem\tTaxa\tType\tLevels\tRegion\tTemporal range\tNodes\tLinks\tConnectance\tMean degree\tSD degree\tHighest Taxonomic Resolution\tReference")
println(join([
    Name, Ecosystem, Taxa, TypeOf, Levels, Region, Temporal_range,
    string(S), string(L),
    string(Connectance),
    string(Mean_prey_per_node),
    string(SD_prey_per_node),
    HighestTaxRes,
    Reference
], "\t"))
