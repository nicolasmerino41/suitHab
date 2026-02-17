using CSV, DataFrames, Statistics

# -----------------------------
# USER SETTINGS
# -----------------------------
path = "Metawebs/CaliforniaRockyIntertidalMetaweb.csv"

allow_self = false           # include i->i links?
directed = true              # food web is directed
include_parasitic = true     # set false to drop parasitic==true rows (if column is usable)
only_trophic = false         # set true to keep only interactionType == "trophic" (edit value as needed)

# -----------------------------
# Helpers
# -----------------------------
stripstr(x) = strip(string(x))

function choose_node_id(df::DataFrame, num_col::Symbol, name_col::Symbol, fallback_col::Symbol)
    if num_col in propertynames(df)
        v = df[!, num_col]
        return [v[i] === missing ? "" : strip(string(v[i])) for i in eachindex(v)]
    elseif name_col in propertynames(df)
        v = df[!, name_col]
        return [v[i] === missing ? "" : strip(string(v[i])) for i in eachindex(v)]
    else
        v = df[!, fallback_col]
        return [v[i] === missing ? "" : strip(string(v[i])) for i in eachindex(v)]
    end
end

function nonempty_pair(a, b)
    keep = trues(length(a))
    @inbounds for i in eachindex(a, b)
        keep[i] = !(isempty(a[i]) || isempty(b[i]))
    end
    return keep
end

# -----------------------------
# Load
# -----------------------------
df = CSV.read(path, DataFrame)

# Optional filters
if !include_parasitic && (:parasitic in propertynames(df))
    # Attempt to interpret parasitic as Bool-ish
    par = df.parasitic
    mask = trues(nrow(df))
    for i in eachindex(par)
        v = par[i]
        s = v === missing ? "" : lowercase(strip(String(v)))
        # treat these as "true"
        istrue = s in ["true", "t", "1", "yes", "y"]
        mask[i] = !istrue
    end
    df = df[mask, :]
end

if only_trophic && (:interactionType in propertynames(df))
    # You may need to change "trophic" to whatever label the dataset uses.
    df = df[stripstr.(df.interactionType) .== "trophic", :]
end

# Choose endpoints (prefer Num columns; else Name; else resource/consumer)
res_id = choose_node_id(df, :resourceNum, :resourceName, :resource)
con_id = choose_node_id(df, :consumerNum, :consumerName, :consumer)

# Drop rows with empty endpoints
keep = nonempty_pair(res_id, con_id)
res_id = res_id[keep]
con_id = con_id[keep]

# Remove self-links if desired
if !allow_self
    mask = res_id .!= con_id
    res_id = res_id[mask]
    con_id = con_id[mask]
end

# Deduplicate realized links
links = DataFrame(resource = res_id, consumer = con_id)
unique!(links)

S = length(unique(vcat(links.resource, links.consumer)))  # Nodes
L = nrow(links)                                           # Links

# Connectance (directed)
Lpossible = directed ? (allow_self ? S^2 : S*(S-1)) :
                       (allow_self ? (S*(S+1)) ÷ 2 : (S*(S-1)) ÷ 2)
Connectance = Lpossible == 0 ? NaN : L / Lpossible

# -----------------------------
# Degree stats
# -----------------------------
# In-degree = number of distinct resources (prey) per consumer node
indeg = combine(groupby(links, :consumer), :resource => (x -> length(unique(x))) => :k).k
# Out-degree = number of distinct consumers per resource node
outdeg = combine(groupby(links, :resource), :consumer => (x -> length(unique(x))) => :k).k

# Ensure we include nodes with zero degree in each distribution
# (important for SD). Build maps then fill for all nodes.
nodes = unique(vcat(links.resource, links.consumer))

indeg_map = Dict{String, Int}()
for (node, k) in zip(combine(groupby(links, :consumer), nrow => :tmp).consumer,
                     indeg)
    indeg_map[node] = k
end

outdeg_map = Dict{String, Int}()
for (node, k) in zip(combine(groupby(links, :resource), nrow => :tmp).resource,
                     outdeg)
    outdeg_map[node] = k
end

indeg_all = [get(indeg_map, node, 0) for node in nodes]
outdeg_all = [get(outdeg_map, node, 0) for node in nodes]
totaldeg_all = indeg_all .+ outdeg_all

mean_in  = mean(Float64.(indeg_all))    # = L/S
sd_in    = std(Float64.(indeg_all))
mean_out = mean(Float64.(outdeg_all))   # = L/S
sd_out   = std(Float64.(outdeg_all))
mean_tot = mean(Float64.(totaldeg_all)) # = 2L/S
sd_tot   = std(Float64.(totaldeg_all))

# -----------------------------
# Print summary
# -----------------------------
println("=== California rocky intertidal metaweb summary ===")
println("Nodes (S):            ", S)
println("Links (L):            ", L)
println("Possible links:       ", Lpossible)
println("Connectance:          ", Connectance)
println()
println("Mean in-degree (prey/node):  ", mean_in)
println("SD in-degree:                ", sd_in)
println("Mean out-degree:             ", mean_out)
println("SD out-degree:               ", sd_out)
println("Mean total degree (in+out):  ", mean_tot)
println("SD total degree:             ", sd_tot)
println()

# -----------------------------
# Ready-to-paste compilation row (tab-separated)
# Put mean degree = mean_in (prey per node) and SD degree = sd_in by default
# -----------------------------
Name  = "California rocky intertidal"
Ecosystem = "California rocky intertidal zone"
Taxa = "parasite-inclusive food web"
Type = "directed trophic network"
Levels = ""                 # unknown -> blank
Region = ""                 # unknown -> blank
Temporal_range = ""         # unknown -> blank
HighestTaxRes = ""          # unknown -> blank
Reference = "Zilz et al. preprint (2024) doi:10.1101/2024.10.14.618335"

println("Name\tEcosystem\tTaxa\tType\tLevels\tRegion\tTemporal range\tNodes\tLinks\tConnectance\tMean degree\tSD degree\tHighest Taxonomic Resolution\tReference")
println(join([
    Name, Ecosystem, Taxa, Type, Levels, Region, Temporal_range,
    string(S), string(L),
    string(Connectance),
    string(mean_in),
    string(sd_in),
    HighestTaxRes,
    Reference
], "\t"))
