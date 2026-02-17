using CSV, DataFrames, Dates, Statistics, Printf

src_path = "Metawebs/AtlanticForestHummingbirds.csv"

stripstr(x) = strip(String(x))

function parse_dates_safe(col)::Vector{Union{Missing,Date}}
    fmts = (
        dateformat"yyyy-mm-dd",
        dateformat"yyyy/mm/dd",
        dateformat"dd/mm/yyyy",
        dateformat"mm/dd/yyyy",
        dateformat"yyyy-mm",
        dateformat"yyyy/mm",
        dateformat"yyyy",
    )
    out = Vector{Union{Missing,Date}}(undef, length(col))
    for (i, v) in pairs(col)
        if v === missing
            out[i] = missing
            continue
        end
        s = strip(String(v))
        if isempty(s)
            out[i] = missing
            continue
        end
        d = nothing
        for f in fmts
            try
                d = Date(s, f)
                break
            catch
            end
        end
        out[i] = d === nothing ? missing : d
    end
    return out
end

# -----------------------------
# Load data
# -----------------------------
df = CSV.read(src_path, DataFrame)

# Endpoints
plants = stripstr.(df.plant_revised_name)
humbs  = stripstr.(df.hummingbird_revised)

# Drop empty endpoints
keep = .!(isempty.(plants) .| isempty.(humbs))
plants = plants[keep]
humbs  = humbs[keep]

# Unique realized links (plant-hummingbird pairs)
links_df = DataFrame(plant = plants, hummingbird = humbs)
unique!(links_df)

P = length(unique(links_df.plant))
H = length(unique(links_df.hummingbird))
Nodes = P + H
Links = nrow(links_df)

# Bipartite connectance (potential links = P*H)
Lpossible = P * H
Connectance = Lpossible == 0 ? NaN : Links / Lpossible

# Degree stats (partners per node; undirected bipartite degrees across all nodes)
deg_plants = combine(groupby(links_df, :plant),
    :hummingbird => (x -> length(unique(x))) => :deg).deg
deg_humbs = combine(groupby(links_df, :hummingbird),
    :plant => (x -> length(unique(x))) => :deg).deg

degrees_all = vcat(Float64.(deg_plants), Float64.(deg_humbs))
Mean_degree = mean(degrees_all)
SD_degree   = std(degrees_all)

# Temporal range best-effort
Temporal_range = ""
if :date in propertynames(df)
    dparsed = parse_dates_safe(df.date)
    dclean = collect(skipmissing(dparsed))
    if !isempty(dclean)
        Temporal_range = string(minimum(dclean), " to ", maximum(dclean))
    end
end

# -----------------------------
# Print results
# -----------------------------
println("=== AtlanticForestsHummingbirds summary ===")
println("Unique plants (P):        ", P)
println("Unique hummingbirds (H):  ", H)
println("Nodes (P+H):              ", Nodes)
println("Links (unique pairs):     ", Links)
println("Connectance (L/(P*H)):    ", Connectance)
println("Mean degree (all nodes):  ", Mean_degree)
println("SD degree (all nodes):    ", SD_degree)
println("Temporal range:           ", isempty(Temporal_range) ? "(blank)" : Temporal_range)
println()

# -----------------------------
# Ready-to-copy compilation row (tab-separated)
# Use blanks for unknown fields
# -----------------------------
Name  = "AtlanticForestsHummingbirds"
Ecosystem = "Atlantic Forest"
Taxa = "plants-hummingbirds"  # plain hyphen avoids Excel encoding headaches
TypeOf = "bipartite interaction network"
Levels = ""
Region = ""
HighestTaxRes = ""
Reference = ""

# Print header (so you can paste into Excel safely)
println("Name\tEcosystem\tTaxa\tType\tLevels\tRegion\tTemporal range\tNodes\tLinks\tConnectance\tMean degree\tSD degree\tHighest Taxonomic Resolution\tReference")
println(join([
    Name, Ecosystem, Taxa, TypeOf, Levels, Region,
    Temporal_range,
    string(Nodes), string(Links),
    @sprintf("%.12g", Connectance),
    @sprintf("%.12g", Mean_degree),
    @sprintf("%.12g", SD_degree),
    HighestTaxRes,
    Reference
], "\t"))
