using CSV, DataFrames

path = "Metawebs/KF_metaweb.csv"

edges = CSV.read(path, DataFrame)

# If your CSV headers are exactly: resource, consumer, Reference
# (Reference is ignored for connectance)
from = strip.(String.(edges.resource))
to   = strip.(String.(edges.consumer))

# Species list (inferred)
species = unique(vcat(from, to))
S = length(species)

# Options
allow_self = false  # set true if you want to allow i→i links
directed   = true   # trophic metawebs are usually directed

# Realized links (unique)
seen = Set{Tuple{String,String}}()
for (u, v) in zip(from, to)
    if !allow_self && u == v
        continue
    end
    if directed
        push!(seen, (u, v))
    else
        push!(seen, u <= v ? (u, v) : (v, u))
    end
end
L = length(seen)

# Possible links
Lpossible = if directed
    allow_self ? S^2 : S*(S-1)
else
    allow_self ? (S*(S+1)) ÷ 2 : (S*(S-1)) ÷ 2
end

C = Lpossible == 0 ? NaN : L / Lpossible

println("S = $S")
println("L = $L")
println("Lpossible = $Lpossible")
println("Connectance C = $C")
