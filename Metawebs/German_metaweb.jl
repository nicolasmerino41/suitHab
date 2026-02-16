using CSV, DataFrames

path = "Metawebs/German_metaweb.csv"  # <- make sure you're in the right folder

df = CSV.read(path, DataFrame)

# Basic cleaning
pred = strip.(String.(df.PREDATOR))
prey = strip.(String.(df.PREY))

# --- options ---
directed   = true
allow_self = false

# Optional filters (uncomment and edit if you need them)
# 1) Filter by area
# df = df[df.area .== "German Bight", :]
# pred = strip.(String.(df.PREDATOR)); prey = strip.(String.(df.PREY))

# 2) Filter by lifestage
# df = df[(df.pred_lifestage .== "adult") .& (df.prey_lifestage .== "adult"), :]
# pred = strip.(String.(df.PREDATOR)); prey = strip.(String.(df.PREY))

# 3) Filter by Quality (depends on encoding; examples: "high", "A", 1..5, etc.)
# df = df[df.Quality .== "high", :]
# pred = strip.(String.(df.PREDATOR)); prey = strip.(String.(df.PREY))

# Species list inferred from endpoints
species = unique(vcat(pred, prey))
S = length(species)

# Realized links: unique prey->predator pairs (trophic direction)
seen = Set{Tuple{String,String}}()
for (p, r) in zip(prey, pred)  # (resource, consumer)
    if !allow_self && p == r
        continue
    end
    if directed
        push!(seen, (p, r))
    else
        push!(seen, p <= r ? (p, r) : (r, p))
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
