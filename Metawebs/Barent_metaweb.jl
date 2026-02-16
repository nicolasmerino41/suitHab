using CSV, DataFrames

path = "Metawebs/Barent_metaweb.txt"  # <- change to your filename

# Read as TSV (tab-separated). `ignorerepeated=true` helps if there are occasional extra tabs.
df = CSV.read(path, DataFrame;
              delim = '\t',
              ignorerepeated = true,
              missingstring = "",
              quotechar = '"')  # quotechar mostly irrelevant here, but harmless

# Expect columns: PREY, PREDATOR, CODE, REFERENCE
# Make it resilient to weird capitalization by normalizing names once:
rename!(df, Symbol.(strip.(String.(names(df)))))

# Grab endpoints
prey = strip.(String.(df.PREY))
pred = strip.(String.(df.PREDATOR))

# --- options ---
directed   = true
allow_self = false

# Optional: filter by CODE if you decide only some interaction codes count as links
# Example: keep only CODE == 1
# df = df[df.CODE .== 1, :]
# prey = strip.(String.(df.PREY)); pred = strip.(String.(df.PREDATOR))

# Species list inferred from endpoints
species = unique(vcat(prey, pred))
S = length(species)

# Realized links: unique pairs
seen = Set{Tuple{String,String}}()
for (u, v) in zip(prey, pred)
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
