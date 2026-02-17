using Statistics

# -----------------------------
# SETTINGS
# -----------------------------
path = "Metawebs/Sierras.csv"

allow_self = false      # only meaningful if square
threshold  = 0.0        # link exists if value > threshold

# Metadata placeholders (edit manually)
Name  = "Sierras"
Ecosystem = ""
Taxa = ""
TypeoOf = "bipartite network (incidence matrix)"
Levels = ""
Region = ""
Temporal_range = ""
HighestTaxRes = ""
Reference = ""

stripstr(x) = strip(String(x))

# -----------------------------
# READ + SPLIT (semicolon CSV from Excel; keeps empty A1)
# -----------------------------
lines = readlines(path)
rows = [split(chomp(l), ';'; keepempty=true) for l in lines]

# Make rectangular
maxlen = maximum(length.(rows))
for r in rows
    if length(r) < maxlen
        append!(r, fill("", maxlen - length(r)))
    end
end

nrowss = length(rows)
ncolss = maxlen
@assert nrowss >= 2 && ncolss >= 2 "File too small."

col_labels = stripstr.(rows[1][2:end])
row_labels = stripstr.([rows[i][1] for i in 2:nrowss])

R = length(row_labels)
C = length(col_labels)

# Numeric block A is R x C
A = Array{Float64}(undef, R, C)
for i in 1:R
    for j in 1:C
        s = stripstr(rows[i+1][j+1])
        A[i, j] = isempty(s) ? 0.0 : parse(Float64, s)
    end
end

# Boolean links
B = A .> threshold

# Optional diagonal removal only if square
if !allow_self && (R == C)
    for i in 1:R
        B[i,i] = false
    end
end

# -----------------------------
# LINKS + CONNECTANCE (bipartite)
# -----------------------------
L = count(B)
Lpossible = R * C
Connectance = Lpossible == 0 ? NaN : L / Lpossible

# -----------------------------
# Degrees
# Row degrees = number of partners in columns (row sums)
# Col degrees = number of partners in rows (col sums)
# -----------------------------
deg_rows = vec(sum(B, dims=2))  # length R
deg_cols = vec(sum(B, dims=1))  # length C

degrees_all = vcat(Float64.(deg_rows), Float64.(deg_cols))
Mean_degree = mean(degrees_all)
SD_degree   = std(degrees_all)

println("=== Sierras matrix summary ===")
println("Matrix type: bipartite (rectangular incidence)")
println("Row nodes (R): ", R)
println("Col nodes (C): ", C)
println("Total nodes (R+C): ", R + C)
println("Links (L): ", L)
println("Connectance L/(R*C): ", Connectance)
println("Mean degree (all nodes): ", Mean_degree)
println("SD degree (all nodes): ", SD_degree)
println()

println("=== Ready-to-paste TSV row ===")
println("Name\tEcosystem\tTaxa\tType\tLevels\tRegion\tTemporal range\tNodes\tLinks\tConnectance\tMean degree\tSD degree\tHighest Taxonomic Resolution\tReference")
println(join([
    Name, Ecosystem, Taxa, TypeOf, Levels, Region, Temporal_range,
    string(R + C), string(L),
    string(Connectance),
    string(Mean_degree),
    string(SD_degree),
    HighestTaxRes,
    Reference
], "\t"))
