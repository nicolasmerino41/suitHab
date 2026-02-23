using CSV, DataFrames
using HTTP, JSON3
using URIs
using Pkg

# -----------------------------
# CONFIG
# -----------------------------
OUTROOT = joinpath(@__DIR__, "RPlots/Plots/outputs_merged_all_ctmax")
edge_path = joinpath(OUTROOT, "edge_table.csv")

worms_base = "https://www.marinespecies.org/rest"

# -----------------------------
# HELPERS
# -----------------------------

# WoRMS query
function worms_records_by_name(name::AbstractString)
    enc = URIs.escapeuri(name)
    url = "$worms_base/AphiaRecordsByName/$enc?like=false&marine_only=false&offset=1"

    try
        r = HTTP.get(url; readtimeout=30)
        r.status != 200 && return nothing
        return JSON3.read(String(r.body))
    catch
        return nothing
    end
end

# Pick accepted record if present
function pick_best_record(recs)
    recs === nothing && return nothing
    isempty(recs) && return nothing

    for rec in recs
        if hasproperty(rec, :status) &&
           lowercase(String(rec.status)) == "accepted"
            return rec
        end
    end

    return recs[1]
end

# Normalize WoRMS flags safely
function habitat_from_worms_record(rec)
    rec === nothing && return "unknown"

    function normflag(x)
        x === nothing && return 0
        x === true && return 1
        x === false && return 0
        return try
            Int(x)
        catch
            0
        end
    end

    isM = hasproperty(rec, :isMarine)      ? normflag(rec.isMarine)      : 0
    isF = hasproperty(rec, :isFreshwater)  ? normflag(rec.isFreshwater)  : 0
    isT = hasproperty(rec, :isTerrestrial) ? normflag(rec.isTerrestrial) : 0

    s = String[]
    isM == 1 && push!(s, "marine")
    isF == 1 && push!(s, "freshwater")
    isT == 1 && push!(s, "terrestrial")

    isempty(s) && return "unknown"
    length(s) == 1 && return s[1]
    return "mixed"
end

# Edge classification rule
function edge_habitat(hpred::AbstractString, hprey::AbstractString)
    (hpred == "unknown" || hprey == "unknown") && return "unknown"
    (hpred == hprey && hpred in ("marine","freshwater","terrestrial")) && return hpred
    return "mixed"
end

# -----------------------------
# MAIN
# -----------------------------
edges = CSV.read(edge_path, DataFrame)

taxa = unique(vcat(String.(edges.pred), String.(edges.prey)))

hab = Dict{String,String}()

# --- WoRMS classification ---
for (i, sp) in enumerate(taxa)
    recs = worms_records_by_name(sp)
    best = pick_best_record(recs)
    hab[sp] = habitat_from_worms_record(best)

    sleep(0.05)

    if i % 200 == 0
        @info "Classified $i / $(length(taxa)) taxa"
    end
end

# -----------------------------
# MANUAL HABITAT FIXES
# -----------------------------

manual_habitat = Dict(

    # Terrestrial vertebrates
    "Erinaceus concolor" => "terrestrial",
    "Podarcis lilfordi" => "terrestrial",
    "Crocidura russula" => "terrestrial",
    "Dipodillus dasyurus" => "terrestrial",
    "Mus macedonicus" => "terrestrial",
    "Acanthodactylus erythrurus" => "terrestrial",
    "Lacerta schreiberi" => "terrestrial",
    "Podarcis bocagei" => "terrestrial",
    "Podarcis muralis" => "terrestrial",
    "Podarcis tiliguerta" => "terrestrial",
    "Psammodromus hispanicus" => "terrestrial",
    "Tarentola boettgeri" => "terrestrial",

    # Freshwater insects & invertebrates
    "Serratella ignita" => "freshwater",
    "Elmis aenea" => "freshwater",
    "Leuctra fusca" => "freshwater",
    "Rhyacophila dorsalis" => "freshwater",
    "Limnius volckmari" => "freshwater",
    "Nemurella pictetii" => "freshwater",
    "Rhithrogena semicolorata" => "freshwater",
    "Hydropsyche pellucidula" => "freshwater",
    "Agabus brunneus" => "freshwater",
    "Deronectes latus" => "freshwater",
    "Sialis lutaria" => "freshwater",
    "Cordulegaster boltonii" => "freshwater",
    "Agabus didymus" => "freshwater",
    "Simulium vittatum" => "freshwater",
    "Heptagenia sulphurea" => "freshwater",
    "Ecdyonurus insignis" => "freshwater",
    "Ischnura pumilio" => "freshwater",
    "Coenagrion puella" => "freshwater",
    "Calopteryx virgo" => "freshwater",
    "Ischnura elegans" => "freshwater",
    "Caenis luctuosa" => "freshwater",
    "Agapetus fuscipes" => "freshwater",
    "Dinocras cephalotes" => "freshwater",
    "Aphelocheirus aestivalis" => "freshwater",
    "Agabus bipustulatus" => "freshwater",
    "Drusus annulatus" => "freshwater",
    "Epeorus longimanus" => "freshwater",
    "Baetis bicaudatus" => "freshwater",
    "Rhyacophila brunnea" => "freshwater",
    "Drunella coloradensis" => "freshwater",
    "Daphnia pulicaria" => "freshwater",
    "Pteronarcys californica" => "freshwater",

    # Freshwater fish
    "Etroplus suratensis" => "freshwater"
)

# Apply overrides
for (sp, hlabel) in manual_habitat
    hab[sp] = hlabel
end

# -----------------------------
# EXPORT TAXON HABITATS
# -----------------------------
hab_df = DataFrame(taxon = taxa, habitat = [hab[t] for t in taxa])
CSV.write(joinpath(OUTROOT, "habitat_by_taxon.csv"), hab_df)

# -----------------------------
# EDGE CLASSIFICATION
# -----------------------------
edges_h = transform(edges,
    :pred => (p -> [hab[String(x)] for x in p]) => :pred_habitat,
    :prey => (p -> [hab[String(x)] for x in p]) => :prey_habitat
)

edges_h.edge_habitat =
    [edge_habitat(edges_h.pred_habitat[i],
                  edges_h.prey_habitat[i]) for i in 1:nrow(edges_h)]

# -----------------------------
# EDGE COUNTS
# -----------------------------
counts = combine(groupby(edges_h, :edge_habitat), nrow => :n_edges)
sort!(counts, :n_edges, rev=true)

CSV.write(joinpath(OUTROOT, "edge_habitat_counts.csv"), counts)

println(counts)

lookup = Dict(row.edge_habitat => row.n_edges for row in eachrow(counts))

println("\nMarine:      ", get(lookup, "marine", 0))
println("Freshwater:  ", get(lookup, "freshwater", 0))
println("Terrestrial: ", get(lookup, "terrestrial", 0))
println("Mixed:       ", get(lookup, "mixed", 0))
println("Unknown:     ", get(lookup, "unknown", 0))