mw_tf = CSV.read(joinpath(@__DIR__, "thermofresh_globi_metaweb_fish_predators.csv"), DataFrame; missingstring="", ntasks=1, pool=false) # CSV.read("February_restart/thermofresh_globi_metaweb_fish_predators.csv", DataFrame; missingstring="", ntasks=1, pool=false)
mw_old = CSV.read(joinpath(@__DIR__, "../outputs_imputed_globi_edges.csv"), DataFrame; missingstring="", ntasks=1, pool=false) # CSV.read("outputs_imputed_globi_edges.csv", DataFrame; missingstring="", ntasks=1, pool=false)

# Canonize into the same column names
mw_tf2 = unique(DataFrame(
    pred = canon_taxon.(string.(mw_tf.predator)),
    prey = canon_taxon.(string.(mw_tf.prey))
))

mw_old2 = unique(DataFrame(
    pred = canon_taxon.(string.(mw_old.pred)),
    prey = canon_taxon.(string.(mw_old.prey))
))

# Intersection (pairs present in both)
overlap = innerjoin(mw_tf2, mw_old2, on=[:pred, :prey])

@info "TF unique edges" nrow(mw_tf2)
@info "Olden unique edges" nrow(mw_old2)
@info "Overlap edges" nrow(overlap)

first(overlap, 50) |> display
CSV.write("overlap_tf_olden_edges.csv", overlap)

# ---- Unique species (nodes) in each metaweb ----
tf_species = Set(vcat(mw_tf2.pred, mw_tf2.prey))
old_species = Set(vcat(mw_old2.pred, mw_old2.prey))

# Species in ThermoFresh×GloBI but NOT in Olden×GloBI
tf_not_in_old = setdiff(tf_species, old_species)

# Species in Olden×GloBI but NOT in ThermoFresh×GloBI
old_not_in_tf = setdiff(old_species, tf_species)

@info "TF unique species (pred ∪ prey)" length(tf_species)
@info "Olden unique species (pred ∪ prey)" length(old_species)
@info "Species in TF but not in Olden" length(tf_not_in_old)
@info "Species in Olden but not in TF" length(old_not_in_tf)

# Optional: inspect / save list
first(sort!(collect(tf_not_in_old)), 50) |> display
CSV.write("tf_species_not_in_olden.csv", DataFrame(species = sort!(collect(tf_not_in_old))))
