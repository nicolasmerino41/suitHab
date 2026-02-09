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
