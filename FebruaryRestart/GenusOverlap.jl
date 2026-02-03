# Build binomial names from GlobalTherm
gt_species = Set{String}()

for i in 1:nrow(gt)
    g = gt.Genus[i]
    s = gt.Species[i]

    if g !== missing && s !== missing
        push!(gt_species, strip(string(g, " ", s)))
    end
end

@info "GlobalTherm species" length(gt_species)

metaweb_species = Set(unique_species)

@info "Metaweb species" length(metaweb_species)

shared_species = intersect(gt_species, metaweb_species)
only_gt        = setdiff(gt_species, metaweb_species)
only_metaweb   = setdiff(metaweb_species, gt_species)

@info "Species overlap summary" (
    gt_total        = length(gt_species),
    metaweb_total   = length(metaweb_species),
    shared          = length(shared_species),
    frac_gt_shared  = length(shared_species) / length(gt_species),
    frac_web_shared = length(shared_species) / length(metaweb_species)
)

# GlobalTherm genera
gt_genera = Set{String}()

for g in gt.Genus
    g === missing && continue
    push!(gt_genera, strip(g))
end

@info "GlobalTherm genera" length(gt_genera)

# Metaweb genera (from species names)
metaweb_genera = Set{String}()

for sp in unique_species
    isempty(sp) && continue
    genus = first(split(sp, " "))
    push!(metaweb_genera, strip(genus))
end

@info "Metaweb genera" length(metaweb_genera)

shared_genera = intersect(gt_genera, metaweb_genera)
only_gt_gen   = setdiff(gt_genera, metaweb_genera)
only_web_gen  = setdiff(metaweb_genera, gt_genera)

@info "Genus overlap summary" (
    gt_total        = length(gt_genera),
    metaweb_total   = length(metaweb_genera),
    shared          = length(shared_genera),
    frac_gt_shared  = length(shared_genera) / length(gt_genera),
    frac_web_shared = length(shared_genera) / length(metaweb_genera)
)
