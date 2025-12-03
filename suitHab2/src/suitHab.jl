module suitHab

# == Load submodules ==
include("Metaweb/Types.jl")
include("Metaweb/MetawebUtils.jl")
include("Metaweb/PPM.jl")
include("Spatial/SpatialUtils.jl")
include("Spatial/Climate.jl")
include("Spatial/Biotic.jl")
include("Spatial/Movement.jl")
include("Niches/Niches.jl")
include("Assembly/Suitability.jl")
include("Assembly/RealizedNetwork.jl")
include("Assembly/LossScenarios.jl")
include("Visualisation/Plotting.jl")

# == Bring functions into namespace ==

using .PPM: build_ppm
using .Climate: make_climate_grid
using .Biotic: build_basal_biotic, propagate_biotic
using .Movement: apply_movement
using .Niches: make_niches
using .Suitability: abiotic_suitability, biotic_suitability,
                    combine_suitability, movement_filter
using .RealizedNetwork: realized_metaweb
using .LossScenarios: apply_loss
using .SpatialUtils
using .Plotting

# == Export clean user API ==
export build_ppm,
       make_climate_grid,
       build_basal_biotic, propagate_biotic,
       apply_movement,
       make_niches,
       abiotic_suitability, combine_suitability, movement_filter,
       realized_metaweb,
       apply_loss,
       plot_grid

end
