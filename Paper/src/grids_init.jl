# src/grids_init.jl
module GridsInit

using ..Grids

export build_default_grids, as_pairs, print_grid_diagnostics

"""
    build_default_grids(nx, ny; seeds=(11,12,13,14))

Builds the four climate grids used everywhere in the paper and
returns a NamedTuple:

(; grad, patch, mosaic, ridge)

Each entry is a `Grids.Grid`.
"""
function build_default_grids(nx::Int, ny::Int; seeds::NTuple{4,Int}=(11,12,13,14))
    sgrad, spatch, smosa, sridge = seeds
    grad   = Grids.make_grid_gradient(nx,ny; seed=sgrad,  texture=0.10)
    patch  = Grids.make_grid_patchy(nx,ny;  seed=spatch,  octaves=3, amp=0.5)
    mosaic = Grids.make_grid_mosaic(nx,ny;  seed=smosa,   cells=6)
    ridge  = Grids.make_grid_ridge(nx,ny;   seed=sridge)
    (; grad, patch, mosaic, ridge)
end

"Return a Vector of (name, grid) pairs for sweep functions."
as_pairs(G) = [("gradient", G.grad),
               ("patchy",   G.patch),
               ("mosaic",   G.mosaic),
               ("ridge",    G.ridge)]

"Print simple grid diagnostics (R² tail index)."
function print_grid_diagnostics(G)
    for (nm, g) in as_pairs(G)
        D = Grids.climate_tail_index(g)
        @info "grid $nm  tail-index R²=$(round(D,digits=3))"
    end
end

end # module
