module Grids

using Random, Statistics, LinearAlgebra

export Grid, make_grid_gradient, make_grid_patchy, make_grid_mosaic, make_grid_ridge,
       climate_tail_index

struct Grid
    nx::Int; ny::Int; C::Int
    xy::Matrix{Float64}      # 2×C
    climate::Vector{Float64} # C in [0,1]
end

# helper: coordinates as (x,y) in [0,1]×[0,1]
function _base_xy(nx::Int, ny::Int)
    xs = range(0, 1; length=nx)
    ys = range(0, 1; length=ny)
    C  = nx*ny
    xy = Matrix{Float64}(undef, 2, C)
    k=1
    for j in 1:ny, i in 1:nx
        xy[:,k] .= (xs[i], ys[j])
        k+=1
    end
    xy
end

# normalize to [0,1]
_norm01!(v) = (v .-= minimum(v); v ./= (maximum(v)+eps()); v)

# 1) smooth gradient + small periodic wiggles
function make_grid_gradient(nx::Int, ny::Int; seed::Int=42, texture::Float64=0.10)
    rng = MersenneTwister(seed)
    xy = _base_xy(nx,ny) 
    C=size(xy,2)
    clim = similar(view(xy,1,:))
    @inbounds for k in 1:C
        x,y = xy[1,k], xy[2,k]
        base   = 0.7x + 0.3y
        wiggle = texture * sin(6π*x) * sin(6π*y)
        clim[k] = base + wiggle
    end
    _norm01!(clim)
    Grid(nx,ny,nx*ny,xy,clim)
end

# 2) patchy noise (fractal-ish)
function make_grid_patchy(nx::Int, ny::Int; seed::Int=43, octaves::Int=3, amp::Float64=0.5)
    rng = MersenneTwister(seed)
    xy  = _base_xy(nx,ny); C=size(xy,2)
    clim = zeros(Float64, C)
    for o in 0:octaves-1
        f = 2.0^o
        @inbounds for k in 1:C
            x,y = xy[1,k], xy[2,k]
            clim[k] += (amp/(f)) * (sin(2π*f*x) + cos(2π*f*y))
        end
    end
    _norm01!(clim)
    Grid(nx,ny,nx*ny,xy,clim)
end

# 3) mosaic steps (checkerboardy)
function make_grid_mosaic(nx::Int, ny::Int; seed::Int=44, cells::Int=6)
    rng = MersenneTwister(seed)
    xy  = _base_xy(nx,ny); C=size(xy,2)
    clim = similar(view(xy,1,:))
    @inbounds for k in 1:C
        x,y = xy[1,k], xy[2,k]
        ix = min(Int(floor(x*cells))+1, cells)
        iy = min(Int(floor(y*cells))+1, cells)
        v  = sin(π*ix/cells) + cos(π*iy/cells)
        clim[k] = v
    end
    _norm01!(clim)
    Grid(nx,ny,nx*ny,xy,clim)
end

# 4) ridge/valley (corridor-like anisotropy)
function make_grid_ridge(nx::Int, ny::Int; seed::Int=45)
    xy  = _base_xy(nx,ny); C=size(xy,2)
    clim = similar(view(xy,1,:))
    @inbounds for k in 1:C
        x,y = xy[1,k], xy[2,k]
        clim[k] = 0.6*sin(4π*x) + 0.4*cos(2π*y) + 0.2*sin(6π*(x+y))
    end
    _norm01!(clim)
    Grid(nx,ny,nx*ny,xy,clim)
end

"Variance explained by best linear combo of x,y (R² of OLS)."
function climate_tail_index(grid::Grid)
    X = hcat(ones(grid.C), vec(grid.xy[1,:]), vec(grid.xy[2,:]))
    β = X \ grid.climate
    ŷ = X*β
    ssr = sum((ŷ .- mean(grid.climate)).^2)
    sst = sum((grid.climate .- mean(grid.climate)).^2) + eps()
    ssr/sst  # ∈ [0,1]; larger ⇒ stronger directional 'tail'
end

end # module
