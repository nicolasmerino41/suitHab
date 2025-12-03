module Climate

export make_climate_grid

using Random, FFTW
using ..SpatialUtils: norm01

function make_climate_grid(nx, ny; kind=:gradient, seed=1)
    rng = MersenneTwister(seed)

    if kind == :gradient
        x = range(0,1; length=nx)
        G = [x[i] for i=1:nx, j=1:ny]
        return norm01(G)
    elseif kind == :ridge
        x = range(0,1; length=nx)
        y = range(0,1; length=ny)
        σ1, σ2 = 0.15, 0.2
        G = [exp(-((x[i]-0.7)^2+(y[j]-0.3)^2)/(2σ1^2)) +
             exp(-((x[i]-0.3)^2+(y[j]-0.7)^2)/(2σ2^2))
             for i in 1:nx, j in 1:ny]
        return norm01(G)
    elseif kind == :fractal
        noise = randn(rng, nx, ny)
        F = fft(noise)
        for i in 1:nx, j in 1:ny
            fx = min(i-1, nx-(i-1))
            fy = min(j-1, ny-(j-1))
            k = sqrt(fx^2 + fy^2)
            F[i,j] /= (k == 0 ? 1 : k)
        end
        C = real(ifft(F))
        return norm01(C)
    end
end

end
