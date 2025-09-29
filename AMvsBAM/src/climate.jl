# ----------------------------------------------------------
# NEW: climate helpers (gradient, ridge, fractal) + resample
# ----------------------------------------------------------
module Climate
using Random, Statistics, Distributions
export make_climate_grid, estimate_corr_length, resample_to

"Normalize array to [0,1]"
_norm01(A) = (A .- minimum(A)) ./ max(eps(), (maximum(A)-minimum(A)))

"Simple bilinear resampler to (nx,ny)"
function resample_to(C::AbstractMatrix{<:Real}, nx::Int, ny::Int)
    sx, sy = size(C)
    xscale = (sx-1) / max(1, nx-1)
    yscale = (sy-1) / max(1, ny-1)
    out = Matrix{Float64}(undef, nx, ny)
    for i in 1:nx, j in 1:ny
        x = (i-1) * xscale + 1
        y = (j-1) * yscale + 1
        x0 = clamp(floor(Int, x), 1, sx-1)
        y0 = clamp(floor(Int, y), 1, sy-1)
        dx = x - x0
        dy = y - y0
        a = C[x0, y0]; b = C[x0+1, y0]
        c = C[x0, y0+1]; d = C[x0+1, y0+1]
        out[i,j] = (1-dx)*(1-dy)*a + dx*(1-dy)*b + (1-dx)*dy*c + dx*dy*d
    end
    _norm01(out)
end

"Create climate grid: :gradient | :ridge | :fractal (β≈1 ~ 1/f)"
function make_climate_grid(nx::Int, ny::Int; kind::Symbol=:gradient, seed::Int=1, β::Real=1.0)
    rng = MersenneTwister(seed)
    if kind == :gradient
        # linear ramp with a touch of smooth noise
        x = range(0, 1; length=nx); y = range(0, 1; length=ny)
        G = [x[i] for i=1:nx, j=1:ny]
        # add gentle 2D cosine noise
        N = [0.05*sin(2π*3*x[i])*cos(2π*2*y[j]) for i=1:nx, j=1:ny]
        return _norm01(G .+ N)
    elseif kind == :ridge
        # two broader ridges + low background (laxer than before)
        x = range(0, 1; length=nx); y = range(0, 1; length=ny)
        σ1, σ2 = 0.10, 0.20       # <- widened from 0.05, 0.07
        base  = 0.15              # <- small background everywhere
        G = [ base +
              exp(-((x[i]-0.70)^2 + (y[j]-0.22)^2)/(2*σ1^2)) +
              exp(-((x[i]-0.30)^2 + (y[j]-0.70)^2)/(2*σ2^2))
              for i=1:nx, j=1:ny ]
        return _norm01(G)
    elseif kind == :fractal
        # spectral synthesis with power ~ k^{-β}
        # build on a somewhat larger field then crop/normalize
        kx = collect(0:(nx-1)); ky = collect(0:(ny-1))
        K = Matrix{Float64}(undef, nx, ny)
        for i in 1:nx, j in 1:ny
            fx = min(kx[i], nx-kx[i])
            fy = min(ky[j], ny-ky[j])
            k = sqrt(fx^2 + fy^2)
            K[i,j] = k == 0 ? 0.0 : 1.0 / (k^β)
        end
        # random phases
        phase = 2π .* rand(rng, nx, ny)
        re = K .* cos.(phase); im = K .* sin.(phase)
        F = Complex.(re, im)
        # inverse FFT via real ifft (naive DFT fallback-free)
        # use FFTW if available; otherwise simple spatial filter:
        # convolve white noise with Gaussian of varying scales
        # ---- fallback without FFTW ----
        noise = randn(rng, nx, ny)
        # multi-scale blur
        function blur(A, σ)
            k = max(3, Int(round(σ*6)))
            k % 2 == 0 && (k += 1)
            # 1D kernel
            xs = collect(-div(k,2):div(k,2))
            g = exp.(-(xs.^2) ./ (2σ^2)); g ./= sum(g)
            # separable blur
            tmp = similar(A)
            # x pass
            for j in 1:size(A,2)
                for i in 1:size(A,1)
                    s = 0.0
                    for (t, w) in zip(xs, g)
                        ii = clamp(i+t, 1, size(A,1))
                        s += w * A[ii, j]
                    end
                    tmp[i,j] = s
                end
            end
            out = similar(A)
            # y pass
            for i in 1:size(A,1)
                for j in 1:size(A,2)
                    s = 0.0
                    for (t, w) in zip(xs, g)
                        jj = clamp(j+t, 1, size(A,2))
                        s += w * tmp[i, jj]
                    end
                    out[i,j] = s
                end
            end
            out
        end
        F1 = blur(noise, 1.0)
        F2 = blur(noise, 2.0)
        F3 = blur(noise, 4.0)
        C = (F1 .+ 0.5F2 .+ 0.25F3)
        return _norm01(C)
    else
        error("Unknown climate kind: $kind")
    end
end

"Estimate correlation length L (in cells) along x using first 1/e drop of the ACF."
function estimate_corr_length(C::AbstractMatrix{<:Real})
    nx, ny = size(C)
    # take row means to get a 1D profile along x
    prof = vec(mean(C; dims=2))
    prof .-= mean(prof)
    denom = sum(prof.^2)
    denom ≈ 0 && return max(1, nx ÷ 10)
    acf = zeros(Float64, nx)
    for lag in 0:nx-1
        num = 0.0
        for i in 1:(nx-lag)
            num += prof[i] * prof[i+lag]
        end
        acf[lag+1] = num / denom
    end
    # find first lag where acf <= exp(-1)
    target = exp(-1)
    for lag in 1:nx
        if acf[lag] <= target
            return lag
        end
    end
    return nx ÷ 3
end

end # module climate