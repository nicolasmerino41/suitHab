# ================================
# Knee: trophic height vs basal richness
# ================================
# 1) Trophic height per cell (longest predator→prey chain among PRESENT spp.)
function trophic_height_per_cell(P::BitMatrix, pool::SpeciesPool)
    S, C = size(P)
    heights = zeros(Int, C)
    ord = sortperm(pool.masses)  # increasing mass (basal first)

    for c in 1:C
        present = findall(@view P[:, c])
        if isempty(present)
            heights[c] = 0
            continue
        end
        h = zeros(Int, S)
        # compute prey first, then predators
        for s in ord
            P[s,c] || continue
            prey = pool.E[s]
            if isempty(prey)
                h[s] = 1                        # basal
            else
                maxprey = 0
                @inbounds for q in prey
                    if P[q,c]
                        maxprey = max(maxprey, h[q])
                    end
                end
                h[s] = max(1, maxprey + 1)      # if no present prey, stays 1
            end
        end
        heights[c] = maximum(h[present])
    end
    return heights
end

# 2) Basal richness per cell
basal_richness_vec(P::BitMatrix, basal::BitVector) = vec(sum(@view P[basal, :]; dims=1))

# 3) Bin x vs y and take means
function bin_mean(x::AbstractVector, y::AbstractVector; nbins::Int=15, min_per_bin::Int=5)
    @assert length(x) == length(y)
    xmin, xmax = extrema(x)
    edges = range(xmin, xmax; length=nbins+1)
    mids  = (edges[1:end-1] .+ edges[2:end]) ./ 2
    μ = similar(mids, Float64); fill!(μ, NaN)
    for i in 1:nbins
        idx = findall((x .>= edges[i]) .& (x .< edges[i+1]))
        if length(idx) >= min_per_bin
            μ[i] = mean(@view y[idx])
        end
    end
    return mids, μ
end

# 4) Two-segment piecewise linear fit to estimate knee (breakpoint)
#    Returns x* (break location) and fitted y-hat (for plotting, optional)
function knee_piecewise(x::Vector{Float64}, y::Vector{Float64}; min_seg::Int=3)
    # keep only finite points and sort by x
    keep = findall(isfinite.(x) .& isfinite.(y))
    x = x[keep]; y = y[keep]
    p = sortperm(x); x = x[p]; y = y[p]
    n = length(x)
    if n < 2*min_seg
        return (NaN, fill(NaN, n))  # not enough points
    end
    # small linear least-squares helper
    linfit_sse(xv, yv) = begin
        X = hcat(ones(length(xv)), xv)
        β = X \ yv
        ŷ = X * β
        sse = sum((yv .- ŷ).^2)
        (β, ŷ, sse)
    end
    best_sse = Inf
    best_i   = 0
    yhat_best = similar(y)
    for i in min_seg:(n-min_seg)
        βL, ŷL, sseL = linfit_sse(x[1:i],     y[1:i])
        βR, ŷR, sseR = linfit_sse(x[i+1:end], y[i+1:end])
        sse = sseL + sseR
        if sse < best_sse
            best_sse = sse
            best_i   = i
            yhat_best = vcat(ŷL, ŷR)
        end
    end
    xstar = x[best_i]
    return (xstar, yhat_best)
end

# 5) Compute knee curve data for a given keep mask
function knee_curve_for_mask(grid::Grid, pool::SpeciesPool; τ=0.5, keepmask, nbins=15)
    Z = climate_pass(pool, grid; τ=τ)
    Z = apply_mask(Z, keepmask)
    P = assemble(Z, pool)
    br = Float64.(basal_richness_vec(P, pool.basal))   # x
    th = Float64.(trophic_height_per_cell(P, pool))    # y
    mids, μ = bin_mean(br, th; nbins=nbins, min_per_bin=10)
    xstar, _ = knee_piecewise(collect(mids), collect(μ); min_seg=3)
    return (mids, μ, xstar)
end

# 6) Plot random vs clustered knee at a chosen loss fraction
function knee_plot(; grid::Grid, pool::SpeciesPool, τ::Float64=0.5, keep_frac::Float64=0.6, nbins::Int=15,
                    seed_r::Int=101, seed_c::Int=202, nseeds_cluster::Int=6)
    km_r = random_mask(grid.C, keep_frac; seed=seed_r)
    km_c = clustered_mask(grid, keep_frac; nseeds=nseeds_cluster, seed=seed_c)

    xR, yR, kR = knee_curve_for_mask(grid, pool; τ=τ, keepmask=km_r, nbins=nbins)
    xC, yC, kC = knee_curve_for_mask(grid, pool; τ=τ, keepmask=km_c, nbins=nbins)

    begin
        fig = Figure(resolution=(820,380))
        ax  = Axis(fig[1,1], xlabel="Basal richness per cell", ylabel="Mean trophic height",
                   title = "Knee under equal-area loss (keep = $(keep_frac))")
        lines!(ax, xR, yR, label="Random")
        lines!(ax, xC, yC, label="Clustered")
        if isfinite(kR); vlines!(ax, [kR], color=:black, linestyle=:dash); end
        if isfinite(kC); vlines!(ax, [kC], color=:red,   linestyle=:dash);  end
        axislegend(ax, position=:lt)
        display(fig)
    end
    return (; random_knee=kR, clustered_knee=kC)
end

# # assuming you already created `grid` and `pool` as before
knees = knee_plot(; grid=grid, pool=pool, τ=0.5, keep_frac=0.6, nbins=15)
println("Estimated knee (random vs clustered): ", knees)

# Probability-of-persistence knee (TH ≥ T)
function knee_prob_for_mask(grid::Grid, pool::SpeciesPool; τ=0.5, keepmask, T::Int=3, nbins::Int=15, min_bin::Int=20)
    Z = climate_pass(pool, grid; τ=τ) |> Z -> apply_mask(Z, keepmask)
    P = assemble(Z, pool)
    br = Float64.(basal_richness_vec(P, pool.basal))
    th = Float64.(trophic_height_per_cell(P, pool))
    y  = th .>= T

    # bin
    xmin, xmax = extrema(br)
    edges = range(xmin, xmax; length=nbins+1)
    mids  = (edges[1:end-1] .+ edges[2:end]) ./ 2
    p̂     = fill(NaN, nbins)
    for i in 1:nbins
        idx = findall((br .>= edges[i]) .& (br .< edges[i+1]))
        if length(idx) ≥ min_bin
            p̂[i] = mean(y[idx])  # probability TH≥T in this bin
        end
    end

    # knee = first mid where p̂ ≥ 0.5
    knee = NaN
    for i in 1:nbins
        if isfinite(p̂[i]) && p̂[i] ≥ 0.5
            knee = mids[i]; break
        end
    end
    return (mids, p̂, knee)
end

# Plot (random vs clustered) at the same keep fraction
function knee_prob_plot(; grid::Grid, pool::SpeciesPool, τ=0.5, keep_frac=0.6, nbins=15, T::Int=3,
                         seed_r=101, seed_c=202, nseeds_cluster=2)
    km_r = random_mask(grid.C, keep_frac; seed=seed_r)
    km_c = clustered_mask(grid, keep_frac; nseeds=nseeds_cluster, seed=seed_c)
    xR, pR, kR = knee_prob_for_mask(grid, pool; τ=τ, keepmask=km_r, T=T, nbins=nbins)
    xC, pC, kC = knee_prob_for_mask(grid, pool; τ=τ, keepmask=km_c, T=T, nbins=nbins)

    begin
        fig = Figure(resolution=(820,380))
        ax  = Axis(fig[1,1], xlabel="Basal richness per cell",
                   ylabel="Pr( trophic height ≥ $(T) )",
                   title="Knee (probability form) at keep=$(keep_frac)")
        lines!(ax, xR, pR, label="Random")
        lines!(ax, xC, pC, label="Clustered")
        if isfinite(kR); vlines!(ax, [kR], color=:black, linestyle=:dash); end
        if isfinite(kC); vlines!(ax, [kC], color=:red,   linestyle=:dash);  end
        hlines!(ax, [0.5], color=(:gray,0.5), linestyle=:dot)
        axislegend(ax, position=:lt)
        display(fig)
    end

    return (; random_knee=kR, clustered_knee=kC)
end
