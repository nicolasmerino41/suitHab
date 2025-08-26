# =========================================
# Minimal habitat-loss demo (Makie-friendly)
# Rule-based assembly: climate filter → ≥1 prey
# No dispersal, no K, δ-free (structural) demo
# =========================================
# ---------- Helpers: grid & indices ----------
struct Grid
    nx::Int
    ny::Int
    C::Int
    xy::Matrix{Float64}  # 2 x C (columns are cells)
    climate::Vector{Float64} # length C
end

# Build a simple grid with a smooth climate gradient (0..1)
function make_grid(nx::Int=50, ny::Int=50; seed=42)
    Random.seed!(seed)
    xs = range(0, 1; length=nx)
    ys = range(0, 1; length=ny)
    xy = Array{Float64}(undef, 2, nx*ny)
    clim = Vector{Float64}(undef, nx*ny)
    k = 1
    for j in 1:ny, i in 1:nx
        xy[:,k] = [xs[i], ys[j]]
        # climate = diagonal gradient + mild sinusoidal texture
        clim[k] = 0.7*xs[i] + 0.3*ys[j] + 0.1*sin(6π*xs[i])*sin(6π*ys[j])
        k += 1
    end
    # normalize to 0..1
    clim .-= minimum(clim); clim ./= (maximum(clim) + eps())
    Grid(nx, ny, nx*ny, xy, clim)
end

# Neighbor indices (4-neighborhood) for clustered masks
function neighbors4(ix::Int, nx::Int, ny::Int)
    i = ((ix - 1) % nx) + 1
    j = ((ix - 1) ÷ nx) + 1
    out = Int[]
    if i > 1;       push!(out, ix - 1);     end
    if i < nx;      push!(out, ix + 1);     end
    if j > 1;       push!(out, ix - nx);    end
    if j < ny;      push!(out, ix + nx);    end
    out
end

# ---------- Species, niches, metaweb ----------
struct SpeciesPool
    S::Int
    masses::Vector{Float64}
    basal::BitVector
    mu::Vector{Float64}   # climate niche centers in [0,1]
    b::Vector{Float64}    # niche breadths
    E::Vector{Vector{Int}}# adjacency list: E[s] = prey indices for predator s
end

# Build a simple allometric metaweb via mass ratio window (DAG by construction)
function build_pool(S::Int=140; basal_frac=0.45, seed=1,
                    Rmin=3.0, Rmax=300.0,  # wider diet window
                    b0=0.12, bspread=0.08) # broader niches
    Random.seed!(seed)

    # masses: log-spread + shuffle
    logm = collect(range(log(1e-2), log(10.0); length=S))
    shuffle!(logm)
    masses = exp.(logm)

    # order by mass (increasing)
    order = sortperm(masses)

    # >>> basal are the LIGHTEST nB species <<<
    nB = round(Int, basal_frac*S)
    basal = falses(S)
    basal[order[1:nB]] .= true

    # climate niches
    mu = rand(S)
    b  = b0 .+ bspread .* rand(S)

    # metaweb via mass-ratio window; ensure each consumer has ≥1 prey
    E = [Int[] for _ in 1:S]
    for ii in eachindex(order)
        s = order[ii]
        if basal[s]; continue; end
        for jj in 1:ii-1
            q = order[jj]
            r = masses[s] / masses[q]
            if r ≥ Rmin && r ≤ Rmax
                push!(E[s], q)
            end
        end
        # fallback: if empty, link to the immediate smaller species
        if isempty(E[s]) && ii > 1
            push!(E[s], order[ii-1])
        end
    end

    return SpeciesPool(S, masses, basal, mu, b, E)
end

# ---------- Climate pass (Z) ----------
# Suitability (Gaussian in climate); pass if ≥ τ
function climate_pass(pool::SpeciesPool, grid::Grid; τ=0.5)
    S, C = pool.S, grid.C
    Z = falses(S, C)
    @inbounds for s in 1:S, c in 1:C
        su = exp(-((grid.climate[c] - pool.mu[s])^2) / (2 * pool.b[s]^2 + 1e-12))
        if su ≥ τ
            Z[s,c] = true
        end
    end
    Z
end

# ---------- Assembly: prune unsupported consumers until fixed point ----------
function assemble(Z::BitMatrix, pool::SpeciesPool)
    S, C = size(Z)
    P = copy(Z) # start with climate pass
    changed = true
    while changed
        changed = false
        @inbounds for c in 1:C, s in 1:S
            if P[s,c] && !pool.basal[s]
                prey = pool.E[s]
                has_pre = false
                @inbounds for q in prey
                    if P[q,c]; has_pre = true; break; end
                end
                if !has_pre
                    P[s,c] = false
                    changed = true
                end
            end
        end
    end
    P
end

# ---------- Metrics ----------
# Biotically Safe Climate Habitat (BSH, r=1) per species (fraction of ALL cells)
function bsh1_per_species(P::BitMatrix, Z::BitMatrix, pool::SpeciesPool)
    S, C = size(P)
    out = zeros(Float64, S)
    @inbounds for s in 1:S
        if isempty(pool.E[s]) # basal: BSH = climate-suitable area
            out[s] = count(Z[s, :]) / C
        else
            ok = 0
            for c in 1:C
                if Z[s,c]
                    # has at least 1 prey present?
                    prey_ct = 0
                    for q in pool.E[s]
                        prey_ct += P[q,c] ? 1 : 0
                        if prey_ct > 0; break; end
                    end
                    ok += (prey_ct > 0) ? 1 : 0
                end
            end
            out[s] = ok / C
        end
    end
    out
end

# ILRF per species: share of climate-suitable cells where prey are absent
function ilrf_per_species(P::BitMatrix, Z::BitMatrix, pool::SpeciesPool)
    S, C = size(P)
    out = fill(Float64(NaN), S)
    @inbounds for s in 1:S
        suitable = 0
        no_prey  = 0
        if pool.basal[s]
            out[s] = 0.0
        else
            for c in 1:C
                if Z[s,c]
                    suitable += 1
                    prey_ct = 0
                    for q in pool.E[s]
                        prey_ct += P[q,c] ? 1 : 0
                        if prey_ct > 0; break; end
                    end
                    if prey_ct == 0; no_prey += 1; end
                end
            end
            out[s] = (suitable == 0) ? NaN : no_prey / suitable
        end
    end
    out
end

# Mean over consumers only
consumer_mask(pool::SpeciesPool) = .!pool.basal

# ---------- Habitat-loss masks ----------
# Apply a keep-mask (true = keep cell; false = removed)
apply_mask(Z::BitMatrix, keep::BitVector) = Z[:, keep]

# Random keep mask for given keep_frac
function random_mask(C::Int, keep_frac::Float64; seed::Int=0)
    Random.seed!(seed)
    keep = falses(C)
    nkeep = round(Int, keep_frac * C)
    keep[randperm(C)[1:nkeep]] .= true
    keep
end

# Clustered removal via BFS growth from multiple seeds until target reached
function clustered_mask(grid::Grid, keep_frac::Float64; nseeds::Int=6, seed::Int=0)
    Random.seed!(seed)
    C = grid.C
    target_remove = C - round(Int, keep_frac*C)
    removed = falses(C)
    q = Int[]
    # pick initial seeds
    seeds = randperm(C)[1:nseeds]
    append!(q, seeds)
    removed_ct = 0
    ptr = 1
    while removed_ct < target_remove
        if ptr > length(q)
            # start a new seed in a not-yet-removed cell
            newseed = first(filter(i -> !removed[i], randperm(C)))
            push!(q, newseed)
        end
        v = q[ptr]; ptr += 1
        if removed[v]; continue; end
        removed[v] = true
        removed_ct += 1
        # enqueue neighbors
        for nb in neighbors4(v, grid.nx, grid.ny)
            if !removed[nb]
                push!(q, nb)
            end
        end
    end
    .!removed
end

# ---------- Experiment runner ----------
struct RunResult
    meanBSH_consumers::Float64
    meanILRF_consumers::Float64
end

function one_run(pool::SpeciesPool, grid::Grid; τ=0.5, keep_mask=nothing)
    Zfull = climate_pass(pool, grid; τ=τ)
    Z = isnothing(keep_mask) ? Zfull : apply_mask(Zfull, keep_mask)
    P = assemble(Z, pool)
    cons = consumer_mask(pool)
    BSH = bsh1_per_species(P, Z, pool)
    ILR = ilrf_per_species(P, Z, pool)
    RunResult(mean(skipmissing(BSH[cons])), mean(skipmissing(ILR[cons])))
end

# ---------- Demo workflow ----------
grid = make_grid(60, 60; seed=11)
heatmap(grid)
pool = build_pool(140; basal_frac=0.4, seed=2)

# Baseline (no loss)
base = one_run(pool, grid; τ=0.5)

# Sweep equal-area loss fractions, compare random vs clustered
loss_fracs = 0.0:0.05:0.6
ΔBSH_random    = Float64[]
ΔBSH_clustered = Float64[]
for (k, f) in enumerate(loss_fracs)
    keep = 1.0 - f
    # random
    km_r = random_mask(grid.C, keep; seed=100 + k)
    r = one_run(pool, grid; τ=0.5, keep_mask=km_r)
    push!(ΔBSH_random,    r.meanBSH_consumers - base.meanBSH_consumers)
    # clustered
    km_c = clustered_mask(grid, keep; nseeds=6, seed=200 + k)
    c = one_run(pool, grid; τ=0.5, keep_mask=km_c)
    push!(ΔBSH_clustered, c.meanBSH_consumers - base.meanBSH_consumers)
end

# ---------- Plot (Makie pattern you prefer) ----------
begin
    fig = Figure(; size = (800, 420))
    ax = Axis(fig[1,1], xlabel = "Area lost (fraction)", ylabel = "ΔBSH (mean over consumers)")
    lines!(ax, collect(loss_fracs), ΔBSH_random, label = "Random loss")
    lines!(ax, collect(loss_fracs), ΔBSH_clustered, label = "Clustered loss")
    axislegend(ax, position = :lb)
    hlines!(ax, [0.0], color = (:gray, 0.5), linestyle = :dash)
    display(fig)
end

# ---------- Optional: quick sanity print ----------
println("Baseline mean BSH (consumers): ", round(base.meanBSH_consumers, digits=3))
println("At 40% loss, ΔBSH random vs clustered: ",
        round(ΔBSH_random[findfirst(x->x≈0.4, loss_fracs)], digits=3), " vs ",
        round(ΔBSH_clustered[findfirst(x->x≈0.4, loss_fracs)], digits=3))

        
