module PlotAbsBSH

using Random, Statistics
using CairoMakie

# -- we only need a minimal interface:
# your code must provide:
#   make_grid(...) -> grid with fields .C  (cell count)
#   build_pool(S; kwargs...) -> pool with fields:
#       S::Int, basal::BitVector, E::Vector{Vector{Int}}
#   climate_pass(pool, grid; τ=0.5) -> BitMatrix Z  (S × C)
#   (AM) assemble(Z, pool) -> BitMatrix P (S × C)   # your existing AM assemble
#   (utility) bsh1_count(P, Z, pool) -> Vector{Int} # counts of BSH cells per species

# ---- BAM assemble with τB (fraction-of-diet present)
# passes a consumer in cell c if (eligible prey / total diet) ≥ τB
function assemble_BAM(Z::BitMatrix, pool; τB::Float64)
    S, C = size(Z)
    P = copy(Z)          # start from climate pass
    changed = true
    while changed
        changed = false
        @inbounds for c in 1:C, s in 1:S
            pool.basal[s] && continue
            if P[s,c]
                prey = pool.E[s]; np = length(prey)
                # total diet size; if someone has 0 (shouldn’t), treat as fail-safe pass
                if np == 0; continue; end
                # count prey that are present in c (under current P)
                elig = 0
                @inbounds for q in prey
                    elig += (P[q,c] ? 1 : 0)
                end
                if (elig / np) < τB
                    P[s,c] = false
                    changed = true
                end
            end
        end
    end
    return P
end

# ---- absolute BSH (per-consumer mean) normalized by ORIGINAL area
# C_full = grid.C even after masking
function mean_BSH_area(P::BitMatrix, Z::BitMatrix, pool, C_full::Int)
    counts = bsh1_count(P, Z, pool)                 # you already have this
    return mean(counts[.!pool.basal]) / C_full
end

# ---- masks
apply_mask(Z::BitMatrix, keep::BitVector) = Z[:, keep]

function random_mask(C::Int, keep_frac::Float64; seed::Int=0)
    Random.seed!(seed)
    keep = falses(C)
    nkeep = round(Int, keep_frac * C)
    keep[randperm(C)[1:nkeep]] .= true
    keep
end

# 4-neighborhood neighbors (helper for clustered mask)
neighbors4(ix::Int, nx::Int, ny::Int) = begin
    i = ((ix - 1) % nx) + 1; j = ((ix - 1) ÷ nx) + 1
    out = Int[]
    if i > 1; push!(out, ix-1); end
    if i < nx; push!(out, ix+1); end
    if j > 1; push!(out, ix - nx); end
    if j < ny; push!(out, ix + nx); end
    out
end

function clustered_mask(grid, keep_frac; nseeds::Int=6, seed::Int=0)
    Random.seed!(seed)
    C = grid.C
    target_remove = C - round(Int, keep_frac*C)
    removed = falses(C)
    q = Int[]
    seeds = randperm(C)[1:nseeds]; append!(q, seeds)
    removed_ct = 0; ptr = 1
    while removed_ct < target_remove
        if ptr > length(q)
            newseed = first(filter(i->!removed[i], randperm(C)))
            push!(q, newseed)
        end
        v = q[ptr]; ptr += 1
        if removed[v]; continue; end
        removed[v] = true; removed_ct += 1
        for nb in neighbors4(v, grid.nx, grid.ny)
            if !removed[nb]; push!(q, nb); end
        end
    end
    .!removed
end

# front-like: keep the “cold” (or “hot”) side fraction
function front_mask(grid, keep_frac; cold_side::Bool=true)
    C = grid.C
    order = sortperm(grid.climate; rev=!cold_side)  # choose side
    keep = falses(C)
    nkeep = round(Int, keep_frac*C)
    keep[order[1:nkeep]] .= true
    keep
end

# ---- main: compute curves for AM and BAM (absolute, vs original area)
function abs_bsh_curves(; grid, S::Int, pool_kwargs=NamedTuple(),
        τA::Float64=0.5, τB::Float64=0.5,
        loss_fracs=0.2:0.1:0.8, geoms=(:random,:clustered,:front),
        pool_seed::Int=1, mask_seeds=1:30, nseeds_cluster::Int=6)

    # one pool for this figure (you can ensemble if you want later)
    pool  = build_pool(S; seed=pool_seed, pool_kwargs...)
    Zfull = climate_pass(pool, grid; τ=τA)
    P0_AM  = assemble(Zfull, pool)
    P0_BAM = assemble_BAM(Zfull, pool; τB=τB)
    Cfull = grid.C

    # struct: geom => (f => (AM_mean, BAM_mean))
    out = Dict{Symbol,NamedTuple}()

    for g in geoms
        AMm = Float64[]; BAMm = Float64[]
        for f in loss_fracs
            keep = 1.0 - f
            vals_AM  = Float64[]
            vals_BAM = Float64[]
            for ms in mask_seeds
                keepmask = g === :random   ? random_mask(Cfull, keep; seed=ms) :
                           g === :clustered ? clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=ms) :
                           g === :front    ? front_mask(grid, keep) :
                                             error("unknown geom $g")
                Z = apply_mask(Zfull, keepmask)
                if size(Z,2) == 0
                    push!(vals_AM, 0.0); push!(vals_BAM, 0.0); continue
                end
                P_AM  = assemble(Z, pool)
                P_BAM = assemble_BAM(Z, pool; τB=τB)
                push!(vals_AM,  mean_BSH_area(P_AM,  Z, pool, Cfull))
                push!(vals_BAM, mean_BSH_area(P_BAM, Z, pool, Cfull))
            end
            push!(AMm,  mean(vals_AM))
            push!(BAMm, mean(vals_BAM))
        end
        out[g] = (loss=collect(loss_fracs), AM=AMm, BAM=BAMm)
    end

    return out
end

# ---- simple plotting helper
function plot_abs_bsh(out::Dict; title::String="")
    fig = Figure(; size=(1050,350))
    geoms = collect(keys(out))
    for (j,g) in enumerate(geoms)
        ax = Axis(fig[1,j], xlabel="area lost (fraction)", ylabel="BSH (mean over consumers / original area)",
                  title=String(g))
        dat = out[g]
        lines!(ax, dat.loss, dat.AM;  linestyle=:dash,  linewidth=3, color=:gray35, label="AM")
        lines!(ax, dat.loss, dat.BAM; linestyle=:solid, linewidth=3, color=:black, label="BAM")
        hlines!(ax, [0], color=(:gray,0.4), linestyle=:dot)
        axislegend(ax, position=:lb)
    end
    Label(fig[0,1:length(geoms)], title, fontsize=18, tellwidth=false)
    return fig
end

end # module
