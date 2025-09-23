# ----------------------- run_coretention_targeting.jl -------------------------
# Narrative 2: A simple co-retention-greedy rule beats random/clustered and is
# competitive with 'front' & A-greedy at the same budget.
# Figures:
#   F1_mean_area_vs_f.png, F2_tailP10_vs_f.png, F3_lift_over_Agreedy.png
# ------------------------------------------------------------------------------

# --- modules
include("../SetUp.jl")
include("src/grids.jl");      using .Grids
include("src/grids_init.jl"); using .GridsInit
include("src/metawebs.jl");   using .Metawebs
include("src/hl.jl");         using .HL
include("src/bsh.jl");        using .BSH

# --- settings
rng        = MersenneTwister(123)
nx, ny     = 120, 120
S          = 200
basal_frac = 0.25
loss_fracs = 0.20:0.05:0.80
outdir     = "Paper/narratives/coretention"; isdir(outdir) || mkpath(outdir)

# pick two grids to show generality
include("src/grids_init.jl"); using .GridsInit
G   = GridsInit.build_default_grids(nx, ny; seeds=(11,12,13,14))
grid = G.grad     # headline (you can loop over more if you want)

# pool & params
pool = Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:mid)
pars = BSH.BAMParams(; τA=0.5, τB=0.5, τocc=0.2, γ=3.0, movement=:component, T=8)

# --- baseline A and AM on full grid (used to score)
A      = BSH.abiotic_matrix(pool, grid; seed=1)
keep0  = trues(grid.C)
P_AM0, _ = BSH.assemble_AM(pool, grid, A, keep0; pars=pars)   # prey presence proxy on intact landscape
cons_inds = findall(!, pool.basal)

# --- scoring functions (vector over cells)
function score_Agreedy(pool, grid, A; τA)
    # count how many consumers are climatically eligible in cell i
    C = grid.C; sc = zeros(Float64, C)
    for i in 1:C
        sc[i] = mean(A[cons_inds, i] .>= τA)
    end
    sc
end

function score_CoRetention(pool, grid, A, P_AM0)
    # for each consumer s and cell i: fraction of its prey present in i under AM on intact landscape
    C = grid.C; sc = zeros(Float64, C)
    for i in 1:C
        acc = 0.0
        for s in cons_inds
            pr = pool.prey[s]
            if isempty(pr); continue; end
            acc += mean(P_AM0[pr, i])   # [0,1]
        end
        sc[i] = acc / length(cons_inds)
    end
    sc
end

scoreA  = score_Agreedy(pool, grid, A; τA=pars.τA)
scoreCR = score_CoRetention(pool, grid, A, P_AM0)

# --- strategies -> keep masks
function keep_by_score(score::AbstractVector{<:Real}, frac_keep::Float64)
    C = length(score); nkeep = clamp(round(Int, frac_keep*C), 0, C)
    idx = sortperm(score; rev=true)[1:nkeep]
    keep = falses(C); keep[idx] .= true; keep
end

function keep_by_geom(geom::Symbol, grid, frac_keep; rng=MersenneTwister(1))
    g = geom
    if g===:random
        return HL.random_mask(rng, grid.C, frac_keep)
    elseif g===:clustered
        return HL.clustered_mask(rng, grid.nx, grid.ny, frac_keep; nseeds=8)
    else
        return HL.front_mask(rng, grid.xy, frac_keep; axis=:x, noise=0.05)
    end
end

# --- evaluation: for each f and strategy compute (i) mean BAM area, (ii) P10 of per-species area
strategies = [:random, :clustered, :front, :Agreedy, :Coretention]
labels     = Dict(:random=>"random", :clustered=>"clustered", :front=>"front",
                  :Agreedy=>"A-greedy", :Coretention=>"co-retention")

meanBAM = Dict(s => Float64[] for s in strategies)
P10     = Dict(s => Float64[] for s in strategies)

for f in loss_fracs
    frac_keep = 1 - f
    keeps = Dict{Symbol,BitVector}()

    keeps[:random]     = keep_by_geom(:random,    grid, frac_keep; rng=MersenneTwister(10))
    keeps[:clustered]  = keep_by_geom(:clustered, grid, frac_keep; rng=MersenneTwister(11))
    keeps[:front]      = keep_by_geom(:front,     grid, frac_keep; rng=MersenneTwister(12))
    keeps[:Agreedy]    = keep_by_score(scoreA,  frac_keep)
    keeps[:Coretention]= keep_by_score(scoreCR, frac_keep)

    for s in strategies
        keep = keeps[s]
        bam = BSH.assemble_BAM(pool, grid, A, keep; pars=pars)
        # (i) mean over consumers / original area
        push!(meanBAM[s], mean(bam.bsh[cons_inds]))
        # (ii) P10 of per-species retained area
        per = [sum(@view bam.P[s, :]) / grid.C for s in cons_inds]
        push!(P10[s], quantile(per, 0.10))
    end
end

# --- Fig 1: mean BAM area vs f (all strategies)
begin
    f1 = Figure(; size=(1000,360))
    ax1 = Axis(
        f1[1,1], title="Mean BAM retained area vs f",
        xlabel="area lost f", ylabel="retained area"
    )
    cols = Dict(:random=>:steelblue, :clustered=>:orange, :front=>:green, :Agreedy=>:black, :Coretention=>:purple)
    for s in strategies
        lines!(ax1, collect(loss_fracs), meanBAM[s], label=labels[s], color=cols[s])
    end
    axislegend(ax1, position=:rt, nbanks=2)
    display(f1)
end
save(joinpath(outdir,"F1_mean_area_vs_f.png"), f1)

# --- Fig 2: tail P10 vs f
begin
    f2 = Figure(; size=(1000,360))
    ax2 = Axis(
        f2[1,1], title="P10 of per-species retained area vs f",
        xlabel="area lost f", ylabel="P10"
    )
    for s in strategies
        lines!(ax2, collect(loss_fracs), P10[s], label=labels[s], color=cols[s])
    end
    axislegend(ax2, position=:rt, nbanks=2)
    display(f2)
end
save(joinpath(outdir,"F2_tailP10_vs_f.png"), f2)

# --- Fig 3: lift of co-retention over A-greedy (mean area)
lift = [meanBAM[:Coretention][k] - meanBAM[:Agreedy][k] for k in eachindex(loss_fracs)]
begin
    f3 = Figure(; size=(800,350))
    ax3 = Axis(
        f3[1,1], title="Lift: co-retention − A-greedy",
        xlabel="area lost f", ylabel="Δ retained area"
    )
    lines!(ax3, collect(loss_fracs), lift, color=:purple)
    hlines!(ax3, [0.0], color=:gray, linestyle=:dash)
    display(f3)
end
save(joinpath(outdir,"F3_lift_over_Agreedy.png"), f3)

@info "Co-retention narrative figs saved in $outdir"
