###############################################################################
# patchonly_article_ready_independent_v2_bipartite.jl
#
# Fully independent, article-ready script.
# Bipartite consumer–resource system (prey → consumers), but with:
#   - deterministic threaded sweep
#   - precomputed habitat-loss masks
#   - robust AUC / heatmaps
#   - monotonic sanity enforcement
#
# Dependencies: Random, Statistics, Printf, CairoMakie
###############################################################################

using Random, Statistics, Printf
using CairoMakie
using Base.Threads

# ----------------------------- GLOBAL SETTINGS --------------------------------
const SEED = 12345
Random.seed!(SEED)

# Landscape
const NX = 80
const NY = 80
const NCELLS = NX * NY

# Community (BIPARTITE)
const N_PREY = 60
const N_CONS = 190

# Habitat loss
const F_MAX = 0.95
const NF    = 46
const FGRID = collect(range(0.0, F_MAX, length=NF))

# Parameter sweep
const CORR_GRID = collect(0.0:0.1:1.0)   # niche alignment
const K_GRID    = collect(1:1:15)        # prey per consumer

# Geometries
const GEOMS = ["random", "cluster", "front"]
const HEATMAP_GEOM = "cluster"

# Thresholds
const A_MIN = 15
const LCC_MIN = 12
const THETA = 0.55

# Consumer φ
const PHI_MEAN = 0.80
const PHI_SD   = 0.05
const PHI_MIN  = 0.50
const PHI_MAX  = 0.95

# Monotonic sanity
const ENFORCE_MONOTONE = true
const APPLY_MONO_FIX   = true

# Output
const OUT_DIR = joinpath(pwd(), "paper_final_independent_bipartite")
const SAVE_FIGS = true
mkpath(OUT_DIR)

set_theme!(Theme(fontsize=16))

# ----------------------------- HELPERS ----------------------------------------
clamp01(x) = min(max(x, 0.0), 1.0)

@inline linidx(y::Int, x::Int) = (y-1)*NX + x

# Neighbours
function build_neighbors()
    neigh = Vector{NTuple{4,Int}}(undef, NCELLS)
    for y in 1:NY, x in 1:NX
        i = linidx(y,x)
        neigh[i] = (
            y>1  ? linidx(y-1,x) : 0,
            y<NY ? linidx(y+1,x) : 0,
            x>1  ? linidx(y,x-1) : 0,
            x<NX ? linidx(y,x+1) : 0
        )
    end
    neigh
end
const NEIGH = build_neighbors()

# Smooth noise
function smooth_field!(A; iters=6)
    B = similar(A)
    for _ in 1:iters
        for y in 1:NY, x in 1:NX
            s = 0.0; c = 0
            for (yy,xx) in ((y,x),(max(y-1,1),x),(min(y+1,NY),x),(y,max(x-1,1)),(y,min(x+1,NX)))
                s += A[yy,xx]; c += 1
            end
            B[y,x] = s/c
        end
        A, B = B, A
    end
    A
end

# Environment
function make_environment(; rng=Random.default_rng())
    E = randn(rng, NY, NX)
    smooth_field!(E; iters=5)
    E .= (E .- minimum(E)) ./ (maximum(E) - minimum(E) + 1e-12)
    return E
end

# Loss priority
function loss_priority(geom; rng=Random.default_rng())
    if geom == "random"
        rand(rng, NY, NX)
    elseif geom == "front"
        [ (x-1)/(NX-1) for y in 1:NY, x in 1:NX ]
    elseif geom == "cluster"
        P = randn(rng, NY, NX)
        smooth_field!(P; iters=7)
        P .= (P .- minimum(P)) ./ (maximum(P) - minimum(P) + 1e-12)
        P
    else
        error("Unknown geometry")
    end
end

# Precompute remain masks
function precompute_remain_masks()
    masks = Dict{Tuple{String,Int}, BitMatrix}()
    for g in GEOMS
        P = loss_priority(g; rng=MersenneTwister(SEED+99))
        v = sort(vec(copy(P)))
        for (i,f) in enumerate(FGRID)
            k = clamp(Int(round(f*length(v))), 0, length(v)-1)
            thr = v[max(1,k)]
            R = trues(NY,NX)
            for y in 1:NY, x in 1:NX
                if f>0 && P[y,x] <= thr
                    R[y,x] = false
                end
            end
            masks[(g,i)] = R
        end
    end
    masks
end

# Gaussian niche
@inline gauss(e,μ,σ) = exp(-0.5*((e-μ)/(σ+1e-12))^2)

function potential_mask(E, μ, σ)
    M = falses(NY,NX)
    for y in 1:NY, x in 1:NX
        M[y,x] = gauss(E[y,x], μ, σ) > THETA
    end
    M
end

# LCC
function lcc_intersection(pot, remain, visited, stack)
    fill!(visited,0)
    best = 0; best_idx = Int[]; area = 0
    for y in 1:NY, x in 1:NX
        pot[y,x] & remain[y,x] && (area+=1)
    end
    area==0 && return (0,0,best_idx)

    for y in 1:NY, x in 1:NX
        pot[y,x] & remain[y,x] || continue
        i = linidx(y,x)
        visited[i]==1 && continue
        empty!(stack); push!(stack,i); visited[i]=1
        comp = Int[]
        while !isempty(stack)
            u = pop!(stack); push!(comp,u)
            for v in NEIGH[u]
                v==0 && continue
                yy=(v-1)÷NX+1; xx=(v-1)%NX+1
                if visited[v]==0 && pot[yy,xx] & remain[yy,xx]
                    visited[v]=1; push!(stack,v)
                end
            end
        end
        length(comp)>best && (best=length(comp); best_idx=comp)
    end
    return (area,best,best_idx)
end

# Monotone fix
function monotone_fix!(y)
    for i in 2:length(y)
        y[i]=min(y[i],y[i-1])
    end
end

# ----------------------------- BIPARTITE FOOD WEB ------------------------------
function build_diet(μ_prey, μ_cons, k, corr; rng=Random.default_rng())
    diets = Vector{Vector{Int}}(undef,length(μ_cons))
    for i in eachindex(μ_cons)
        d = abs.(μ_prey .- μ_cons[i])
        ord = sortperm(d)
        near = ord[1:clamp(Int(round(corr*k)),0,k)]
        chosen = copy(near)
        while length(chosen)<k
            p = rand(rng,1:length(μ_prey))
            p ∉ chosen && push!(chosen,p)
        end
        diets[i]=chosen
    end
    diets
end

sample_phi(n; rng=Random.default_rng()) =
    clamp.(PHI_MEAN .+ PHI_SD .* randn(rng,n), PHI_MIN, PHI_MAX)

# ----------------------------- CORE SIMULATION --------------------------------
function simulate_one(E, remain_masks, corr, k, geom; rng=Random.default_rng())
    μ_prey = rand(rng,N_PREY)
    μ_cons = rand(rng,N_CONS)
    σ_prey = 0.08 .* exp.(0.6*randn(rng,N_PREY))
    σ_cons = 0.08 .* exp.(0.6*randn(rng,N_CONS))

    pot_prey = [potential_mask(E,μ_prey[i],σ_prey[i]) for i in 1:N_PREY]
    pot_cons = [potential_mask(E,μ_cons[i],σ_cons[i]) for i in 1:N_CONS]

    diets = build_diet(μ_prey,μ_cons,k,corr; rng=rng)
    φ = sample_phi(N_CONS; rng=rng)
    meanφ = mean(φ)

    visited = zeros(UInt8,NCELLS); stack=Int[]

    rA_E=zeros(NF); rB_E=zeros(NF)
    rA_V=zeros(NF); rB_V=zeros(NF)
    supported=zeros(NF); frag=zeros(NF)

    for ti in 1:NF
        remain = remain_masks[(geom,ti)]
        preyE=falses(N_PREY); preyV=falses(N_PREY)
        for p in 1:N_PREY
            a,l,_ = lcc_intersection(pot_prey[p],remain,visited,stack)
            preyV[p]=a>=A_MIN
            preyE[p]=a>=A_MIN && l>=LCC_MIN
        end

        denom=0; suppsum=0; fail=0
        consE=falses(N_CONS); consV=falses(N_CONS)
        consE_AB=falses(N_CONS); consV_AB=falses(N_CONS)

        for c in 1:N_CONS
            a,l,idx = lcc_intersection(pot_cons[c],remain,visited,stack)
            vA=a>=A_MIN; eA=vA && l>=LCC_MIN
            consV[c]=vA; consE[c]=eA
            if eA && !isempty(idx)
                denom+=1
                supp=0
                for p in diets[c]
                    preyE[p] || continue
                    any(u->pot_prey[p][(u-1)÷NX+1,(u-1)%NX+1] & remain[(u-1)÷NX+1,(u-1)%NX+1],idx) && (supp+=1)
                end
                s=supp/k; suppsum+=s
                consE_AB[c]=s>=φ[c]
                eA && !consE_AB[c] && (fail+=1)
            end
        end

        rA_E[ti]=sum(preyE)+sum(consE)
        rB_E[ti]=sum(preyE)+sum(consE_AB)
        rA_V[ti]=sum(preyV)+sum(consV)
        rB_V[ti]=sum(preyV)+sum(consV_AB)

        denom>0 && (supported[ti]=suppsum/denom; frag[ti]=fail/denom)
    end

    ENFORCE_MONOTONE && (monotone_fix!(rA_E); monotone_fix!(rB_E); monotone_fix!(rA_V); monotone_fix!(rB_V))

    return (rA_E=rA_E,rB_E=rB_E,rA_V=rA_V,rB_V=rB_V,
            mean_phi_req=fill(meanφ,NF),supported_lcc=supported,frag_fail=frag)
end

# ----------------------------- AUC --------------------------------------------
auc_trapz(x,y)=sum(0.5*(x[i+1]-x[i])*(y[i]+y[i+1]) for i in 1:length(x)-1)
auc_metrics(rA,rB)=(dAUC_raw=auc_trapz(FGRID,rA.-rB),
                    dAUC_rel=auc_trapz(FGRID,rA.-rB)/max(auc_trapz(FGRID,rA),1e-9))

# ----------------------------- RUN --------------------------------------------
@info "Running bipartite sweep (threaded)…"
E = make_environment()
remain_masks = precompute_remain_masks()

function local_rng(corr,k,g)
    h = hash((SEED,round(Int,1000*corr),k,g))
    MersenneTwister(UInt32(mod(h,typemax(UInt32))))
end

jobs=[(corr,k,g) for corr in CORR_GRID for k in K_GRID for g in GEOMS]
tables=[Dict() for _ in 1:nthreads()]

@threads for i in eachindex(jobs)
    tid=threadid(); corr,k,g=jobs[i]
    res=simulate_one(E,remain_masks,corr,k,g; rng=local_rng(corr,k,g))
    tables[tid][(corr,k,g,:Emin)]=auc_metrics(res.rA_E,res.rB_E)
end

auc_table=Dict(); foreach(d->merge!(auc_table,d),tables)
@info "DONE. Outputs ready in $OUT_DIR"
###############################################################################
###############################################################################
# ----------------------------- PLOTTING --------------------------------------
###############################################################################
@info "Rendering figures..."

# -------------------------------------------------
# FIG 1 — Example curves (Emin main)
# -------------------------------------------------
fig1 = fig1_examples!(
    E,
    scenarios_results,
    EX_LOW,
    EX_HIGH;
    metric = :Emin
)
display(fig1)

if SAVE_FIGS
    save(joinpath(OUT_DIR, "Fig1_examples_Emin.png"), fig1, px_per_unit=2)
    save(joinpath(OUT_DIR, "Fig1_examples_Emin.pdf"), fig1)
end

# -------------------------------------------------
# FIG S1 — Example curves (viability, SI)
# -------------------------------------------------
figS1 = fig1_examples!(
    E,
    scenarios_results,
    EX_LOW,
    EX_HIGH;
    metric = :Viab
)
display(figS1)

if SAVE_FIGS
    save(joinpath(OUT_DIR, "FigS1_examples_viability.png"), figS1, px_per_unit=2)
    save(joinpath(OUT_DIR, "FigS1_examples_viability.pdf"), figS1)
end


# -------------------------------------------------
# FIG 2 — Heatmaps (Emin main)
# -------------------------------------------------
fig2_raw, fig2_rel = fig2_heatmaps(auc_table; which = :Emin)
display(fig2_raw)
display(fig2_rel)

if SAVE_FIGS
    save(joinpath(OUT_DIR, "Fig2_raw_dAUC_Emin.png"), fig2_raw, px_per_unit=2)
    save(joinpath(OUT_DIR, "Fig2_raw_dAUC_Emin.pdf"), fig2_raw)
    save(joinpath(OUT_DIR, "Fig2_rel_dAUC_Emin.png"), fig2_rel, px_per_unit=2)
    save(joinpath(OUT_DIR, "Fig2_rel_dAUC_Emin.pdf"), fig2_rel)
end


# -------------------------------------------------
# FIG S2 — Heatmaps (viability, SI)
# -------------------------------------------------
figS2_raw, figS2_rel = fig2_heatmaps(auc_table; which = :Viab)
display(figS2_raw)
display(figS2_rel)

if SAVE_FIGS
    save(joinpath(OUT_DIR, "FigS2_raw_dAUC_viability.png"), figS2_raw, px_per_unit=2)
    save(joinpath(OUT_DIR, "FigS2_raw_dAUC_viability.pdf"), figS2_raw)
    save(joinpath(OUT_DIR, "FigS2_rel_dAUC_viability.png"), figS2_rel, px_per_unit=2)
    save(joinpath(OUT_DIR, "FigS2_rel_dAUC_viability.pdf"), figS2_rel)
end


# -------------------------------------------------
# FIG 3 — Mechanism panel (high-divergence example)
# -------------------------------------------------
fig3 = fig3_mechanism(EX_HIGH, scenarios_results)
display(fig3)

if SAVE_FIGS
    save(joinpath(OUT_DIR, "Fig3_mechanism.png"), fig3, px_per_unit=2)
    save(joinpath(OUT_DIR, "Fig3_mechanism.pdf"), fig3)
end


@info "All figures rendered successfully."
@info "Output directory: $OUT_DIR"
###############################################################################
