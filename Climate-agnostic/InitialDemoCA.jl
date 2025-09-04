# ============ 1) Patchy spatial fields (no packages) ============
# Build smooth 'hotspot' fields as sums of radial bumps.
########## 0) Utilities ##########
function rbf_field(W, H; ncenters=6, scale=10.0, seed=1)
    rng = MersenneTwister(seed)
    xs = collect(0:W-1); ys = collect(0:H-1)
    X = repeat(xs', H, 1); Y = repeat(ys, 1, W)
    cx = rand(rng, 0:W-1, ncenters); cy = rand(rng, 0:H-1, ncenters)
    w  = rand(rng, ncenters) .+ 0.5
    σ2 = scale^2
    F = zeros(Float64, H, W)
    @inbounds for k in 1:ncenters
        F .+= w[k] .* @. exp(-((X - cx[k])^2 + (Y - cy[k])^2) / (2σ2))
    end
    F .-= minimum(F); F ./= maximum(F) + 1e-12
    return F
end

# Convenience: flatten H×W to C vector aligned with your Grid order
flatten_field(F) = vec(permutedims(F, (2,1)))  # (W fast) → C = W*H

# ============ 2) Community-first assembly ============

"""
Build a pool with consumers tied to prey guilds (low redundancy, high synchrony).
- S: total spp; basal_frac controls # basal spp
- G: number of basal guilds
- diet_k: typical #prey per consumer (small → low redundancy)
- sync_bias in [0,1]: prob a consumer draws all prey from a *single* guild (synchrony)
Returns (pool, guild_of_basal::Vector{Int})
"""
function build_pool_cf(S; basal_frac=0.35, G=3, diet_k=3, sync_bias=0.8, seed=1)
    rng = MersenneTwister(seed)
    nb = round(Int, basal_frac*S); nc = S - nb
    basal = falses(S); basal[1:nb] .= true

    g_basal = [rand(rng, 1:G) for _ in 1:nb]              # guild of each basal sp

    # diets
    E = [Int[] for _ in 1:S]
    for s in nb+1:S
        if rand(rng) < sync_bias
            g = rand(rng, 1:G)
            cand = findall(i->i==g, g_basal)
            E[s] = rand(rng, cand, min(diet_k, length(cand)))
        else
            E[s] = rand(rng, 1:nb, min(diet_k, nb))
        end
    end

    mu = fill(0.5, S); b = fill(0.2, S)
    masses = range(0.5, 10.0; length=S) |> collect
    return SpeciesPool(S, masses, basal, mu, b, E), g_basal
end

"""
Assemble communities from fields (no climate used yet):
- Fabi: abiotic field in [0,1]
- Fg[g]: basal guild fields in [0,1] (vector of H×W matrices)
Rules:
  basal i in guild g present if Fg[g] ≥ t_basal
  consumer present if (# present prey) / diet_k ≥ prey_ratio_min AND Fabi ≥ t_cons
Returns P (S×C BitMatrix) and climate vector Aclim0 = Fabi (flattened, 0..1).
"""
function assemble_from_fields(pool::SpeciesPool, Fabi::Matrix{Float64}, Fg::Vector{Matrix{Float64}};
                              g_basal::Vector{Int}, t_basal=0.55, t_cons=0.35, prey_ratio_min=0.5)
    H, W = size(Fabi); C = W*H
    P = falses(pool.S, C)
    Aclim0 = flatten_field(Fabi)

    nb = count(pool.basal)
    # basal presence by guild threshold
    for (i,g) in enumerate(g_basal)
        present = flatten_field(Fg[g] .>= t_basal)
        P[i, :] .= present
    end
    # consumers
    for s in nb+1:pool.S
        prey = pool.E[s]; diet = max(1, length(prey))
        @inbounds for c in 1:C
            if Aclim0[c] ≥ t_cons
                k = 0; for q in prey; k += P[q,c] ? 1 : 0; end
                if k / diet ≥ prey_ratio_min
                    P[s,c] = true
                end
            end
        end
    end
    return P, Aclim0
end

# ============ 3) Fit "climate niche" from occupancy ============
# Simple Gaussian fit to occupancy vs abiotic field (using mean & sd of used cells)
function fit_niche_from_P(Ps::BitVector, Aclim0::Vector{Float64}; min_b=0.04)
    vals = Aclim0[Ps .== true]
    if isempty(vals); return (mu=0.5, b=min_b); end
    μ = clamp(mean(vals), 0, 1); σ = max(std(vals), min_b)
    return (mu=μ, b=σ)
end

"""
Infer niches for all spp from P & abiotic field, then build Z0 (S×C BitMatrix)
with a simple threshold: |A - μ| ≤ k*b
"""
function climate_from_fitted(P::BitMatrix, Aclim0::Vector{Float64}, pool::SpeciesPool; k=1.2)
    S, C = size(P)
    μ = similar(pool.mu); b = similar(pool.b)
    for s in 1:S
        t = fit_niche_from_P(BitVector(view(P,s,:)), Aclim0)
        μ[s] = t.mu; b[s] = t.b
    end
    Z0 = falses(S, C)
    @inbounds for s in 1:S, c in 1:C
        Z0[s,c] = abs(Aclim0[c] - μ[s]) ≤ k*b[s]
    end
    return Z0, μ, b
end

# ============ 4) Viability decomposition (reuse) ============
# Use the viability-based decomp you already have; pasting for completeness:
viability_bits(A,p;A_min=0.10,p_min=0.10) = map(x->x[1] ≥ A_min && x[2] ≥ p_min,
                                                zip(A,p)) |> BitVector
mean_viability(A,p;A_min=0.10,p_min=0.10) = mean(viability_bits(A,p;A_min=A_min,p_min=p_min))

# --- NEW: viability decomposition that takes a precomputed climate mask Z0 ---
"""
decomp_viability_fromZ(Z0, pool; keepmask, A_min, p_min)
- Z0: S×C BitMatrix, the 'before' climate-suitable mask YOU fitted
- keepmask: C-length BitVector, kept cells (true = kept)
Returns (clim, inter, syn, U0, U1)
"""
function decomp_viability_fromZ(Z0::BitMatrix, pool::SpeciesPool;
                                keepmask::BitVector, A_min::Float64, p_min::Float64)
    # assemble with this Z0 (not calling climate_pass)
    P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask)
    P1 = assemble(Z1, pool)

    cons = .!pool.basal
    A0 = climate_area_per_species(Z0)[cons]
    A1 = climate_area_per_species(Z1)[cons]
    p0 = bsh1_cond_per_species(P0, Z0, pool)[cons]
    p1 = bsh1_cond_per_species(P1, Z1, pool)[cons]

    # undefined conditional → 0
    for x in (p0,p1); @inbounds for i in eachindex(x); if isnan(x[i]); x[i]=0.0; end; end; end

    V00 = mean(viability_bits(A0,p0;A_min=A_min,p_min=p_min))
    V10 = mean(viability_bits(A1,p0;A_min=A_min,p_min=p_min))   # climate-only change
    V01 = mean(viability_bits(A0,p1;A_min=A_min,p_min=p_min))   # interaction-only change
    V11 = mean(viability_bits(A1,p1;A_min=A_min,p_min=p_min))   # full change

    clim = V10 - V00
    inter = V01 - V00
    syn = V11 - V10 - V01 + V00
    return (clim=clim, inter=inter, syn=syn, U0=V00, U1=V11)
end

# (B) ΔBSH, fraction-basis (vs original area), from Z0  ← this is what you asked to “go back to”
"""
decomp_BSH_fraction_fromZ(Z0, pool; keepmask)
Returns mean over consumers of climate-only, interaction-only, synergy, plus B0,B1.
"""
function decomp_BSH_fraction_fromZ(Z0::BitMatrix, pool::SpeciesPool; keepmask::BitVector)
    P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask); P1 = assemble(Z1, pool)
    cons = .!pool.basal
    C0 = size(Z0, 2)  # original cell count; ensures 'vs original area'

    Afrac(Z) = vec(sum(Z; dims=2)) ./ C0

    A0 = Afrac(Z0)[cons]
    A1 = Afrac(Z1)[cons]
    p0 = bsh1_cond_per_species(P0, Z0, pool)[cons]
    p1 = bsh1_cond_per_species(P1, Z1, pool)[cons]
    for x in (p0,p1); @inbounds for i in eachindex(x); if isnan(x[i]); x[i]=0.0; end; end; end

    B0 = A0 .* p0
    B1 = A1 .* p1
    Abar = 0.5 .* (A0 .+ A1)
    pbar = 0.5 .* (p0 .+ p1)
    dA_only  = pbar .* (A1 .- A0)
    dInt_only= Abar .* (p1 .- p0)
    synergy  = (B1 .- B0) .- dA_only .- dInt_only

    (; clim=mean(dA_only), inter=mean(dInt_only), syn=mean(synergy),
       B0=mean(B0), B1=mean(B1))
end

# --- Grid & fields
W,H = 60,60
grid = make_grid(W,H; seed=11)
Fabi = rbf_field(W,H; ncenters=2, scale=80, seed=1)
Fg   = [rbf_field(W,H; ncenters=5, scale=8, seed=s) for s in 2:4]

# --- Pool and assembly
pool, g_basal = build_pool_cf(200; basal_frac=0.35, G=3, diet_k=3, sync_bias=0.85, seed=5)
P0, Aclim0    = assemble_from_fields(pool, Fabi, Fg; g_basal=g_basal,
                                     t_basal=0.60, t_cons=0.35, prey_ratio_min=0.67)

# --- Fit niches + climate mask from what assembled
Z0, μ, b = climate_from_fitted(P0, Aclim0, pool; k=1.2)
pool.mu .= μ; pool.b .= b  # make other routines use the fitted niche

# --- Masks across loss fractions
losses = collect(0.15:0.05:0.80)

# Example: ΔBSH (vs original area) for Random / Clustered
dR = Float64[]; dC = Float64[]; dH = Float64[]
clR = Float64[]; inR = Float64[]; syR = Float64[]
clC = Float64[]; inC = Float64[]; syC = Float64[]

for keep in (1 .- losses)
    kmR = random_mask(grid.C, keep; seed=101)
    kmC = clustered_mask(grid, keep; nseeds=1, seed=202)

    r = decomp_BSH_fraction_fromZ(Z0, pool; keepmask=kmR)
    c = decomp_BSH_fraction_fromZ(Z0, pool; keepmask=kmC)

    push!(clR, r.clim); push!(inR, r.inter); push!(syR, r.syn); push!(dR, r.B1 - r.B0)
    push!(clC, c.clim); push!(inC, c.inter); push!(syC, c.syn); push!(dC, c.B1 - c.B0)
end

function pick_thresholds(A0::Vector{Float64}, p0::Vector{Float64}; qA=0.4, qP=0.4)
    A_min = quantile(A0, qA)
    p_pos = filter(!isnan, p0[p0 .> 0])   # only defined/positive p
    p_min = isempty(p_pos) ? 0.2 : quantile(p_pos, qP)
    return A_min, p_min
end

begin
    fig = Figure(; size=(980,320))
    for (j,(ttl,cl,intr,syn)) in enumerate((("Random",clR,inR,syR),("Clustered",clC,inC,syC)))
        ax = Axis(fig[1,j], title=ttl, xlabel="Area lost (fraction)", ylabel=j==1 ? "ΔBSH decomposition" : "")
        lines!(ax, losses, cl;  label="Climate-only")
        lines!(ax, losses, intr;label="Interaction-only")
        lines!(ax, losses, syn; label="Synergy")
        hlines!(ax, [0], color=(:gray,0.5), linestyle=:dash)
        if j==2; axislegend(ax, position=:lt); end
    end
    display(fig)
end
# ---------- PLOTS ----------

# 3A) Excess curves (clustered vs random)
begin
    fig = Figure(; size=(920,340))

    ax = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Excess ΔU (clustered − random)",
              title="Viability decomposition — community-first, thresholds from data")
    lines!(ax, losses, ex_clim_C, label="Climate-only")
    lines!(ax, losses, ex_int_C,  label="Interaction-only")
    lines!(ax, losses, ex_syn_C,  label="Synergy")
    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax, position=:lt)
    display(fig)
end

begin
    using GeometryBasics
    function stackedbars!(ax, x::AbstractVector{<:Real}, cols::Vector{Symbol}, parts::Matrix{Float64}; width=0.65)
        @assert size(parts,2) == length(x)
        y0 = zeros(length(x))
        for k in 1:size(parts,1)
            yk = view(parts,k,:)
            for i in eachindex(x)
                y = yk[i]
                low  = y >= 0 ? y0[i] : y0[i] + y
                high = y >= 0 ? y0[i] + y : y0[i]
                poly!(ax, Rect(x[i]-width/2, low, width, high-low), color=cols[k])
                y0[i] += y
            end
        end
    end

    # Example: bars at 50% loss comparing Random vs Clustered
    i50 = argmin(abs.(losses .- 0.5))
    parts = [clR[i50] clC[i50];
            inR[i50] inC[i50];
            syR[i50] syC[i50]]
    fig = Figure(; size=(550,380)); ax = Axis(fig[1,1], xlabel="Scenario", ylabel="ΔBSH @ 50%")
    stackedbars!(ax, [1,2], [:dodgerblue,:orange,:forestgreen], parts; width=0.7)
    ax.xticks = ([1,2], ["Random","Clustered"]); hlines!(ax,[0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax, [PolyElement(color=:dodgerblue), PolyElement(color=:orange), PolyElement(color=:forestgreen)],
                ["Climate","Interaction","Synergy"]; position=:lt)
    display(fig)
end

# 3B) One barplot at 50% loss (stacked)
begin
    # same dR, dC as above...
    x = [1, 2]; labels = ["Random","Clustered"]
    parts = [(:dodgerblue,  dR.clim,  dC.clim,  "Climate"),
             (:orange,      dR.inter, dC.inter, "Interaction"),
             (:forestgreen, dR.syn,   dC.syn,   "Synergy")]

    fig = Figure(; size=(720, 600))
    ax  = Axis(fig[1,1], xlabel="Scenario", ylabel="ΔU at 50% loss",
               title="Stacked viability contributions")

    # draw stacks manually using barplot in cumulative fashion
    w = 0.55
    base = zeros(length(x))  # current top per column
    for (col, yR, yC, _) in parts
        ys = [yR, yC]
        # draw the positive parts on top of base
        pos = map(max, ys, zeros(2))
        barplot!(ax, x, base .+ pos; width=w, color=col)
        # update base by adding the *full* contribution (pos + neg), so next layer aligns
        base .+= ys
    end

    ax.xticks = (x, labels)
    ax.xticklabelrotation[] = π/18
    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax,
        [PolyElement(color=:dodgerblue), PolyElement(color=:orange), PolyElement(color=:forestgreen)],
        ["Climate","Interaction","Synergy"]; position=:lt)

    display(fig)
end

# For the same i50
begin
    A1 = climate_area_per_species(apply_mask(Z0, kmC))[cons]
    P1 = assemble(apply_mask(Z0, kmC), pool)
    p1 = bsh1_cond_per_species(P1, apply_mask(Z0, kmC), pool)[cons]
    for x in (p1,); @inbounds for i in eachindex(x); if isnan(x[i]); x[i]=0.0; end; end; end

    fig = Figure(; size=(860,300))
    ax1 = Axis(fig[1,1], title="A0 vs A1 (consumers)", xlabel="A0", ylabel="A1")
    scatter!(ax1, A0, A1); vlines!(ax1, [A_min]); hlines!(ax1, [A_min])

    ax2 = Axis(fig[1,2], title="p0 vs p1 (consumers)", xlabel="p0", ylabel="p1")
    scatter!(ax2, p0, p1); vlines!(ax2, [p_min]); hlines!(ax2, [p_min])
    display(fig)
end
