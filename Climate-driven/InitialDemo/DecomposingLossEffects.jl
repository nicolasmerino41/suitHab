############################
# A. COMPONENTS ON AREA BASIS
############################
# A_s = suitable / C_full
# p_s = (suitable & has-prey) / suitable  (0 if suitable==0)
# B_s = A_s * p_s
function components_area(P::BitMatrix, Z::BitMatrix, pool::SpeciesPool, Cfull::Int)
    S, C = size(P); @assert C <= Cfull
    A = zeros(Float64, S)
    p = zeros(Float64, S)
    B = zeros(Float64, S)

    @inbounds for s in 1:S
        suitable = count(@view Z[s, :])                # number of climate-suitable cells (after mask)
        A[s] = suitable / Cfull

        if pool.basal[s]
            p[s] = 1.0
            B[s] = A[s]
        else
            if suitable == 0
                p[s] = 0.0
                B[s] = 0.0
            else
                ok = 0
                prey = pool.E[s]
                # in each climate-suitable cell, check if ≥1 prey is present
                for c in 1:size(Z,2)
                    if Z[s,c] == 1
                        for q in prey
                            if P[q,c] == 1
                                ok += 1
                                break
                            end
                        end
                    end
                end
                p[s] = ok / suitable
                B[s] = A[s] * p[s]
            end
        end
    end
    return A, p, B
end

############################
# B. MIDPOINT DECOMPOSITION
############################
struct Decomp{T}
    dB::T
    dA_only::T
    dInt_only::T
    synergy::T
end

# single-species (scalars)
@inline function decompose_midpoint(A0, p0, B0, A1, p1, B1)
    dA = A1 - A0
    dp = p1 - p0
    Abar = 0.5*(A0 + A1)
    pbar = 0.5*(p0 + p1)
    dA_only   = pbar * dA
    dInt_only = Abar * dp
    dB        = B1 - B0
    synergy   = dB - dA_only - dInt_only
    return Decomp(dB, dA_only, dInt_only, synergy)
end

############################
# C. ONE-SHOT DECOMP AT A MASK
############################
"""
decomp_at_mask(pool, grid; τ, keepmask, consumers_only=true)

Returns:
  - species-level vectors: dB, dA_only, dInt_only, synergy
  - simple aggregates (means over selected species)
Everything on **area basis** (original Cfull denominator).
"""
function decomp_at_mask(pool::SpeciesPool, grid::Grid; τ::Float64=0.5, keepmask::BitVector,
                        consumers_only::Bool=true)
    Cfull = grid.C

    # BEFORE
    Z0 = climate_pass(pool, grid; τ=τ)
    P0 = assemble(Z0, pool)
    A0, p0, B0 = components_area(P0, Z0, pool, Cfull)

    # AFTER (apply mask to climate, then assemble)
    Z1 = apply_mask(Z0, keepmask)
    P1 = assemble(Z1, pool)
    A1, p1, B1 = components_area(P1, Z1, pool, Cfull)

    # choose species set
    sel = consumers_only ? .!pool.basal : trues(pool.S)

    dB   = similar(B0);    dB  .= NaN
    dAcl = similar(B0);    dAcl .= NaN
    dInt = similar(B0);    dInt .= NaN
    dSyn = similar(B0);    dSyn .= NaN

    @inbounds for s in 1:pool.S
        sel[s] || continue
        d = decompose_midpoint(A0[s], p0[s], B0[s], A1[s], p1[s], B1[s])
        dB[s]   = d.dB
        dAcl[s] = d.dA_only
        dInt[s] = d.dInt_only
        dSyn[s] = d.synergy
    end

    # aggregates (drop NaNs)
    μ(x) = mean(filter(!isnan, x))
    agg = ( dB   = μ(dB),
            dAcl = μ(dAcl),
            dInt = μ(dInt),
            dSyn = μ(dSyn) )

    return (dB=dB, dAcl=dAcl, dInt=dInt, dSyn=dSyn, agg=agg,
            A0=A0, p0=p0, B0=B0, A1=A1, p1=p1, B1=B1)
end

############################
# D. CONVENIENCE: RANDOM vs CLUSTERED AT ONE LOSS
############################
"""
compare_decomp_once(grid; keep_frac, τ, nseeds_cluster=1, pool_kwargs...)

Builds a pool with `build_pool(S; pool_kwargs...)`, computes area-basis
decomposition at the same loss for `:random` and `:clustered`, and plots bars.
"""
function compare_decomp_once(; grid::Grid, S::Int, basal_frac::Float64,
        keep_frac::Float64=0.6, τ::Float64=0.5, nseeds_cluster::Int=1, seed_pool::Int=1,
        consumers_only::Bool=true, pool_kwargs...)

    pool = build_pool(S; basal_frac=basal_frac, seed=seed_pool, pool_kwargs...)

    km_r = random_mask(grid.C, keep_frac; seed=123)
    km_c = clustered_mask(grid, keep_frac; nseeds=nseeds_cluster, seed=456)

    dr = decomp_at_mask(pool, grid; τ=τ, keepmask=km_r, consumers_only=consumers_only)
    dc = decomp_at_mask(pool, grid; τ=τ, keepmask=km_c, consumers_only=consumers_only)

    # --- Plot
    begin
        fig = Figure(size=(820,320))
        ax  = Axis(fig[1,1],
            title = "ΔBSH decomposition at loss = $(round(1-keep_frac, digits=2)) (area basis)",
            ylabel= "Contribution",
            xticks = (1:3, ["Climate-only","Interaction-only","Synergy"])
        )
        x   = 1:3
        off = 0.18
        barplot!(ax, x .- off, [dr.agg.dAcl, dr.agg.dInt, dr.agg.dSyn];
                 width=0.35, color=:dodgerblue, label="Random")
        barplot!(ax, x .+ off, [dc.agg.dAcl, dc.agg.dInt, dc.agg.dSyn];
                 width=0.35, color=:orange,    label="Clustered")
        hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
        axislegend(ax, position=:rb)
        display(fig)
    end

    # quick consistency check (should be ~0)
    err_r = dr.agg.dB - (dr.agg.dAcl + dr.agg.dInt + dr.agg.dSyn)
    err_c = dc.agg.dB - (dc.agg.dAcl + dc.agg.dInt + dc.agg.dSyn)
    println("Check: random ΔB - parts = ", err_r, "   clustered = ", err_c)

    return (; random=dr, clustered=dc)
end

# Example parameters (tweak freely later)
res = compare_decomp_once(
    ; grid=grid, S=200, basal_frac=0.4,
    keep_frac=0.6, τ=0.55, nseeds_cluster=1, seed_pool=3,
    consumers_only=false,
    # diet / redundancy
    sigma=0.30, density=0.20, pmax=0.85, R0_mean=12.0, R0_sd=0.4,
    # climate / synchrony
    niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.04,
    b0_basal=0.09, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.05
)

# species-level ΔA, Δp, synergy for a given mask
function deltas_by_species(pool, grid; τ, keepmask, consumers_only=true)
    Cfull = grid.C
    Z0 = climate_pass(pool, grid; τ=τ);  P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask);       P1 = assemble(Z1, pool)
    A0,p0,B0 = components_area(P0, Z0, pool, Cfull)
    A1,p1,B1 = components_area(P1, Z1, pool, Cfull)
    sel = consumers_only ? .!pool.basal : trues(pool.S)
    ΔA = A1 .- A0; Δp = p1 .- p0; syn = (A1.-A0) .* (p1.-p0)
    return (ΔA=ΔA[sel], Δp=Δp[sel], syn=syn[sel])
end

# Example scatter for the same pool+loss used in your figure
pool = build_pool(200; basal_frac=0.35,
    sigma=0.22, density=0.12, pmax=0.70, R0_mean=10.0, R0_sd=0.25,
    niche_mode=:bimodal, mu_basal_centers=(0.25,0.75), mu_basal_sd=0.04,
    b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04)

keep = 0.6; τ = 0.55
km_r = random_mask(grid.C, keep; seed=1)
km_c = clustered_mask(grid, keep; nseeds=1, seed=2)

dr = deltas_by_species(pool, grid; τ=τ, keepmask=km_r)
dc = deltas_by_species(pool, grid; τ=τ, keepmask=km_c)

begin
    fig = Figure(; size=(800,360))
    ax1 = Axis(fig[1,1], title="Random",  xlabel="ΔA", ylabel="Δp")
    ax2 = Axis(fig[1,2], title="Clustered", xlabel="ΔA", ylabel="Δp")
    scatter!(ax1, dr.ΔA, dr.Δp; markersize=4, color=:dodgerblue)
    scatter!(ax2, dc.ΔA, dc.Δp; markersize=4, color=:orange)
    display(fig)
end

function decomp_vs_loss(grid, pool; τ=0.55, kind=:random, loss_fracs=0.0:0.05:0.8, nseeds=1)
    dAcl = Float64[]; dInt = Float64[]; dSyn = Float64[]
    for f in loss_fracs
        keep = 1.0 - f
        km = kind==:random ? random_mask(grid.C, keep; seed=1) :
                             clustered_mask(grid, keep; nseeds=nseeds, seed=2)
        r = decomp_at_mask(pool, grid; τ=τ, keepmask=km)
        push!(dAcl, r.agg.dAcl); push!(dInt, r.agg.dInt); push!(dSyn, r.agg.dSyn)
    end
    return (loss=collect(loss_fracs), dAcl=dAcl, dInt=dInt, dSyn=dSyn)
end

decomp_vs_loss(grid, pool; τ=0.55, kind=:random, loss_fracs=0.0:0.05:0.8, nseeds=1)

begin
    res = decomp_vs_loss(grid, pool; τ=0.55, kind=:random, loss_fracs=0.0:0.05:0.8, nseeds=1)

    fig = Figure(size=(800,400))
    ax = Axis(fig[1,1],
        xlabel="Area lost (fraction)",
        ylabel="ΔBSH decomposition"
    )

    lines!(ax, res.loss, res.dAcl; label="Climate-only", color=:dodgerblue)
    lines!(ax, res.loss, res.dInt; label="Interaction-only", color=:darkorange)
    lines!(ax, res.loss, res.dSyn; label="Synergy", color=:seagreen)

    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax; position=:rb)

    display(fig)
end
