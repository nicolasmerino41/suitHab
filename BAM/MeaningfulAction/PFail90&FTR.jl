# ---------- Tail exposure by a front (loss-aware) ----------
# We use the fitted climate axis t in [0,1]. A front that removes loss_frac=f
# can be idealized as removing the lowest f (or highest f) of t. We report the
# *maximum* overlap with either climate tail so direction errors don't kill signal.

function _climate_axis_proj(grid::Grid)
    C = grid.C
    X = hcat(ones(C), grid.xy[1,:], grid.xy[2,:])  # [1 x y]
    β = X \ grid.climate
    a, b = β[2], β[3]
    t = a .* grid.xy[1,:] .+ b .* grid.xy[2,:]
    t .-= minimum(t); t ./= (maximum(t) + eps())
    return t
end

"""
FTR = fraction of the climate tail removed by a front at loss f.
Returns (FTR, side) where side ∈ (:low, :high).
"""
function front_tail_exposure(grid::Grid, f::Float64; tail_q::Float64=0.20)
    t = _climate_axis_proj(grid)
    # tail masks
    lo = t .<= quantile(t, tail_q)
    hi = t .>= quantile(t, 1 - tail_q)

    # ideal front: remove the lowest f OR highest f of t
    thr_low  = quantile(t, f)
    thr_high = quantile(t, 1 - f)

    removed_low  = t .<= thr_low      # removing the lower side
    removed_high = t .>= thr_high     # removing the upper side

    F_lo = count(removed_low  .& lo) / max(count(lo), 1)
    F_hi = count(removed_high .& hi) / max(count(hi), 1)

    if F_lo ≥ F_hi
        return (F_lo, :low)
    else
        return (F_hi, :high)
    end
end

# ---------- Tail prey-failure probability (baseline) ----------
# For each consumer s, among climate-suitable tail cells, compute the fraction
# of cells with *zero* prey present (based on baseline assembly with movement off).
# Return the 90th percentile across consumers (captures the "weakest link" effect).

const PF_Q = 0.90

function prey_failure_p90!(rng::AbstractRNG, pool::SpeciesPool, grid::Grid,
                           A::Matrix{Float64}; τA::Float64, B_level::Symbol, τocc::Float64)
    # baseline assembly with movement off (we only care about prey layers)
    base_keep = trues(grid.C)
    bam  = bam_from_axes(; B_level=B_level, M_level=:off, τA=τA, τocc=τocc).bam
    mp   = MovementParams(; mode=:none, T=8)
    P0, _, _, _ = assemble_BAM(pool, grid, A, base_keep; bam=bam, mp=mp)

    # climate pass boolean
    Z = falses(pool.S, grid.C)
    @inbounds for s in 1:pool.S, c in 1:grid.C
        Z[s,c] = A[s,c] ≥ τA
    end

    # climate axis tails
    t = _climate_axis_proj(grid)
    lo = t .<= quantile(t, 0.20)
    hi = t .>= quantile(t, 0.80)

    cons = .!pool.basal
    pfails = Float64[]
    for s in 1:pool.S
        cons[s] || continue
        # choose the threatened tail as the one where s has *more* suitable cells
        n_lo = count(i -> Z[s,i] && lo[i], 1:grid.C)
        n_hi = count(i -> Z[s,i] && hi[i], 1:grid.C)
        tailmask = (n_lo ≥ n_hi) ? lo : hi

        denom = 0
        zeros = 0
        for i in 1:grid.C
            if Z[s,i] && tailmask[i]
                denom += 1
                # count prey presence in that cell
                has_prey = false
                for q in pool.E[s]
                    if P0[q,i]; has_prey = true; break; end
                end
                if !has_prey; zeros += 1; end
            end
        end
        if denom > 0
            push!(pfails, zeros / denom)
        end
    end

    if isempty(pfails)
        return 0.0
    else
        return quantile(pfails, PF_Q)
    end
end

# ---------- Loss-aware predictors (extended) ----------
const PC_SITE = 0.5927

"""
compute_lossaware_predictors_ext(grid, S, basal_frac; A_level,B_level,M_level, τA, τocc, T_frac_on, loss_pick)

Returns:
  CRIl           :: loss-aware connectivity risk (percolation proxy)
  FTR            :: fraction of climate tail removed by front at loss_pick
  Pfail_p90      :: 90th percentile of tail prey-failure across consumers
  side           :: :low or :high (which tail was considered for FTR)
"""
function compute_lossaware_predictors_ext(grid::Grid, S::Int, basal_frac::Float64;
        A_level::Symbol, B_level::Symbol, M_level::Symbol,
        τA::Float64, τocc::Float64, T_frac_on::Float64, loss_pick::Float64,
        sim_seed::Int=1234, pool_seed::Int=1)

    rng = MersenneTwister(hash((sim_seed, :pool, pool_seed)))
    pool = build_pool_from_axes(rng; S, basal_frac, A_level, B_level)
    A    = abiotic_matrix(pool, grid)

    # CRIℓ
    Z = falses(pool.S, grid.C)
    @inbounds for s in 1:pool.S, c in 1:grid.C
        Z[s,c] = A[s,c] ≥ τA
    end
    cons = .!pool.basal
    p0 = [mean(@view Z[s,:]) for s in 1:pool.S if cons[s]]
    keep = 1 - loss_pick
    CRIl = mean(keep .* p0 .< PC_SITE)

    # front tail exposure (loss-aware, map-only)
    FTR, side = front_tail_exposure(grid, loss_pick; tail_q=0.20)

    # tail prey-failure p90 (baseline metaweb, movement off)
    Pfail_p90 = prey_failure_p90!(rng, pool, grid, A; τA=τA, B_level=B_level, τocc=τocc)

    return (; CRIl, FTR, Pfail_p90, side)
end

# ---------- Prediction rule (simple, falsifiable cutoffs) ----------
"""
predict_geometry_ranking_lossaware_ext(ind; movement_on, B_level)
Rule:
  if movement_on && ind.CRIl > τ_M         ⇒ worst=random
  elseif (B_level != :none) &&
         (ind.FTR > τ_FTR) && (ind.Pfail_p90 > τ_P)  ⇒ worst=front
  else ⇒ worst=clustered
"""

function predict_geometry_ranking_lossaware_ext(ind; movement_on::Bool, B_level::Symbol,
        τ_M::Float64=0.30, τ_FTR::Float64=0.35, τ_P::Float64=0.50)
    if movement_on && ind.CRIl > τ_M
        return (; worst=:random, best=:clustered,
                rationale=@sprintf("CRIℓ=%.2f>%.2f ⇒ fragmentation risk", ind.CRIl, τ_M))
    elseif (B_level != :none) && (ind.FTR > τ_FTR) && (ind.Pfail_p90 > τ_P)
        return (; worst=:front, best=:random,
                rationale=@sprintf("FTR=%.2f & Pfail,90%%=%.2f ⇒ last-prey risk", ind.FTR, ind.Pfail_p90))
    else
        return (; worst=:clustered, best=:random,
                rationale="area-dominated; weak geometry effect")
    end
end

# ---------- Wire into your audit ----------
function audit_phase_rule_lossaware_ext(; nx::Int=60, ny::Int=60, S::Int=150,
        basal_frac::Float64=0.45, loss_pick::Float64=0.6,
        seeds_pool=1:5, seeds_mask=1:5, sim_seed::Int=1234,
        τA::Float64=0.5, τocc::Float64=0.35, T_frac_on::Float64=0.02)

    grid = make_grid(nx, ny; seed=42)
    out = NamedTuple[]

    for (nm, axes) in COMBOS
        # observed ΔBSH at loss_pick
        d = Dict{Symbol,Float64}()
        for g in (:random, :clustered, :front)
            rr = sweep_dBSH_axes(; grid,S,basal_frac, A_level=axes.A_level, B_level=axes.B_level, M_level=axes.M_level,
                                 hl_kind=g, loss_fracs=[loss_pick], seeds_pool, seeds_mask, sim_seed)
            d[g] = rr.y[1]
        end
        obs_order = sort(collect(keys(d)); by=x->d[x]) # most negative = worst

        # loss-aware indicators + prediction
        ind = compute_lossaware_predictors_ext(grid, S, basal_frac; A_level=axes.A_level, B_level=axes.B_level,
                                               M_level=axes.M_level, τA, τocc, T_frac_on, loss_pick)
        pred = predict_geometry_ranking_lossaware_ext(ind; movement_on=(axes.M_level===:on), B_level=axes.B_level)

        push!(out, (; name=nm, A=axes.A_level, B=axes.B_level, M=axes.M_level,
                    d_random=d[:random], d_clustered=d[:clustered], d_front=d[:front],
                    obs_order=obs_order,
                    CRIl=ind.CRIl, FTR=ind.FTR, Pfail_p90=ind.Pfail_p90, tail_side=ind.side,
                    pred_worst=pred.worst, pred_best=pred.best, rationale=pred.rationale))
    end
    return out
end

"""
plot_audit_phase_rule_lossaware_ext(audit; loss_pick=0.6, τ_M=0.30, τ_FTR=0.35, τ_P=0.50)

Parameters
- audit: Vector of NamedTuples from `audit_phase_rule_lossaware_ext`.
- loss_pick, τ_M, τ_FTR, τ_P: shown on the figure; used to draw threshold lines.

Panels
1) Confusion matrix: observed worst geometry vs predicted worst (extended rule).
2) Indicator map for front risk: x = FTR (tail fraction removed), y = Pfail_p90 (tail prey-failure p90),
   points colored by observed worst geometry; vertical/horizontal lines at τ_FTR and τ_P.
3) Geometry sensitivity bars: spread = max(ΔBSH) − min(ΔBSH), sorted by spread.
"""
function plot_audit_phase_rule_lossaware_ext(audit;
        loss_pick::Float64=0.6, τ_M::Float64=0.30, τ_FTR::Float64=0.35, τ_P::Float64=0.50)

    geoms   = (:random, :clustered, :front)
    glabel  = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")
    gindex  = Dict(:random=>1, :clustered=>2, :front=>3)
    gcolor  = Dict(:random=>RGBf(0.25,0.45,0.85), :clustered=>RGBf(0.95,0.60,0.15), :front=>RGBf(0.20,0.65,0.35))

    # ---- crunch audit arrays
    obs_worst   = [a.obs_order[1]  for a in audit]
    pred_worst  = [a.pred_worst    for a in audit]
    FTR         = [a.FTR           for a in audit]
    Pfail_p90   = [a.Pfail_p90     for a in audit]
    CRIl        = [a.CRIl          for a in audit]
    names       = [Symbol("A_"*String(a.A)*"_B_"*String(a.B)*"_M_"*String(a.M)) for a in audit]
    spread      = [maximum((a.d_random,a.d_clustered,a.d_front)) -
                   minimum((a.d_random,a.d_clustered,a.d_front)) for a in audit]

    # confusion matrix
    M = zeros(Int, 3, 3)  # rows = observed, cols = predicted
    for i in eachindex(audit)
        M[gindex[obs_worst[i]], gindex[pred_worst[i]]] += 1
    end

    # ---- figure
    fig = Figure(; size=(1400, 650))
    Label(fig[0,:],
          @sprintf("Phase rule audit (loss-aware, extended) — loss=%.2f  |  τ_M=%.2f, τ_FTR=%.2f, τ_P=%.2f",
                   loss_pick, τ_M, τ_FTR, τ_P);
          fontsize=18, padding=(0,0,8,0))

    # (1) Confusion
    ax1 = Axis(fig[1,1], title="Observed worst vs Predicted worst",
               xlabel="Predicted", ylabel="Observed",
               xticks=(1:3, [glabel[g] for g in geoms]),
               yticks=(1:3, [glabel[g] for g in geoms]))
    heatmap!(ax1, M; colormap=:Blues)
    for r in 1:3, c in 1:3
        text!(ax1, c, r, text=string(M[r,c]), align=(:center,:center), color=:black)
    end

    # (2) Indicator map for front risk
    ax2 = Axis(fig[1,2], title="Front-risk indicators",
               xlabel="FTR (fraction of climate tail removed)",
               ylabel="Pfail,90% (tail prey-failure probability)")
    # thresholds
    vlines!(ax2, [τ_FTR]; color=RGBAf(0,0,0,0.4), linestyle=:dash)
    hlines!(ax2, [τ_P];    color=RGBAf(0,0,0,0.4), linestyle=:dash)
    # scatter by observed worst geometry
    for g in geoms
        idx = findall(i -> obs_worst[i] === g, eachindex(audit))
        isempty(idx) && continue
        scatter!(ax2, FTR[idx], Pfail_p90[idx];
                 color=gcolor[g], markersize=10, label=glabel[g])
    end
    axislegend(ax2; position=:cb, framevisible=false, labelsize=10, padding=(2,2,2,2))

    # annotate CRIℓ info (as text; we keep the panel uncluttered)
    meanCRIl = mean(CRIl)
    text!(ax2, 0.02, 0.98, text=@sprintf("Mean CRIℓ across combos: %.2f (τ_M=%.2f)", meanCRIl, τ_M),
          align=(:left,:top), space=:relative, fontsize=10, color=RGBAf(0,0,0,0.7))

    # (3) Geometry sensitivity bars
    ord = sortperm(spread; rev=true)
    ax3 = Axis(fig[2,1:2], title="Geometry sensitivity (spread = max−min of ΔBSH)",
               ylabel="Spread")
    barplot!(ax3, 1:length(ord), spread[ord])
    ax3.xticks = (1:length(ord), string.(names[ord]))
    ax3.xticklabelrotation = π/6

    display(fig)
    return fig
end

A = audit_phase_rule_lossaware_ext(; loss_pick=0.4)
plot_audit_phase_rule_lossaware_ext(A; loss_pick=0.4)