# --- loss-aware predictors (no full HL run needed) --------------------------

const PC_SITE = 0.5927  # site percolation threshold (4-neigh square lattice)

"""
compute_lossaware_predictors(grid, S, basal_frac; A_level,B_level,M_level, τA, τocc,
                             T_frac_on, loss_pick, sim_seed=1234, pool_seed=1)
Returns:
  CRIℓ   = mean over consumers of I[ keep*p0 < p_c ]
  D_tail = climate divergence using **basal** suitability (tail-contrast)
  R_mean = mean diet size (consumers)
  LPRI_tail_min = min prey count in low tail (baseline)
"""
function compute_lossaware_predictors(grid::Grid, S::Int, basal_frac::Float64;
        A_level::Symbol, B_level::Symbol, M_level::Symbol,
        τA::Float64, τocc::Float64, T_frac_on::Float64,
        loss_pick::Float64, sim_seed::Int=1234, pool_seed::Int=1)

    rng_pool = MersenneTwister(hash((sim_seed, :pool, pool_seed)))
    pool = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
    A    = abiotic_matrix(pool, grid)
    keep = 1 - loss_pick

    # climate pass Z
    Z = falses(pool.S, grid.C)
    @inbounds for s in 1:pool.S, c in 1:grid.C
        Z[s,c] = A[s,c] ≥ τA
    end

    # p0 for consumers
    cons = .!pool.basal
    p0 = [mean(@view Z[s,:]) for s in 1:pool.S if cons[s]]
    CRIl = mean(keep .* p0 .< PC_SITE)

    # D_tail from BASAL suitability only
    tails = climate_tails(grid)
    dens_basal = vec(mean(Z[pool.basal, :]; dims=1))
    sum_lo  = sum(dens_basal[tails.low]); sum_hi = sum(dens_basal[tails.high])
    D_tail = abs(sum_hi - sum_lo) / max(sum_hi + sum_lo, 1e-12)

    # diet redundancy
    R_mean = mean([length(pool.E[s]) for s in 1:pool.S if !pool.basal[s]])

    # LPRI in low tail (baseline assemble)
    pars = bam_from_axes(; B_level, M_level, τA, τocc)
    bam  = pars.bam
    mp   = MovementParams(; mode=:none, T=8)  # movement off for baseline prey count
    base_keep = trues(grid.C)
    P0,_,_,_ = assemble_BAM(pool, grid, A, base_keep; bam=bam, mp=mp)

    function prey_count_in_tail(s::Int, tailmask::BitVector)
        idx = findall(i -> Z[s,i] && tailmask[i], 1:grid.C)
        if isempty(idx); return 0.0; end
        acc = 0.0
        for i in idx
            pc = 0
            for q in pool.E[s]; pc += (P0[q,i] ? 1 : 0); end
            acc += pc
        end
        return acc / length(idx)
    end
    LPRI_tail_min = minimum([prey_count_in_tail(s, tails.low) for s in 1:pool.S if cons[s]]; init=0.0)

    return (; CRIl, D_tail, R_mean, LPRI_tail_min)
end

"""
predict_geometry_ranking_lossaware(indices; movement_on)
Heuristic:
  if movement_on & CRIl > 0.30  ⇒ worst=random
  elseif D_tail > 0.20 & R_mean ≤ 3.0 & LPRI_tail_min < 1.0 ⇒ worst=front
  else ⇒ worst=clustered
"""
function predict_geometry_ranking_lossaware(ind; movement_on::Bool)
    if movement_on && ind.CRIl > 0.30
        return (; worst=:random, best=:clustered,
                rationale=@sprintf("CRIℓ=%.2f > 0.30 ⇒ fragmentation risk", ind.CRIl))
    elseif (ind.D_tail > 0.20) && (ind.R_mean ≤ 3.0) && (ind.LPRI_tail_min < 1.0)
        return (; worst=:front, best=:random,
                rationale=@sprintf("D_tail=%.2f, low redundancy ⇒ last-prey risk", ind.D_tail))
    else
        return (; worst=:clustered, best=:random,
                rationale="area loss dominates; weak geometry effect")
    end
end

"""
audit_phase_rule_lossaware(...; loss_pick=0.6)
Computes *observed* ΔBSH per geometry at loss_pick and two predictions:
  - baseline prediction (existing rule, via predict_geometry_ranking)
  - loss-aware prediction (CRIl/D_tail/LPRI_tail_min rule)
"""
function audit_phase_rule_lossaware(; nx::Int=60, ny::Int=60, S::Int=150,
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
        obs_order = sort(collect(keys(d)); by=x->d[x])  # most negative = worst

        # baseline indices + prediction
        base = baseline_indicators(grid, S, basal_frac; A_level=axes.A_level, B_level=axes.B_level,
                                   M_level=axes.M_level, τA, τocc, T_frac_on)
        pred_base = predict_geometry_ranking(base.D, base.R, base.Θ; movement_on=(axes.M_level===:on))

        # loss-aware indices + prediction
        ind = compute_lossaware_predictors(grid, S, basal_frac; A_level=axes.A_level, B_level=axes.B_level,
                                           M_level=axes.M_level, τA, τocc, T_frac_on, loss_pick)
        pred_la = predict_geometry_ranking_lossaware(ind; movement_on=(axes.M_level===:on))

        push!(out, (; name=nm, A=axes.A_level, B=axes.B_level, M=axes.M_level,
                    d_random=d[:random], d_clustered=d[:clustered], d_front=d[:front],
                    obs_order=obs_order,
                    base_D=base.D, base_R_mean=base.R.mean, base_Θ=base.Θ,
                    pred_base_worst=pred_base.worst, pred_base_best=pred_base.best,
                    la_CRIl=ind.CRIl, la_D_tail=ind.D_tail, la_R_mean=ind.R_mean, la_LPRI_min=ind.LPRI_tail_min,
                    pred_la_worst=pred_la.worst, pred_la_best=pred_la.best))
    end
    return out
end

"""
plot_audit_phase_rule_enhanced(audit_la; loss_pick=0.6)
Shows baseline vs loss-aware confusion side-by-side + spread bars.
"""
function plot_audit_phase_rule_enhanced(audit_la; loss_pick=0.6)
    geoms = (:random,:clustered,:front)
    glabel = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")
    idx = Dict(:random=>1,:clustered=>2,:front=>3)

    # confusion matrices
    M_base = zeros(Int,3,3); M_la = zeros(Int,3,3)
    for a in audit_la
        ow = a.obs_order[1]
        M_base[idx[ow], idx[a.pred_base_worst]] += 1
        M_la[idx[ow], idx[a.pred_la_worst]]     += 1
    end

    # spread bars
    spread = [maximum((a.d_random,a.d_clustered,a.d_front)) - minimum((a.d_random,a.d_clustered,a.d_front)) for a in audit_la]
    labels = [Symbol("A_"*String(a.A)*"_B_"*String(a.B)*"_M_"*String(a.M)) for a in audit_la]
    ord = sortperm(spread; rev=true)

    fig = Figure(; size=(1400, 600))
    Label(fig[0,:], @sprintf("Phase rule audit at loss=%.2f — baseline vs loss-aware predictors", loss_pick);
          fontsize=18, padding=(0,0,8,0))

    # Baseline confusion
    ax1 = Axis(fig[1,1], title="Observed worst vs Predicted worst (baseline)",
               xticks=(1:3, [glabel[g] for g in geoms]), yticks=(1:3,[glabel[g] for g in geoms]),
               xlabel="Predicted", ylabel="Observed")
    heatmap!(ax1, M_base; colormap=:Blues)
    for r in 1:3, c in 1:3
        text!(ax1, c, r, text=string(M_base[r,c]), align=(:center,:center))
    end

    # Loss-aware confusion
    ax2 = Axis(fig[1,2], title="Observed worst vs Predicted worst (loss-aware)",
               xticks=(1:3, [glabel[g] for g in geoms]), yticks=(1:3,[glabel[g] for g in geoms]),
               xlabel="Predicted", ylabel="Observed")
    heatmap!(ax2, M_la; colormap=:Greens)
    for r in 1:3, c in 1:3
        text!(ax2, c, r, text=string(M_la[r,c]), align=(:center,:center))
    end

    # Geometry spread bar
    ax3 = Axis(
        fig[2,1:2],
        title="Geometry sensitivity (spread = max-min of ΔBSH)", ylabel="Spread",
        xticklabelrotation=π/6
        )
    barplot!(ax3, 1:length(ord), spread[ord])
    ax3.xticks = (1:length(ord), string.(labels[ord]))


    display(fig)
    return fig
end

A = audit_phase_rule_lossaware(; loss_pick=0.6)

plot_audit_phase_rule_enhanced(A)