include("BAMRunnerHelpers.jl")

# List all combo names
keys(COMBOS)  # 18 names like :A_divergent__B_strong__M_on

# Run one combo (produces the 2×3 figure for Random/Clustered/Front)
_ = run_combo_name(
    :A_intermediate__B_strong__M_on; S=150, loss_fracs=0.2:0.1:0.8,
    align_keep_frac=0.5, T=8, seeds_pool=1:6, seeds_mask=1:6
)

for combo in keys(COMBOS)
    _ = run_combo_name(
        combo; S=150, loss_fracs=0.2:0.1:0.8,
        align_keep_frac=0.5, T=8, seeds_pool=1:6, seeds_mask=1:6
    )
end

# 1) Compute summaries at a fixed loss (say 0.6)
S = summarize_all_combos(; nx=60, ny=60, S=150, loss_pick=0.6, T_frac_on=0.05)

# 2) Plot one figure with ternary shares (top) + geometry sensitivity bars (bottom)
_ = plot_ternary_summary(S; title="Shares & Geometry sensitivity at 60% loss")

#######################################################
#######################################################
# 1) Compute elasticities at, say, 60% loss (threaded across all 18 combos)
E = elasticity_summarize_all(
    ; nx=60, ny=60, S=150, loss_pick=0.6,
    seeds_pool=1:6, seeds_mask=1:6,
    T_frac_on=0.02
)

# 2) One figure with ternary *elasticity shares* (top) and an
#    "elasticity geometry sensitivity" bar chart (bottom).
_ = plot_elasticity_summary(E; title="Elasticity shares at 60% loss")
#######################################################
#######################################################
# Compute both summaries at the same loss (e.g., 0.6)
raw = summarize_all_combos(; nx=60, ny=60, S=150, loss_pick=0.8,
                                     seeds_pool=1:6, seeds_mask=1:6)
elas = elasticity_summarize_all(; nx=60, ny=60, S=150, loss_pick=0.8,
                                          seeds_pool=1:6, seeds_mask=1:6)

# Side-by-side comparison
_ = compare_raw_vs_elasticity(
    raw; sum_elast=elas,
    title="Raw ΔF vs Elasticity shares (all combos)"
)