include("BAMRunner.jl")
using .BAMRunner

# List all combo names
keys(BAMRunner.COMBOS)  # 18 names like :A_divergent__B_strong__M_on

# Run one combo (produces the 2×3 figure for Random/Clustered/Front)
_ = BAMRunner.run_combo_name(:A_divergent__B_strong__M_on; S=150, loss_fracs=0.2:0.1:0.8,
                             align_keep_frac=0.5, T=8, seeds_pool=1:6, seeds_mask=1:6)
