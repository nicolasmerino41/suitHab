curves = Metrics.abs_bsh_vs_loss(; rng, pool=pool_low, grid=grid_grad, pars=pars,
                                 A_fn=A_fn_punch, loss_fracs=0.2:0.1:0.8, seed_A=1)
figABS = Figs.fig_abs_bsh_vs_loss(curves; title="Absolute BSH vs loss — gradient")
save("Paper/figs/ABS_BSH_vs_loss.png", figABS)

###############################

realA = run_realities_Aonly(; rng, grid=grid_grad,
         loss_fracs=0.2:0.1:0.8, S=S, basal_frac=basal_frac, seed_A=1,
         Npools=12, geoms=(:random,:clustered,:front), τA=pars.τA)

figR = Figs.fig_realities_Aonly(realA; title="A-only response under ABM / MAB / BAM realities")
save("Paper/figs/Realities_Aonly.png", figR)

