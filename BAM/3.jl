"""
run_combo_with_alignment(combo; ...) -> (res_by_hl, fig)

Top row: ΔF component lines for Random / Clustered / Front (mean over consumers).
Bottom row: alignment scatter ΔBSH vs ΔF (consumers) at a chosen keep_frac.
Works with any movement choice in the combo (e.g., :component or :access).
"""
function run_combo_with_alignment(combo::Symbol;
        nx::Int=60, ny::Int=50, S::Int=120, basal_frac::Float64=0.45,
        loss_fracs = 0.2:0.1:0.8,
        seed_grid::Int=42, seeds_pool = 1:5, seeds_mask = 1:5,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        sim_seed::Int=1234,
        # alignment panel options
        align_keep_frac::Float64=0.5, align_pool_seed::Int=1, align_mask_seed::Int=1)

    @assert haskey(BAM_COMBOS, combo) "Unknown combo $(combo). Keys: $(keys(BAM_COMBOS))"
    p = BAM_COMBOS[combo]

    grid = make_grid(nx, ny; seed=seed_grid)

    # --- compute top-row sweeps (threaded) ---
    res = Dict{Symbol,Any}()
    for hk in (:random, :clustered, :front)
        res[hk] = run_sweep_threaded(; grid, S, basal_frac,
            τA=0.5, τocc=0.35,
            α=p.α, β=p.β, μ=p.μ, γ=p.γ,
            move_mode=p.move_mode, λ=p.λ, T=p.T,
            hl_kind=hk, nseeds_cluster=nseeds_cluster, front_axis=front_axis, front_noise=front_noise,
            loss_fracs=loss_fracs,
            seeds_pool=seeds_pool, seeds_mask=seeds_mask,
            sim_seed=sim_seed)
    end

    # --- compute bottom-row alignment (single pool & mask for clarity) ---
    rng_pool = MersenneTwister(hash((sim_seed, :pool, align_pool_seed)))
    pool     = build_pool(S; rng=rng_pool, basal_frac=basal_frac)
    A        = abiotic_matrix(pool, grid)
    base_keep = trues(grid.C)

    function one_alignment(hk::Symbol)
        rng_mask = MersenneTwister(hash((sim_seed, :mask, align_mask_seed, hk)))
        keep = hk === :random    ? random_mask(rng_mask, grid.C, 1 - align_keep_frac) :
               hk === :clustered ? clustered_mask(rng_mask, grid, 1 - align_keep_frac; nseeds=nseeds_cluster) :
               hk === :front     ? frontlike_mask(rng_mask, grid, 1 - align_keep_frac; axis=front_axis, noise=front_noise) :
               error("Unknown HL $(hk)")
        d = deltas_F_vs_BSH(pool, grid, A, base_keep, keep;
                            τA=0.5, τocc=0.35, α=p.α, β=p.β, μ=p.μ, γ=p.γ, T=p.T)
        x = d.dF[d.cons]; y = d.dBSH[d.cons]
        r = (std(x)>0 && std(y)>0) ? cor(x,y) : NaN
        slope = sum((x .- mean(x)).*(y .- mean(y))) / sum((x .- mean(x)).^2)
        (; x, y, r, slope)
    end

    A_rand = one_alignment(:random)
    A_cl   = one_alignment(:clustered)
    A_fr   = one_alignment(:front)

    # --- make the 2×3 figure ---
    fig = Figure(; size=(1260, 760))
    labels = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")

    # top row: ΔF components
    for (col, hk) in enumerate((:random, :clustered, :front))
        r = res[hk]
        ax = Axis(fig[1,col], xlabel="Area lost (fraction)", ylabel="ΔF (mean cons.)",
                  title=labels[hk])
        lines!(ax, r.x, r.dF,  label="Total ΔF")
        lines!(ax, r.x, r.dA,  label="Abiotic")
        lines!(ax, r.x, r.dB,  label="Biotic")
        lines!(ax, r.x, r.dM,  label="Movement")
        lines!(ax, r.x, r.dSyn,label="Synergy")
        if col == 1; axislegend(ax; position=:lt, framevisible=false); end
    end

    # bottom row: alignment
    for (col, (hk, Ahl)) in enumerate(zip((:random,:clustered,:front), (A_rand, A_cl, A_fr)))
        ax = Axis(fig[2,col], xlabel="ΔF per species", ylabel=(col==1 ? "ΔBSH per species" : ""),
                  title = labels[hk]*@sprintf("  (r=%.2f, slope=%.2f)", Ahl.r, Ahl.slope))
        scatter!(ax, Ahl.x, Ahl.y; markersize=5)
        lo = minimum(vcat(Ahl.x, Ahl.y)); hi = maximum(vcat(Ahl.x, Ahl.y))
        lines!(ax, [lo, hi], [lo, hi]; linestyle=:dash)
    end

    Label(fig[0, :], p.title; fontsize=16, padding=(0,0,10,0))
    display(fig)
    return (res_by_hl=res, fig=fig)
end

_ = run_combo_with_alignment(
    :threshold_component;
    nx=60, ny=60, S=150,
    loss_fracs=0.2:0.1:0.8,
    align_keep_frac=0.5
)

_ = run_combo_with_alignment(:baseline)
_ = run_combo_with_alignment(:strongB_access)

