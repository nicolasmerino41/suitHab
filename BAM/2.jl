"""
Compute (ΔF_s, ΔBSH_s) per species between base_keep and keep, for movement=:component.
Returns a NamedTuple with vectors (dF, dBSH) and a mask of consumers.
"""
function deltas_F_vs_BSH(pool::SpeciesPool, grid::Grid, A::Matrix{Float64},
                         base_keep::BitVector, keep::BitVector;
                         τA=0.5, τocc=0.35, α=1.0, β=1.0, μ=0.6, γ=3.0, T=8)

    bam = BAMParams(α=α, β=β, μ=μ, γ=γ, τA=τA, τocc=τocc)
    mp  = MovementParams(mode=:component, T=T)

    # baseline
    P0, B0, _, M0 = assemble_BAM(pool, grid, A, base_keep; bam=bam, mp=mp)
    S0 = species_stats(pool, grid, A, base_keep, P0, B0, M0, bam)

    # scenario
    P1, B1, _, M1 = assemble_BAM(pool, grid, A, keep; bam=bam, mp=mp)
    S1 = species_stats(pool, grid, A, keep, P1, B1, M1, bam)

    dF   = S1.F   .- S0.F
    dBSH = S1.BSH .- S0.BSH
    cons = consumer_mask(pool)
    (; dF, dBSH, cons)
end

"""
Make a 1×3 Figure: Random / Clustered / Front.
Each panel: ΔBSH (y) vs ΔF (x) for consumers, with identity line, r and slope.
Reuses your RNG-safe builders from earlier.
"""
function plot_alignment_three(; grid::Grid, S::Int=120, basal_frac=0.45,
                              τA=0.5, τocc=0.35, α=1.0, β=1.5, μ=0.6, γ=3.0, T=8,
                              keep_frac=0.5,
                              nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
                              sim_seed::Int=2025, pool_seed::Int=1, mask_seed::Int=1)

    # one pool & A to isolate alignment (you can loop seeds if you want CIs)
    rng_pool = MersenneTwister(hash((sim_seed, :pool, pool_seed)))
    pool     = build_pool(S; rng=rng_pool, basal_frac=basal_frac)
    A        = abiotic_matrix(pool, grid)
    base_keep = trues(grid.C)

    hl_list = (:random, :clustered, :front)
    labels  = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")

    fig = Figure(; size=(1260, 360))

    for (col, hk) in enumerate(hl_list)
        rng_mask = MersenneTwister(hash((sim_seed, :mask, mask_seed, hk)))
        keep = hk === :random    ? random_mask(rng_mask, grid.C, keep_frac) :
               hk === :clustered ? clustered_mask(rng_mask, grid, keep_frac; nseeds=nseeds_cluster) :
               hk === :front     ? frontlike_mask(rng_mask, grid, keep_frac; axis=front_axis, noise=front_noise) :
               error("Unknown HL $(hk)")

        d = deltas_F_vs_BSH(pool, grid, A, base_keep, keep;
                            τA=τA, τocc=τocc, α=α, β=β, μ=μ, γ=γ, T=T)
        x = d.dF[d.cons]; y = d.dBSH[d.cons]

        ax = Axis(fig[1,col], title=labels[hk], xlabel="ΔF per species", ylabel=(col==1 ? "ΔBSH per species" : ""))
        scatter!(ax, x, y; markersize=5)
        # identity line
        lo = minimum(vcat(x,y)); hi = maximum(vcat(x,y))
        lines!(ax, [lo, hi], [lo, hi]; linestyle=:dash)
        # stats
        r = isfinite(std(x)) && isfinite(std(y)) && std(x)>0 && std(y)>0 ? cor(x,y) : NaN
        slope = (x -> sum((x .- mean(x)).*(y .- mean(y)))/sum((x .- mean(x)).^2))(x)
        text!(ax, 0.02, 0.95, text= @sprintf("r = %.2f\nslope = %.2f", r, slope))
    end

    display(fig)
    return fig
end

grid = make_grid(60, 50; seed=42)

# Fixed movement choice for the paper:
#   movement = :component, μ=0.6, T=8
plot_alignment_three(; grid, β=1.5, μ=0.6, T=8, keep_frac=0.5)
