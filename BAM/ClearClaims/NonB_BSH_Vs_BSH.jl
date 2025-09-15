"""
sweep_AM_vs_ABM(; grid, S, basal_frac, A_level, B_level, M_level,
                 hl_kind, loss_fracs, seeds_pool, seeds_mask, sim_seed,
                 nseeds_cluster, front_axis, front_noise, τA, τocc, T_frac_on)

Returns: (x, mean_AM, mean_ABM, gap_biotic)
  x           :: loss fractions
  mean_AM     :: mean consumer BSH with Biotic OFF (A×M; B ignored)
  mean_ABM    :: mean consumer BSH with Biotic ON  (A×B×M)
  gap_biotic  :: mean_AM - mean_ABM  (underestimation if B is ignored)
"""
function sweep_AM_vs_ABM(; grid::Grid, S::Int, basal_frac::Float64,
        A_level::Symbol, B_level::Symbol, M_level::Symbol,
        hl_kind::Symbol, loss_fracs,
        seeds_pool, seeds_mask, sim_seed::Int,
        nseeds_cluster::Int, front_axis::Symbol, front_noise::Float64,
        τA::Float64, τocc::Float64, T_frac_on::Float64)

    # local helper: make a keepmask with the right RNG-first signature
    make_keepmask = function(rng::AbstractRNG, kind::Symbol, keep_frac::Float64)
        if kind === :random
            return random_mask(rng, grid.C, keep_frac)
        elseif kind === :clustered
            return clustered_mask(rng, grid, keep_frac; nseeds=nseeds_cluster)
        elseif kind === :front
            return frontlike_mask(rng, grid, keep_frac; axis=front_axis, noise=front_noise)
        else
            error("Unknown hl_kind=$kind")
        end
    end

    xs   = Float64[]
    AMs  = Float64[]
    ABMs = Float64[]
    GAPs = Float64[]

    for f in loss_fracs
        keepfrac = 1 - f
        vals_AM  = Float64[]
        vals_ABM = Float64[]

        for ps in seeds_pool, ms in seeds_mask
            # pool + abiotic
            rng_pool = MersenneTwister(hash((sim_seed, :pool, ps)))
            pool = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
            A    = abiotic_matrix(pool, grid)

            # movement/biotic knobs (your API)
            pars = bam_from_axes(; B_level, M_level, τA, τocc)
            bam, mp = pars.bam, pars.mp
            # scale T by T_frac_on if you do that elsewhere; here we keep mp.T as provided.

            # mask with RNG up front
            rng_mask = MersenneTwister(hash((sim_seed, :mask, ms, hl_kind, f)))
            keep = make_keepmask(rng_mask, hl_kind, keepfrac)

            # --- FULL BSH (A×B×M)
            Pfull, Bsup_full, _, _ = assemble_BAM(pool, grid, A, keep; bam=bam, mp=mp)
            bsh_full = bsh1_per_species(Pfull, (A .>= τA), pool)
            mean_ABM_f = mean(bsh_full[consumer_mask(pool)])

            # --- CLIMATE-ONLY (A×M): Biotic OFF but same movement rule
            pars_AM = bam_from_axes(; B_level=:none, M_level=M_level, τA=τA, τocc=τocc)
            P_am, _, _, _ = assemble_BAM(pool, grid, A, keep; bam=pars_AM.bam, mp=pars_AM.mp)
            bsh_AM = bsh1_per_species(P_am, (A .>= τA), pool)
            mean_AM_f = mean(bsh_AM[consumer_mask(pool)])

            push!(vals_ABM, mean_ABM_f)
            push!(vals_AM,  mean_AM_f)
        end

        push!(xs,  f)
        μAM  = mean(vals_AM)
        μABM = mean(vals_ABM)
        push!(AMs,  μAM)
        push!(ABMs, μABM)
        push!(GAPs, μAM - μABM)
    end

    return (; x=xs, mean_AM=AMs, mean_ABM=ABMs, gap_biotic=GAPs)
end

"""
plot_biotic_gap_dashboard(; grid, S, basal_frac, A_level, B_level, M_level,
                           loss_fracs, seeds_pool, seeds_mask, sim_seed,
                           nseeds_cluster, front_axis, front_noise, τA, τocc, T_frac_on)

Top row: A×M vs A×B×M (Random / Clustered / Front-like)
Bottom row: biotic penalty (A×M − A×B×M) for the three geometries.
"""
function plot_biotic_gap_dashboard(; grid::Grid, S::Int=150, basal_frac::Float64=0.45,
        A_level::Symbol=:intermediate, B_level::Symbol=:soft, M_level::Symbol=:on,
        loss_fracs=0.2:0.1:0.8, seeds_pool=1:6, seeds_mask=1:6, sim_seed::Int=1234,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        τA::Float64=0.5, τocc::Float64=0.35, T_frac_on::Float64=0.02)

    # sweep per geometry with correct mask calls
    R = sweep_AM_vs_ABM(; grid,S,basal_frac,A_level,B_level,M_level,
                        hl_kind=:random,   loss_fracs,seeds_pool,seeds_mask,sim_seed,
                        nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)

    C = sweep_AM_vs_ABM(; grid,S,basal_frac,A_level,B_level,M_level,
                        hl_kind=:clustered,loss_fracs,seeds_pool,seeds_mask,sim_seed,
                        nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)

    F = sweep_AM_vs_ABM(; grid,S,basal_frac,A_level,B_level,M_level,
                        hl_kind=:front,    loss_fracs,seeds_pool,seeds_mask,sim_seed,
                        nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)

    fig = Figure(; size=(1400, 900))
    Label(fig[0,:],
          "Climate-only vs full BSH — A=$(A_level), B=$(B_level), M=$(M_level)";
          fontsize=18, padding=(0,0,8,0))

    # reusable panel
    function panel!(gl, S; title="")
        ax = Axis(gl; xlabel="Area lost (fraction)",
                     ylabel="Mean BSH (consumers)",
                     title=title)
        h1 = lines!(ax, S.x, S.mean_AM;  color=RGBf(0.30,0.45,0.85), linewidth=2,
                    label="Climate-only (A×M)")
        h2 = lines!(ax, S.x, S.mean_ABM; color=RGBf(0.20,0.65,0.35), linewidth=2,
                    label="Full BSH (A×B×M)")
        axislegend(ax, [h1,h2],
                   ["Climate-only (A×M)","Full BSH (A×B×M)"];
                   position=:lb, framevisible=false, labelsize=10, padding=(2,2,2,2))
        return ax
    end

    # top row
    panel!(fig[1,1], R; title="Random")
    panel!(fig[1,2], C; title="Clustered")
    panel!(fig[1,3], F; title="Front-like")

    # bottom: biotic penalty
    ax = Axis(fig[2,1:3], xlabel="Area lost (fraction)",
              ylabel="Biotic penalty  (A×M − A×B×M)")
    lines!(ax, R.x, R.gap_biotic; color=:steelblue, linewidth=2, label="Random")
    lines!(ax, C.x, C.gap_biotic; color=:orange,    linewidth=2, label="Clustered")
    lines!(ax, F.x, F.gap_biotic; color=:seagreen,  linewidth=2, label="Front-like")
    axislegend(ax; position=:lt, framevisible=false, labelsize=10)

    # placebo curves (same colors, dashed)
    Rp = sweep_AM_vs_ABM_placebo(; grid,S,basal_frac,A_level,B_level,M_level,
                                 hl_kind=:random,   loss_fracs,seeds_pool,seeds_mask,sim_seed,
                                 nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)
    Cp = sweep_AM_vs_ABM_placebo(; grid,S,basal_frac,A_level,B_level,M_level,
                                 hl_kind=:clustered,loss_fracs,seeds_pool,seeds_mask,sim_seed,
                                 nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)
    Fp = sweep_AM_vs_ABM_placebo(; grid,S,basal_frac,A_level,B_level,M_level,
                                 hl_kind=:front,    loss_fracs,seeds_pool,seeds_mask,sim_seed,
                                 nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)

    lines!(ax, Rp.x, Rp.gap_biotic; color=:steelblue, linestyle=:dash, label="Random (placebo)")
    lines!(ax, Cp.x, Cp.gap_biotic; color=:orange,    linestyle=:dash, label="Clustered (placebo)")
    lines!(ax, Fp.x, Fp.gap_biotic; color=:seagreen,  linestyle=:dash, label="Front-like (placebo)")
    axislegend(ax; position=:lt, framevisible=false, labelsize=10)

    display(fig)
    return (; fig, random=R, clustered=C, front=F)
end

grid = make_grid(60, 60; seed=42)
# (1) Dashboard with placebo overlay
_ = plot_biotic_gap_dashboard(; grid,
       S=150, basal_frac=0.45,
       A_level=:divergent, B_level=:strong, M_level=:on,
       loss_fracs=0.2:0.1:0.8,
       seeds_pool=1:6, seeds_mask=1:6, sim_seed=1234,
       nseeds_cluster=6, front_axis=:x, front_noise=0.04,
       τA=0.5, τocc=0.35, T_frac_on=0.02)

# --- small, degree-preserving rewire (lower-mass prey only)
function rewire_metaweb_inplace!(rng::AbstractRNG, pool::SpeciesPool)
    order = sortperm(pool.masses)  # light→heavy
    rank  = zeros(Int, pool.S); rank[order] .= 1:pool.S
    for s in 1:pool.S
        pool.basal[s] && continue
        k = length(pool.E[s])
        k == 0 && continue
        # candidate prey: all species with lower mass than s
        cand = order[1:rank[s]-1]
        if isempty(cand)
            pool.E[s] = Int[]; continue
        end
        # sample exactly k distinct prey (or all if not enough)
        if length(cand) <= k
            pool.E[s] = copy(cand)
        else
            pool.E[s] = sample(rng, cand, k; replace=false)
        end
    end
    return pool
end

# --- placebo sweep (identical to your real sweep, but rewires after building the pool)
function sweep_AM_vs_ABM_placebo(; grid::Grid, S::Int, basal_frac::Float64,
        A_level::Symbol, B_level::Symbol, M_level::Symbol,
        hl_kind::Symbol, loss_fracs,
        seeds_pool, seeds_mask, sim_seed::Int,
        nseeds_cluster::Int, front_axis::Symbol, front_noise::Float64,
        τA::Float64, τocc::Float64, T_frac_on::Float64)

    make_keepmask = function(rng::AbstractRNG, kind::Symbol, keep_frac::Float64)
        kind === :random   && return random_mask(rng, grid.C, keep_frac)
        kind === :clustered&& return clustered_mask(rng, grid, keep_frac; nseeds=nseeds_cluster)
        kind === :front    && return frontlike_mask(rng, grid, keep_frac; axis=front_axis, noise=front_noise)
        error("Unknown hl_kind=$kind")
    end

    xs   = Float64[]; AMs = Float64[]; ABMs = Float64[]; GAPs = Float64[]
    for f in loss_fracs
        keepfrac = 1 - f
        vals_AM  = Float64[]; vals_ABM = Float64[]

        for ps in seeds_pool, ms in seeds_mask
            rng_pool = MersenneTwister(hash((sim_seed, :pool, ps)))
            pool     = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
            # REWIRE after build (placebo)
            rewire_metaweb_inplace!(MersenneTwister(hash((sim_seed,:rewire,ps))), pool)

            A        = abiotic_matrix(pool, grid)
            pars     = bam_from_axes(; B_level, M_level, τA, τocc)
            bam, mp  = pars.bam, pars.mp

            rng_mask = MersenneTwister(hash((sim_seed,:mask,ms,hl_kind,f)))
            keep     = make_keepmask(rng_mask, hl_kind, keepfrac)

            # Full BSH with rewired metaweb
            Pfull, _, _, _ = assemble_BAM(pool, grid, A, keep; bam=bam, mp=mp)
            bsh_full       = bsh1_per_species(Pfull, (A .>= τA), pool)
            mean_ABM_f     = mean(bsh_full[consumer_mask(pool)])

            # Climate-only (same as real)
            pars_AM = bam_from_axes(; B_level=:none, M_level=M_level, τA=τA, τocc=τocc)
            P_am, _, _, _ = assemble_BAM(pool, grid, A, keep; bam=pars_AM.bam, mp=pars_AM.mp)
            bsh_AM        = bsh1_per_species(P_am, (A .>= τA), pool)
            mean_AM_f     = mean(bsh_AM[consumer_mask(pool)])

            push!(vals_ABM, mean_ABM_f); push!(vals_AM, mean_AM_f)
        end

        push!(xs, f)
        μAM, μABM = mean(vals_AM), mean(vals_ABM)
        push!(AMs, μAM); push!(ABMs, μABM); push!(GAPs, μAM - μABM)
    end
    return (; x=xs, mean_AM=AMs, mean_ABM=ABMs, gap_biotic=GAPs)
end
