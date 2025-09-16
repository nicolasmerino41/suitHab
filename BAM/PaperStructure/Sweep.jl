# --- config ---
GEOMS    = (:random, :clustered, :front)
ALEVELS  = (:neutral, :intermediate, :divergent)
BLEVELS  = (:none, :soft, :strong)

function biotic_gap_at_fstar(; grid, S=120, basal_frac=0.45,
                             A_level=:neutral, B_level=:soft, M_level=:on,
                             hl_kind=:random, fstar=0.6,
                             seeds_pool=1:3, seeds_mask=1:3,
                             nseeds_cluster=6, front_axis=:x, front_noise=0.04,
                             τA=0.5, τocc=0.35, sim_seed=1234)

    keepfrac = 1 - fstar
    vals_gap = Float64[]

    for ps in seeds_pool, ms in seeds_mask
        rng_pool = MersenneTwister(hash((sim_seed,:pool,ps,A_level,B_level)))
        pool     = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
        A        = abiotic_matrix(pool, grid)

        pars_ABM = bam_from_axes(; B_level, M_level, τA, τocc)
        pars_AM  = bam_from_axes(; B_level=:none, M_level, τA, τocc)

        rng_mask = MersenneTwister(hash((sim_seed,:mask,ms,hl_kind,fstar)))
        keep     = keepmask_for(hl_kind, grid, keepfrac; rng=rng_mask,
                                nseeds_cluster=nseeds_cluster,
                                front_axis=front_axis, front_noise=front_noise)

        # assemble ABM
        P_abm,B_abm,_,M_abm = assemble_BAM(pool, grid, A, keep; bam=pars_ABM.bam, mp=pars_ABM.mp)
        stats_abm = species_stats(pool, grid, A, keep, P_abm, B_abm, M_abm, pars_ABM.bam)
        bsh_abm = mean(stats_abm.BSH[consumer_mask(pool)])

        # assemble AM
        P_am,B_am,_,M_am = assemble_BAM(pool, grid, A, keep; bam=pars_AM.bam, mp=pars_AM.mp)
        stats_am = species_stats(pool, grid, A, keep, P_am, B_am, M_am, pars_AM.bam)
        bsh_am = mean(stats_am.BSH[consumer_mask(pool)])

        push!(vals_gap, bsh_am - bsh_abm)
    end

    return mean(vals_gap)
end

function sweep_biotic_gap(; nx=60, ny=60, S=120, basal_frac=0.45,
                          A_levels=ALEVELS, B_levels=BLEVELS,
                          geoms=GEOMS, fstar=0.6)
    grid = make_grid(nx,ny; seed=42)
    results = Dict{Tuple{Symbol,Symbol,Symbol},Float64}()

    for Alev in A_levels, Blev in B_levels, geom in geoms
        gap = biotic_gap_at_fstar(; grid,S,basal_frac,
                                   A_level=Alev, B_level=Blev, M_level=:on,
                                   hl_kind=geom, fstar=fstar)
        results[(Alev,Blev,geom)] = Float64(gap)  # force scalar
    end
    return results
end

function sweep_biotic_gap_array(; nx=60, ny=60, S=120, basal_frac=0.45,
                                A_levels=ALEVELS, B_levels=BLEVELS,
                                geoms=GEOMS, fstar=0.6)

    grid = make_grid(nx,ny; seed=42)
    nA, nB, nG = length(A_levels), length(B_levels), length(geoms)

    # allocate directly: rows=A, cols=B, depth=geom
    gaps = Array{Float64}(undef, nA, nB, nG)

    for (ia,Alev) in enumerate(A_levels),
        (ib,Blev) in enumerate(B_levels),
        (ig,geom) in enumerate(geoms)

        gap = biotic_gap_at_fstar(; grid,S,basal_frac,
                                   A_level=Alev, B_level=Blev, M_level=:on,
                                   hl_kind=geom, fstar=fstar)
        gaps[ia, ib, ig] = gap
    end

    return (gaps=gaps, A_levels=A_levels, B_levels=B_levels, geoms=geoms, fstar=fstar)
end

# --- visualization ---
function plot_gap_heatmap(results; fstar=0.6)
    fig = Figure(; size=(1000, 400))

    for (k, geom) in enumerate(GEOMS)
        ax = Axis(fig[1,k], title="HL=$(geom)",
                  xlabel="B-level", ylabel="A-level")

        # Safe numeric conversion
        mat = [try
                   Float64(results[(A,B,geom)])
               catch
                   NaN
               end
               for A in ALEVELS, B in BLEVELS]

        heatmap!(ax, 1:length(BLEVELS), 1:length(ALEVELS), mat;
                 colormap=:balance, colorrange=(-0.2,0.2))

        ax.xticks = (1:length(BLEVELS), string.(BLEVELS))
        ax.yticks = (1:length(ALEVELS), string.(ALEVELS))

        Colorbar(fig[1,k], limits=(-0.2,0.2),
                 label="Gap (AM−ABM)", vertical=false)
    end

    Label(fig[0,:],
          @sprintf("Biotic gap at f* = %.2f", fstar);
          fontsize=18)

    display(fig)
    return fig
end

function plot_gap_heatmap_array(sweepres)
    gaps, A_levels, B_levels, geoms, fstar =
        sweepres.gaps, sweepres.A_levels, sweepres.B_levels, sweepres.geoms, sweepres.fstar

    fig = Figure(; size=(1000, 400))

    for (ig, geom) in enumerate(geoms)
        ax = Axis(fig[1,ig], title="HL=$(geom)", xlabel="B-level", ylabel="A-level")

        mat = gaps[:,:,ig]  # slice for this geometry

        heatmap!(ax, 1:length(B_levels), 1:length(A_levels), mat;
                 colormap=:balance, colorrange=(-0.2,0.2))

        ax.xticks = (1:length(B_levels), string.(B_levels))
        ax.yticks = (1:length(A_levels), string.(A_levels))

        Colorbar(fig[2,ig], limits=(-0.2,0.2), label="Gap (AM−ABM)", vertical=false)
    end

    Label(fig[0,:], @sprintf("Biotic gap at f* = %.2f", fstar); fontsize=18)

    display(fig)
    return fig
end


# --- run ---
res = sweep_biotic_gap_array()
fig = plot_gap_heatmap_array(res)
