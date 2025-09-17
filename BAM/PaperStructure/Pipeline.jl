# ---------- helpers you likely already have ----------
# make_grid, random_mask, clustered_mask, frontlike_mask
# build_pool_from_axes(rng; S, basal_frac, A_level, B_level)
# bam_from_axes(; B_level, M_level, α=1, τA=0.35, τocc=0.42)  -> (bam, mp)
# abiotic_matrix(pool, grid)
# assemble_BAM(pool, grid, A, keep; bam, mp) -> (P, Bsup, Score, M)
# species_stats(pool, grid, A, keep, P, Bsup, M, bam) -> BAMStats
# shapley_per_species(pool, grid, A, keep0, keep1; bam, mp) -> (Shapley3 per sp, S0, S1)
# consumer_mask(pool)

# ---------- placebo (rewire) ----------
"""
placebo_rewire(pool; rng)
Rewire each consumer's prey list to a random set of same length drawn from
all lighter species (preserves diet length distribution and basal fraction).
"""
function placebo_rewire(pool; rng=MersenneTwister(0))
    S = pool.S
    E2 = [Int[] for _ in 1:S]
    order = sortperm(pool.masses) # light -> heavy
    lighter = [order[1:i-1] for i in 1:S]
    for s in 1:S
        if pool.basal[s] || isempty(pool.E[s])
            E2[s] = Int[]
        else
            L = length(pool.E[s])
            cand = lighter[s]
            if isempty(cand)
                E2[s] = Int[]
            elseif L >= length(cand)
                E2[s] = copy(cand)
            else
                E2[s] = sample(rng, cand, L; replace=false)
            end
        end
    end
    SpeciesPool(S, pool.masses, pool.basal, pool.mu, pool.b, E2)
end

# ---------- common evaluation primitives ----------
const GEOMS = (:random, :clustered, :front)
const GLABEL = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")

function keepmask_for(hl::Symbol, grid, keep_frac; rng, nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04)
    hl === :random    && return random_mask(rng, grid.C, keep_frac)
    hl === :clustered && return clustered_mask(rng, grid, keep_frac; nseeds=nseeds_cluster)
    hl === :front     && return frontlike_mask(rng, grid, keep_frac; axis=front_axis, noise=front_noise)
    error("unknown HL $hl")
end

"""
ΔBSH curves for a SPECIFIC pool and BAM settings (we keep pool fixed to compare lenses).
Returns Dict(hl => (x, Δ̄(f))) where Δ̄ is mean ΔBSH over consumers vs f.
"""
function curves_ΔBSH(grid, pool; bam, mp, loss_fracs=0.2:0.1:0.8,
                     seeds_mask=1:6, nseeds_cluster=6, front_axis=:x, front_noise=0.04,
                     sim_seed=1234)
    A = abiotic_matrix(pool, grid)
    base_keep = trues(grid.C)
    P0, B0, _, M0 = assemble_BAM(pool, grid, A, base_keep; bam=bam, mp=mp)
    S0 = species_stats(pool, grid, A, base_keep, P0, B0, M0, bam)
    cons = .!pool.basal

    out = Dict{Symbol, NamedTuple}()
    xs = collect(loss_fracs)
    for hl in GEOMS
        ys = Float64[]
        for f in xs
            keepfrac = 1 - f
            vals = Float64[]
            for ms in seeds_mask
                rng_mask = MersenneTwister(hash((sim_seed,:mask,ms,hl,f)))
                keep = keepmask_for(hl, grid, keepfrac; rng=rng_mask,
                                    nseeds_cluster=nseeds_cluster, front_axis=front_axis, front_noise=front_noise)
                P1, B1, _, M1 = assemble_BAM(pool, grid, A, keep; bam=bam, mp=mp)
                S1 = species_stats(pool, grid, A, keep, P1, B1, M1, bam)
                push!(vals, mean(S1.BSH[cons] .- S0.BSH[cons]))
            end
            push!(ys, mean(vals))
        end
        out[hl] = (x=xs, y=ys)
    end
    out
end

# ---------- Figure A: AM vs ABM with placebo + DiD/WU inset ----------
function figure_A_AM_vs_ABM(; nx=60, ny=60, S=150, basal_frac=0.45,
        A_level=:intermediate, B_level=:strong, M_level=:on,
        loss_fracs=0.2:0.1:0.8, seeds_pool=1:4, seeds_mask=1:6,
        nseeds_cluster=6, front_axis=:x, front_noise=0.04,
        τA=0.5, τocc=0.35, T=6, sim_seed=1234)

    grid = make_grid(nx,ny; seed=42)
    rng_pool = MersenneTwister(hash((sim_seed,:pool,1)))
    pool_ABM = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
    pool_AM  = pool_ABM
    pool_PL  = placebo_rewire(pool_ABM; rng=MersenneTwister(hash((sim_seed,:placebo))))

    (bam_abm, mp_abm) = bam_from_axes(; B_level, M_level, τA=τA, τocc=τocc)
    mp_abm = MovementParams(; mode=mp_abm.mode, T=T)
    (bam_am,  mp_am)  = bam_from_axes(; B_level=:none, M_level, τA=τA, τocc=τocc)
    mp_am  = MovementParams(; mode=mp_am.mode, T=T)

    cur_AM   = curves_ΔBSH(grid, pool_AM; bam=bam_am,  mp=mp_am,
                           loss_fracs, seeds_mask, nseeds_cluster, front_axis, front_noise, sim_seed)
    cur_ABM  = curves_ΔBSH(grid, pool_ABM; bam=bam_abm, mp=mp_abm,
                           loss_fracs, seeds_mask, nseeds_cluster, front_axis, front_noise, sim_seed)
    cur_PL   = curves_ΔBSH(grid, pool_PL;  bam=bam_abm, mp=mp_abm,
                           loss_fracs, seeds_mask, nseeds_cluster, front_axis, front_noise, sim_seed)

    xs = cur_AM[:random].x
    DiD = Dict(g => cur_ABM[g].y .- cur_AM[g].y for g in GEOMS)
    WU  = map(i -> maximum([cur_ABM[g].y[i] for g in GEOMS]) -
                   maximum([cur_AM[g].y[i] for g in GEOMS]), eachindex(xs))

    fig = Figure(; size=(1800, 700))
    Label(fig[0,1:2],
          @sprintf("HL effect with vs without biotic — A=%s, B=%s, M=%s  (placebo dashed)",
                   String(A_level), String(B_level), String(M_level));
          fontsize=20)

    # Left: AM
    axL = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="ΔBSH (mean consumers)",
               title="Climate-only (AM: A×M)")
    for g in GEOMS
        lines!(axL, xs, cur_AM[g].y; label=GLABEL[g])
    end
    axislegend(axL; position=:lb, framevisible=false)

    # Right: ABM + placebo
    axR = Axis(fig[1,2], xlabel="Area lost (fraction)", ylabel="ΔBSH (mean consumers)",
               title="Full ABM (A×B×M)  +  placebo (dashed)")
    cols = Dict(:random=>:royalblue, :clustered=>:orange, :front=>:seagreen)
    for g in GEOMS
        lines!(axR, xs, cur_ABM[g].y; color=cols[g], label=GLABEL[g])
        lines!(axR, xs, cur_PL[g].y; color=cols[g], linestyle=:dash)
    end
    axislegend(axR; position=:lb, framevisible=false)

    # Inset: DiD + WU
    axI = Axis(fig[2,1], halign=:right, valign=:top, width=320, height=220,
               title="Biotic amplification: DiD(f) and WU(f)")
    for g in GEOMS
        lines!(axI, xs, DiD[g]; color=cols[g], label=GLABEL[g])
    end
    lines!(axI, xs, WU; color=:black, linestyle=:dot, label="WU (worst-case uplift)")
    axislegend(axI; position=:rb, framevisible=false, nbanks=1)

    display(fig)
    return (; fig, xs, cur_AM, cur_ABM, cur_PL, DiD, WU, grid, pool_ABM, bam_abm, mp_abm, bam_am, mp_am)
end

# ---------- Figure B: Elasticity ternary at f* ----------
# ternary coordinate
@inline function ternary_xy(sA, sB, sM)
    x = sB + 0.5*sM
    y = (sqrt(3)/2) * sM
    return x, y
end

function figure_B_elasticity_ternary(; grid, pool, bam, mp, fstar=0.6,
        nseeds_cluster=6, front_axis=:x, front_noise=0.04, τA=bam.τA, τocc=bam.τocc,
        sim_seed=1234)

    A = abiotic_matrix(pool, grid)
    base_keep = trues(grid.C)
    xs = [fstar]; keepfrac = 1 - fstar

    shares = Dict{Symbol,Tuple{Float64,Float64,Float64}}()
    for hl in GEOMS
        rng_mask = MersenneTwister(hash((sim_seed,:mask,hl,fstar)))
        keep = keepmask_for(hl, grid, keepfrac; rng=rng_mask,
                            nseeds_cluster=nseeds_cluster, front_axis=front_axis, front_noise=front_noise)
        shap, _, _ = shapley_per_species(pool, grid, A, base_keep, keep; bam=bam, mp=mp)
        dA = mean(abs.(getfield.(shap, :dA)))
        dB = mean(abs.(getfield.(shap, :dB)))
        dM = mean(abs.(getfield.(shap, :dM)))
        S = dA + dB + dM + 1e-12
        shares[hl] = (dA/S, dB/S, dM/S)
    end

    fig = Figure(; size=(900, 600))
    Label(fig[0,1], @sprintf("Elasticity shares at f* = %.2f", fstar); fontsize=18)

    ax = Axis(fig[1,1], title="A/B/M shares by geometry", xlabel="Biotic →", ylabel="Movement ↑")
    # draw triangle frame
    xsf = [0.0, 1.0, 0.5, 0.0]; ysf = [0.0, 0.0, sqrt(3)/2, 0.0]; lines!(ax, xsf, ysf; color=:black)

    cols = Dict(:random=>:royalblue, :clustered=>:orange, :front=>:seagreen)
    for hl in GEOMS
        sA, sB, sM = shares[hl]
        x, y = ternary_xy(sA, sB, sM)
        scatter!(ax, [x], [y]; color=cols[hl], markersize=16, label=GLABEL[hl])
        text!(ax, x, y, text=@sprintf("B=%.2f", sB), align=(:left,:bottom), fontsize=10, color=cols[hl])
    end
    axislegend(ax; position=:rb, framevisible=false)
    display(fig)
    return (; fig, shares)
end

# ---------- Figure C: Rank-shift grid ----------
"""
Rank-shift grid at loss f (shows AM worst vs ABM worst per (A,B) cell).
Minimal nudge pre-baked: T=6 and β bumps come from bam_from_axes you adjusted.
"""
function figure_C_rank_shifts(; nx=60, ny=60, S=150, basal_frac=0.45,
        A_levels=(:neutral,:intermediate,:divergent), B_levels=(:none,:soft,:strong),
        M_level=:on, f=0.5, τA=0.5, τocc=0.35, T=6, nseeds_cluster=6, sim_seed=1234)

    grid = make_grid(nx,ny; seed=42)
    keepfrac = 1 - f

    fig = Figure(; size=(1200, 500))
    Label(fig[0,1], @sprintf("Rank shifts at loss=%.2f — M=%s", f, String(M_level)); fontsize=20)

    # layout 3x3
    for (i,Alev) in enumerate(A_levels), (j,Blev) in enumerate(B_levels)
        ax = Axis(fig[i, j], title=@sprintf("A=%s, B=%s", String(Alev), String(Blev)))
        hidedecorations!(ax); hidespines!(ax)
        # build pool
        rng_pool = MersenneTwister(hash((sim_seed,:pool,Alev,Blev)))
        pool = build_pool_from_axes(rng_pool; S, basal_frac, A_level=Alev, B_level=Blev)
        A = abiotic_matrix(pool, grid)

        # BAMs
        (bam_abm, mp_abm) = bam_from_axes(; B_level=Blev, M_level, τA, τocc); mp_abm = MovementParams(; mode=mp_abm.mode, T=T)
        (bam_am,  mp_am)  = bam_from_axes(; B_level=:none, M_level, τA, τocc);   mp_am  = MovementParams(; mode=mp_am.mode,  T=T)

        # find worst geometry
        function worst_for(bam, mp)
            base_keep = trues(grid.C)
            P0,B0,_,M0 = assemble_BAM(pool, grid, A, base_keep; bam, mp)
            S0 = species_stats(pool, grid, A, base_keep, P0,B0,M0, bam)
            cons = .!pool.basal
            y = Dict{Symbol,Float64}()
            for hl in GEOMS
                rng_mask = MersenneTwister(hash((sim_seed,:mask,hl,Alev,Blev)))
                keep = keepmask_for(hl, grid, keepfrac; rng=rng_mask, nseeds_cluster=nseeds_cluster)
                P1,B1,_,M1 = assemble_BAM(pool, grid, A, keep; bam, mp)
                S1 = species_stats(pool, grid, A, keep, P1,B1,M1, bam)
                y[hl] = mean(S1.BSH[cons] .- S0.BSH[cons])
            end
            # most negative is worst
            return argmin([y[:random], y[:clustered], y[:front]]) == 1 ? :random :
                   argmin([y[:random], y[:clustered], y[:front]]) == 2 ? :clustered : :front
        end

        g_off = worst_for(bam_am,  mp_am)
        g_on  = worst_for(bam_abm, mp_abm)

        # draw markers (left = OFF, right = ON)
        x0, x1, yv = 0.3, 0.7, 0.5
        mshape(g) = g==:random ? :circle : (g==:clustered ? :rect : :utriangle)
        scatter!(ax, [x0], [yv]; marker=mshape(g_off), color=:gray, markersize=14)
        scatter!(ax, [x1], [yv]; marker=mshape(g_on),  color=:black, markersize=14)
        text!(ax, x0, yv-0.18, text="B OFF", align=(:center,:top))
        text!(ax, x1, yv-0.18, text="B ON",  align=(:center,:top))
        if g_off != g_on
            arrows!(ax, [x0], [yv], [x1-x0], [0.0]; arrowsize=10, linewidth=2, color=:black)
        end
    end
    display(fig)
    fig
end

# ---------- Figure D: Management toy policy ----------
"""
Score cells at baseline and keep top-K to achieve loss f; compute ΔBSH.
policy=:random uses random keep; :AM uses mean baseline A×M; :ABM uses mean baseline A×B×M.
Optionally compare ABM with placebo-rewired metaweb.
"""
function policy_curves(; grid, pool, bam_abm, mp_abm, bam_am, mp_am,
        loss_fracs=0.2:0.1:0.8, seeds_mask=1:6, nseeds_cluster=6, sim_seed=1234)

    A = abiotic_matrix(pool, grid)
    base_keep = trues(grid.C)
    cons = .!pool.basal

    # baseline stats for scoring
    P0a,B0a,_,M0a = assemble_BAM(pool, grid, A, base_keep; bam=bam_am,  mp=mp_am)
    S0a = species_stats(pool, grid, A, base_keep, P0a,B0a,M0a, bam_am)
    score_AM  = mean(S0a.A .* (S0a.M .^ (mp_am.mode==:none ? 0.0 : 1.0))) # not used directly; do cellwise
    # better: compute per-cell mean scores
    function cell_scores(bam, mp)
        # approximate mean over consumers by counting fraction passing each factor
        # use Score from assemble at baseline
        _, _, Score, _ = assemble_BAM(pool, grid, A, base_keep; bam=bam, mp=mp)
        mean(Score; dims=1)[:]   # 1×C → Vector{Float64}
    end
    S_AM_cells  = cell_scores(bam_am,  mp_am)   # climate-only score
    S_ABM_cells = cell_scores(bam_abm, mp_abm)  # interaction-aware score

    # helper to evaluate ΔBSH when keeping top-K by given score
    function delta_by_score(scorecells)
        xs = Float64[]; ys = Float64[]
        for f in loss_fracs
            keepfrac = 1 - f
            K = clamp(round(Int, keepfrac*grid.C), 0, grid.C)
            keep_idx = sortperm(scorecells; rev=true)[1:K]
            keep = falses(grid.C); keep[keep_idx] .= true
            # evaluate ABM outcome under this keep
            P0,B0,_,M0 = assemble_BAM(pool, grid, A, base_keep; bam=bam_abm, mp=mp_abm)
            S0 = species_stats(pool, grid, A, base_keep, P0,B0,M0, bam_abm)
            P1,B1,_,M1 = assemble_BAM(pool, grid, A, keep;      bam=bam_abm, mp=mp_abm)
            S1 = species_stats(pool, grid, A, keep, P1,B1,M1, bam_abm)
            push!(xs, f)
            push!(ys, mean(S1.BSH[cons] .- S0.BSH[cons]))
        end
        (x=xs, y=ys)
    end

    cur_rand = Dict{Symbol,NamedTuple}()
    # random baseline (averaged over masks)
    xs = collect(loss_fracs); yr = zeros(length(xs))
    for (k,f) in enumerate(xs)
        keepfrac = 1 - f
        vals = Float64[]
        for ms in seeds_mask
            rng = MersenneTwister(hash((sim_seed,:prand,ms,f)))
            keep = random_mask(rng, grid.C, keepfrac)
            P0,B0,_,M0 = assemble_BAM(pool, grid, A, base_keep; bam=bam_abm, mp=mp_abm)
            S0 = species_stats(pool, grid, A, base_keep, P0,B0,M0, bam_abm)
            P1,B1,_,M1 = assemble_BAM(pool, grid, A, keep;      bam=bam_abm, mp=mp_abm)
            S1 = species_stats(pool, grid, A, keep, P1,B1,M1, bam_abm)
            push!(vals, mean(S1.BSH[cons] .- S0.BSH[cons]))
        end
        yr[k] = mean(vals)
    end
    cur_rand = (x=xs, y=yr)

    cur_AM  = delta_by_score(S_AM_cells)
    cur_ABM = delta_by_score(S_ABM_cells)

    # placebo policy (interaction-aware on rewired metaweb)
    pool_PL  = placebo_rewire(pool; rng=MersenneTwister(hash((sim_seed,:pl)))) 
    (bam_pl, mp_pl) = (bam_abm, mp_abm) # same parameters, different metaweb
    # score on placebo baseline
    function cell_scores_pool(poolX)
        A = abiotic_matrix(poolX, grid)
        base_keep = trues(grid.C)
        _,_,Score,_ = assemble_BAM(poolX, grid, A, base_keep; bam=bam_pl, mp=mp_pl)
        mean(Score; dims=1)[:]
    end
    S_PL_cells = cell_scores_pool(pool_PL)
    # evaluate Δ using ABM **on placebo** too (apple-to-apple within placebo world)
    function delta_by_score_placebo(scorecells)
        Apl = abiotic_matrix(pool_PL, grid)
        xs = Float64[]; ys = Float64[]
        base_keep = trues(grid.C)
        P0,B0,_,M0 = assemble_BAM(pool_PL, grid, Apl, base_keep; bam=bam_pl, mp=mp_pl)
        S0 = species_stats(pool_PL, grid, Apl, base_keep, P0,B0,M0, bam_pl)
        cons_pl = .!pool_PL.basal
        for f in loss_fracs
            keepfrac = 1 - f
            K = clamp(round(Int, keepfrac*grid.C), 0, grid.C)
            keep_idx = sortperm(scorecells; rev=true)[1:K]
            keep = falses(grid.C); keep[keep_idx] .= true
            P1,B1,_,M1 = assemble_BAM(pool_PL, grid, Apl, keep;      bam=bam_pl, mp=mp_pl)
            S1 = species_stats(pool_PL, grid, Apl, keep, P1,B1,M1, bam_pl)
            push!(xs, f)
            push!(ys, mean(S1.BSH[cons_pl] .- S0.BSH[cons_pl]))
        end
        (x=xs, y=ys)
    end
    cur_PL = delta_by_score_placebo(S_PL_cells)

    # plot
    fig = Figure(; size=(1600, 550))
    Label(fig[0,1], "Policy-value: ABM outcome under different prioritizations"; fontsize=20)

    ax = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="ΔBSH (mean consumers)")
    lines!(ax, cur_rand.x, cur_rand.y, label="Random")
    lines!(ax, cur_AM.x,  cur_AM.y,  label="Climate-only score (A×M)")
    lines!(ax, cur_ABM.x, cur_ABM.y, label="Interaction-aware (A×B×M)")
    lines!(ax, cur_PL.x,  cur_PL.y;  linestyle=:dash, color=:gray, label="Interaction-aware (placebo)")
    axislegend(ax; position=:lb, framevisible=false)

    display(fig)
    return (; fig, cur_rand, cur_AM, cur_ABM, cur_PL)
end

# -------------------- RUN: generate the four figures --------------------
# Canonical regime + minimal nudge
FIGA = figure_A_AM_vs_ABM(; A_level=:divergent, B_level=:soft, M_level=:on,
                                T=6, loss_fracs=0.2:0.1:0.8, seeds_pool=1:4, seeds_mask=1:6)

# Use the pool and BAM returned by Figure A to keep objects consistent
FIGB = figure_B_elasticity_ternary(; grid=FIGA.grid, pool=FIGA.pool_ABM,
                                         bam=FIGA.bam_abm, mp=FIGA.mp_abm, fstar=0.60)

FIGC = figure_C_rank_shifts(; A_levels=(:neutral,:intermediate,:divergent),
                                  B_levels=(:none,:soft,:strong), M_level=:on,
                                  f=0.50, T=6)

FIGD = policy_curves(; grid=FIGA.grid, pool=FIGA.pool_ABM,
                            bam_abm=FIGA.bam_abm, mp_abm=FIGA.mp_abm,
                            bam_am=FIGA.bam_am,   mp_am=FIGA.mp_am,
                            loss_fracs=0.2:0.1:0.8)
