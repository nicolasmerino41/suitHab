##############  ELASTICITY ANALYSIS  #########################################

# --- per-species elasticities from two states (baseline vs mask) -----------
# Inputs: species-level stats S0, S1 (fields A,B,M,F,BSH) and BAM params
# Returns: eA,eB,eM vectors (one per species)
function _elasticities_per_species(S0::BAMStats, S1::BAMStats, bam::BAMParams)
    S = length(S0.A)
    eA = zeros(Float64, S)
    eB = zeros(Float64, S)
    eM = zeros(Float64, S)
    ϵ = 1e-12
    @inbounds for s in 1:S
        A0=A1=0.0; B0=B1=0.0; M0=M1=0.0
        A0 = S0.A[s]; A1 = S1.A[s]
        B0 = S0.B[s]; B1 = S1.B[s]
        M0 = S0.M[s]; M1 = S1.M[s]
        Abar = 0.5*(A0+A1); Bbar = 0.5*(B0+B1); Mbar = 0.5*(M0+M1)
        dA = A1 - A0; dB = B1 - B0; dM = M1 - M0
        eA[s] = abs(dA) / max(Abar, ϵ) * 1.0
        eB[s] = abs(dB) / max(Bbar, ϵ) * bam.β
        eM[s] = abs(dM) / max(Mbar, ϵ) * bam.μ
    end
    return eA, eB, eM
end

# --- summarize elasticities for ONE mask (one HL geometry) ------------------
# Averages over pool & mask seeds; returns mean(E_A), mean(E_B), mean(E_M) for consumers.
function _elasticity_for_geometry(; grid::Grid, S::Int, basal_frac::Float64,
        A_level::Symbol, B_level::Symbol, M_level::Symbol,
        hl_kind::Symbol, keepfrac::Float64,
        seeds_pool, seeds_mask, sim_seed::Int,
        nseeds_cluster::Int, front_axis::Symbol, front_noise::Float64,
        τA::Float64, τocc::Float64, T_frac_on::Float64)

    # BAM + movement (scale T with map size when on)
    pars = bam_from_axes(; B_level, M_level, τA, τocc)
    bam = pars.bam
    T = M_level === :on ? max(8, round(Int, T_frac_on * grid.C)) : 8
    mp = MovementParams(; mode=(M_level===:on ? :component : :none), T=T)

    base_keep = trues(grid.C)
    EA = Float64[]; EB = Float64[]; EM = Float64[]

    for ps in seeds_pool, ms in seeds_mask
        rng_pool = thread_rng(sim_seed, :pool, ps, hl_kind)
        rng_mask = thread_rng(sim_seed, :mask, ms, hl_kind)

        pool = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
        A    = abiotic_matrix(pool, grid)

        keep = hl_kind === :random    ? random_mask(rng_mask, grid.C, keepfrac) :
               hl_kind === :clustered ? clustered_mask(rng_mask, grid, keepfrac; nseeds=nseeds_cluster) :
               hl_kind === :front     ? frontlike_mask(rng_mask, grid, keepfrac; axis=front_axis, noise=front_noise) :
               error("Unknown hl_kind")

        # species-level stats baseline vs mask
        P0,B0,_,M0 = assemble_BAM(pool, grid, A, base_keep; bam=bam, mp=mp)
        S0 = species_stats(pool, grid, A, base_keep, P0,B0,M0, bam)
        P1,B1,_,M1 = assemble_BAM(pool, grid, A, keep;       bam=bam, mp=mp)
        S1 = species_stats(pool, grid, A, keep,      P1,B1,M1, bam)

        eA,eB,eM = _elasticities_per_species(S0, S1, bam)
        cons = consumer_mask(pool)
        push!(EA, mean(eA[cons])); push!(EB, mean(eB[cons])); push!(EM, mean(eM[cons]))
    end

    E_A = mean(EA); E_B = mean(EB); E_M = mean(EM)
    E_sum = max(E_A + E_B + E_M, 1e-12)
    return (; E_A, E_B, E_M, S_A=E_A/E_sum, S_B=E_B/E_sum, S_M=E_M/E_sum)
end

# --- summarize ALL 18 combos at a chosen loss, for all 3 geometries ---------
function elasticity_summarize_all(; nx::Int=60, ny::Int=60, S::Int=150, basal_frac::Float64=0.45,
        loss_pick::Float64=0.6, seeds_pool=1:5, seeds_mask=1:5, sim_seed::Int=1234,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        τA::Float64=0.5, τocc::Float64=0.35, T_frac_on::Float64=0.02)

    grid = make_grid(nx, ny; seed=42)
    keepfrac = 1 - loss_pick
    names = collect(keys(COMBOS))
    results = Vector{NamedTuple}(undef, length(names))

    Threads.@threads for i in eachindex(names)
        nm   = names[i]
        axes = COMBOS[nm]
        per_g = Dict{Symbol,NamedTuple}()
        for g in (:random, :clustered, :front)
            per_g[g] = _elasticity_for_geometry(; grid, S, basal_frac,
                            A_level=axes.A_level, B_level=axes.B_level, M_level=axes.M_level,
                            hl_kind=g, keepfrac, seeds_pool, seeds_mask, sim_seed,
                            nseeds_cluster, front_axis, front_noise,
                            τA, τocc, T_frac_on)
        end
        results[i] = (; name=nm, A_level=axes.A_level, B_level=axes.B_level, M_level=axes.M_level, per_hl=per_g)
    end
    return (; grid, names, results, loss=loss_pick)
end

# --- ternary plot of elasticity SHARES + bar of "elasticity geometry sensitivity" ---
@inline function _ternary_xy(sA, sB, sM)
    x = sB + 0.5*sM
    y = (√3/2)*sM
    return x, y
end

function plot_elasticity_summary(sumres; title="Elasticity shares at fixed loss")
    names   = sumres.names
    results = sumres.results
    geoms   = (:random, :clustered, :front)
    labels  = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")

    # encodings
    colA = Dict(:neutral=>RGBAf(0.55,0.55,0.55,1.0),
                :intermediate=>RGBAf(0.25,0.45,0.85,1.0),
                :divergent=>RGBAf(0.85,0.35,0.25,1.0))
    markB = Dict(:none=>:circle, :soft=>:utriangle, :strong=>:rect)

    # compute "elasticity geometry sensitivity" (EGS) = range of E_sum across geometries
    EGS = Float64[]
    labels_ord = Symbol[]

    fig = Figure(; size=(1250, 760))
    Label(fig[0, :], "$(title) — loss=$(round(sumres.loss; digits=2))"; fontsize=16, padding=(0,0,10,0))

    # Top: ternary per geometry
    for (col,g) in enumerate(geoms)
        ax = Axis(fig[1,col], title=labels[g], xlabel="Biotic →", ylabel="Movement ↑")
        xs = [0.0, 1.0, 0.5, 0.0]; ys = [0.0, 0.0, √3/2, 0.0]; lines!(ax, xs, ys)
        text!(ax, 0.02, -0.03, text="Abiotic", align=(:left,:bottom))
        text!(ax, 0.98, -0.03, text="Biotic",  align=(:right,:bottom))
        text!(ax, 0.50,  √3/2 + 0.03, text="Movement", align=(:center,:bottom))

        # also collect E_sum to scale marker sizes consistently
        Esums = [ (r.per_hl[:random].E_A + r.per_hl[:random].E_B + r.per_hl[:random].E_M,
                   r.per_hl[:clustered].E_A + r.per_hl[:clustered].E_B + r.per_hl[:clustered].E_M,
                   r.per_hl[:front].E_A + r.per_hl[:front].E_B + r.per_hl[:front].E_M) for r in results]
        maxE = maximum([maximum(t) for t in Esums]) + 1e-12

        for r in results
            p = r.per_hl[g]
            x,y = _ternary_xy(p.S_A, p.S_B, p.S_M)
            Esum = p.E_A + p.E_B + p.E_M
            ms = 6 + 40*Esum / maxE
            color = colA[r.A_level]; marker = markB[r.B_level]
            strokew = r.M_level === :on ? 1.5 : 0.5
            scatter!(ax, [x], [y]; markersize=ms, marker=marker,
                     color=color, strokecolor=:black, strokewidth=strokew)
        end
        hidespines!(ax, :r, :t)
    end

    custom_legend(fig, colA, markB)

    # Bottom: EGS bars
    ax2 = Axis(
        fig[2,1:3], xlabel="Combo", ylabel="EGS = range of ∑ elasticities",
        title="Geometry sensitivity (elasticity)",
        xticklabelrotation=π/6
    )
    es = Float64[]; nms = Symbol[]
    for r in results
        eR = r.per_hl[:random];    sR = eR.E_A + eR.E_B + eR.E_M
        eC = r.per_hl[:clustered]; sC = eC.E_A + eC.E_B + eC.E_M
        eF = r.per_hl[:front];     sF = eF.E_A + eF.E_B + eF.E_M
        push!(es, maximum((sR,sC,sF)) - minimum((sR,sC,sF)))
        push!(nms, r.name)
    end
    ord = sortperm(es; rev=true)
    barpositions = 1:length(ord)
    barplot!(ax2, barpositions, es[ord])
    ax2.xticks = (barpositions, string.(nms[ord]))

    display(fig)
    return fig
end

# If not already in scope:
@inline function _ternary_xy(sA, sB, sM)
    x = sB + 0.5*sM
    y = (√3/2)*sM
    return x, y
end

"""
compare_raw_vs_elasticity(sum_raw, sum_elast; title)
- sum_raw   = result from summarize_all_combos(...)
- sum_elast = result from elasticity_summarize_all(...)
Produces a 2×3 figure: Row 1 = raw ΔF shares, Row 2 = elasticity shares, for
Random / Clustered / Front. Encodings:
- Color = Abiotic level (neutral/intermediate/divergent)
- Marker = Biotic level (none/soft/strong)
- Thick outline = Movement:on
- Marker size ∝ magnitude (|ΔF| in row 1, ∑ elasticities in row 2)
"""
function compare_raw_vs_elasticity(sum_raw; sum_elast=sum_raw, title="Raw ΔF vs Elasticity shares")
    # encodings
    geoms  = (:random, :clustered, :front)
    glabel = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")
    colA = Dict(:neutral=>RGBAf(0.55,0.55,0.55,1.0),
                :intermediate=>RGBAf(0.25,0.45,0.85,1.0),
                :divergent=>RGBAf(0.85,0.35,0.25,1.0))
    markB = Dict(:none=>:circle, :soft=>:utriangle, :strong=>:rect)

    # Helper to get max scaling per row
    function max_scale_row_raw(results)
        maximum([maximum(abs.([r.per_hl[g].dF for g in geoms])) for r in results]) + 1e-12
    end
    function max_scale_row_elast(results)
        maximum([maximum([r.per_hl[g].E_A + r.per_hl[g].E_B + r.per_hl[g].E_M for g in geoms])
                 for r in results]) + 1e-12
    end

    Rraw   = sum_raw.results
    Relast = sum_elast.results
    maxΔF  = max_scale_row_raw(Rraw)
    maxE   = max_scale_row_elast(Relast)

    fig = Figure(; size=(1320, 650))
    Label(fig[0, :], title * @sprintf(" — loss=%.2f", get(sum_elast, :loss, NaN)); fontsize=16, padding=(0,0,10,0))

    # Draw a small ternary frame
    function frame!(ax)
        xs = [0.0, 1.0, 0.5, 0.0]; ys = [0.0, 0.0, √3/2, 0.0]
        lines!(ax, xs, ys)
        text!(ax, 0.02, -0.03, text="Abiotic", align=(:left,:bottom))
        text!(ax, 0.98, -0.03, text="Biotic",  align=(:right,:bottom))
        text!(ax, 0.50,  √3/2 + 0.03, text="Movement", align=(:center,:bottom))
        hidespines!(ax, :r, :t)
    end

    # Row 1: RAW ΔF shares
    for (col,g) in enumerate(geoms)
        ax = Axis(fig[1,col], title=glabel[g]*" — raw ΔF shares", xlabel="Biotic →", ylabel="Movement ↑")
        frame!(ax)
        for r in Rraw
            p = r.per_hl[g]
            sA, sB, sM = p.sA, p.sB, p.sM
            x,y = _ternary_xy(sA, sB, sM)
            mag = abs(p.dF)
            ms  = 6 + 40 * mag / maxΔF
            color = colA[r.A_level]; marker = markB[r.B_level]
            strokew = (r.M_level === :on) ? 1.5 : 0.5
            scatter!(ax, [x], [y]; markersize=ms, marker=marker,
                     color=color, strokecolor=:black, strokewidth=strokew)
        end
    end

    # Row 2: ELASTICITY shares
    for (col,g) in enumerate(geoms)
        ax = Axis(fig[2,col], title=glabel[g]*" — elasticity shares", xlabel="Biotic →", ylabel="Movement ↑")
        frame!(ax)
        for r in Relast
            p = r.per_hl[g]
            x,y = _ternary_xy(p.S_A, p.S_B, p.S_M)
            Esum = p.E_A + p.E_B + p.E_M
            ms  = 6 + 40 * Esum / maxE
            color = colA[r.A_level]; marker = markB[r.B_level]
            strokew = (r.M_level === :on) ? 1.5 : 0.5
            scatter!(ax, [x], [y]; markersize=ms, marker=marker,
                     color=color, strokecolor=:black, strokewidth=strokew)
        end
    end

    custom_legend(fig, colA, markB)

    display(fig)
    return fig
end