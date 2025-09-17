# Summarize one combo at a chosen loss (e.g., 0.6) for all 3 geometries
function summarize_combo_at_loss(name::Symbol; grid::Grid, S::Int=120, basal_frac::Float64=0.45,
                                 loss_pick::Float64=0.6,
                                 seeds_pool=1:5, seeds_mask=1:5, sim_seed::Int=1234,
                                 nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
                                 T_frac_on::Float64=0.02, τA::Float64=0.35, τocc::Float64=0.35)

    @assert haskey(COMBOS, name) "Unknown combo $(name)"
    axes = COMBOS[name]
    pars = bam_from_axes(; B_level=axes.B_level, M_level=axes.M_level, τA, τocc)
    bam = pars.bam

    # Scale T with grid size only if movement is on
    T = axes.M_level === :on ? max(8, round(Int, T_frac_on * grid.C)) : 8
    mp = MovementParams(; mode=(axes.M_level===:on ? :component : :none), T=T)

    keepfrac = 1 - loss_pick
    base_keep = trues(grid.C)

    out = Dict{Symbol,NamedTuple}()

    for hk in (:random, :clustered, :front)
        Δ  = Float64[]; AΔ = Float64[]; BΔ = Float64[]; MΔ = Float64[]; SΔ = Float64[]
        for ps in seeds_pool, ms in seeds_mask
            rng_pool = thread_rng(sim_seed, :pool, ps, hk)
            rng_mask = thread_rng(sim_seed, :mask, ms, hk)

            pool = build_pool_from_axes(rng_pool; S, basal_frac,
                                        A_level=axes.A_level, B_level=axes.B_level)
            A    = abiotic_matrix(pool, grid)
            keep = hk === :random    ? random_mask(rng_mask, grid.C, keepfrac) :
                   hk === :clustered ? clustered_mask(rng_mask, grid, keepfrac; nseeds=nseeds_cluster) :
                   hk === :front     ? frontlike_mask(rng_mask, grid, keepfrac; axis=front_axis, noise=front_noise) :
                   error("HL")

            shap,_,_ = shapley_per_species(pool, grid, A, base_keep, keep; bam=bam, mp=mp)
            cons = consumer_mask(pool)
            push!(Δ,  mean(getfield.(shap[cons], :dF)))
            push!(AΔ, mean(getfield.(shap[cons], :dA)))
            push!(BΔ, mean(getfield.(shap[cons], :dB)))
            push!(MΔ, mean(getfield.(shap[cons], :dM)))
            push!(SΔ, mean(getfield.(shap[cons], :synergy)))
        end
        dF  = mean(Δ);  dA = mean(AΔ);  dB = mean(BΔ);  dM = mean(MΔ);  syn = mean(SΔ)
        tot = max(abs(dF), 1e-12)
        sA  = abs(dA)/tot;  sB = abs(dB)/tot;  sM = abs(dM)/tot
        out[hk] = (; dF, dA, dB, dM, syn, sA, sB, sM)
    end

    return (; A_level=axes.A_level, B_level=axes.B_level, M_level=axes.M_level, T=T, loss=loss_pick, per_hl=out)
end

# Run ALL 18 combos (threaded over combos) and collect a flat vector of results
function summarize_all_combos(; nx::Int=60, ny::Int=60, S::Int=150, basal_frac::Float64=0.45,
                              loss_pick::Float64=0.6, T_frac_on::Float64=0.02,
                              seeds_pool=1:5, seeds_mask=1:5, sim_seed::Int=1234)
    grid = make_grid(nx, ny; seed=42)
    names = collect(keys(COMBOS))
    results = Vector{NamedTuple}(undef, length(names))
    Threads.@threads for i in eachindex(names)
        results[i] = summarize_combo_at_loss(names[i]; grid, S, basal_frac,
                                             loss_pick, T_frac_on,
                                             seeds_pool, seeds_mask, sim_seed)
    end
    return (grid=grid, names=names, results=results)
end

# Map shares (sA,sB,sM) to 2D in an equilateral triangle (ternary)
@inline function _ternary_xy(sA, sB, sM)
    # vertices: A=(0,0), B=(1,0), M=(0.5, √3/2)
    x = sB + 0.5*sM
    y = (√3/2)*sM
    return x, y
end

function plot_ternary_summary(sumres; title="BAM×HL summary at fixed loss")
    names   = sumres.names
    results = sumres.results

    # Prepare arrays per geometry
    geoms = (:random, :clustered, :front)
    labels = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")

    # color & marker encodings
    colA = Dict(:neutral=>RGBAf(0.55,0.55,0.55,1.0),
                :intermediate=>RGBAf(0.25,0.45,0.85,1.0),
                :divergent=>RGBAf(0.85,0.35,0.25,1.0))
    markB = Dict(:none=>:circle, :soft=>:utriangle, :strong=>:rect)

    # Geometry Sensitivity Index per combo
    GSI = Float64[]; combo_labels = Symbol[]

    # Figure
    fig = Figure(; size=(1250, 760))
    Label(fig[0, :], title; fontsize=16, padding=(0,0,10,0))

    # ── TOP: ternary, one axis per geometry ────────────────────────────────
    for (col, g) in enumerate(geoms)
        ax = Axis(fig[1, col], title=labels[g], xlabel="Biotic →", ylabel="Movement ↑")
        # draw triangle frame
        xs = [0.0, 1.0, 0.5, 0.0]; ys = [0.0, 0.0, √3/2, 0.0]
        lines!(ax, xs, ys)
        text!(ax, 0.02, -0.03, text="Abiotic", align=(:left,:bottom))
        text!(ax, 0.98, -0.03, text="Biotic",  align=(:right,:bottom))
        text!(ax, 0.50,  √3/2 + 0.03, text="Movement", align=(:center,:bottom))

        for (i, name) in enumerate(names)
            r = results[i]
            p = r.per_hl[g]
            sA, sB, sM = p.sA, p.sB, p.sM
            x,y = _ternary_xy(sA, sB, sM)
            Δmag = abs(p.dF)
            ms = 6 + 40*Δmag / (maximum([abs(results[j].per_hl[h].dF) for j in eachindex(results), h in geoms]) + 1e-12)
            color = colA[r.A_level]
            marker = markB[r.B_level]
            strokew = r.M_level === :on ? 1.5 : 0.5
            scatter!(ax, [x], [y]; markersize=ms, marker=marker, color=color, strokecolor=:black, strokewidth=strokew, transparency=false)
        end
        hidespines!(ax, :r, :t)
    end

    custom_legend(fig, colA, markB)

    # ── BOTTOM: geometry sensitivity index (bars) ───────────────────────────
    ax2 = Axis(
        fig[2, 1:3], 
        xlabel="Combo", ylabel="GSI = max(ΔF) - min(ΔF)", title="Geometry sensitivity",
        xticklabelrotation=π/6
    )
    xs = Float64[]; ys = Float64[]
    for i in eachindex(results)
        dFs = [results[i].per_hl[g].dF for g in geoms]
        push!(GSI, maximum(dFs) - minimum(dFs))
        push!(combo_labels, names[i])
    end
    ord = sortperm(GSI; rev=true)
    barpositions = 1:length(ord)
    barplot!(ax2, barpositions, GSI[ord])
    ax2.xticks = (barpositions, string.(combo_labels[ord]))
    # rotate!(ax2.xticklabelrotation, π/2)
    display(fig)
    return fig
end

function custom_legend(fig, colA, markB)
    leggrid = fig[1,4] = GridLayout()

    # Colors for A-levels
    for (i,(lev,col)) in enumerate(colA)
        ax = Axis(leggrid[i,1]; width=20, height=20)
        hidedecorations!(ax); hidespines!(ax)
        scatter!(ax, [0],[0]; color=col, markersize=12)
        Label(leggrid[i,2], String(lev); halign=:left)
    end

    offset = length(colA)

    # Markers for B-levels
    for (j,(lev,mk)) in enumerate(markB)
        ax = Axis(leggrid[offset+j,1]; width=20, height=20)
        hidedecorations!(ax); hidespines!(ax)
        scatter!(ax, [0],[0]; color=:gray, marker=mk, markersize=12)
        Label(leggrid[offset+j,2], String(lev); halign=:left)
    end

    offset += length(markB)

    # Stroke for M-levels
    ax = Axis(leggrid[offset+1,1]; width=20, height=20)
    hidedecorations!(ax); hidespines!(ax)
    scatter!(ax, [0],[0]; color=:white, strokecolor=:black, strokewidth=1.5, markersize=14)
    Label(leggrid[offset+1,2], "movement:on (thick outline)"; halign=:left)

    return leggrid
end

# --- quick diagnostics ---
# prey stats per consumer (mean±sd and distribution tail)
function prey_stats(pool)
    L = [length(pool.E[s]) for s in 1:pool.S if !pool.basal[s]]
    (; mean=mean(L), sd=std(L), p95=quantile(L,0.95), min=minimum(L), max=maximum(L))
end

# movement percolation: fraction of consumers with any cells failing the T rule
function movement_stats(pool, grid, A, keep; mp, bam)
    P, Bsup, _, M = assemble_BAM(pool, grid, A, keep; bam=bam, mp=mp)
    cons = .!pool.basal
    fails = [mean(@view M[s, :]) < 1.0 for s in 1:pool.S if cons[s]]
    (; frac_impacted=mean(fails))
end

# abiotic bite: mean As loss vs. B and M changes (baseline→mask)
function delta_ABM(pool, grid, A, base_keep, keep; bam, mp)
    P0, B0, _, M0 = assemble_BAM(pool, grid, A, base_keep; bam=bam, mp=mp)
    S0 = species_stats(pool, grid, A, base_keep, P0,B0,M0, bam)
    P1, B1, _, M1 = assemble_BAM(pool, grid, A, keep; bam=bam, mp=mp)
    S1 = species_stats(pool, grid, A, keep, P1,B1,M1, bam)
    cons = .!pool.basal
    dA = mean(S1.A[cons] .- S0.A[cons])
    dB = mean(S1.B[cons] .- S0.B[cons])
    dM = mean(S1.M[cons] .- S0.M[cons])
    (; dA, dB, dM)
end

# ---- climate axis (fit plane Climate ~ a*x + b*y, use (a,b) as axis) ----
function _climate_axis_proj(grid::Grid)
    C = grid.C
    X = hcat(ones(C), grid.xy[1,:], grid.xy[2,:])          # [1 x y]
    β = X \ grid.climate                                   # OLS β
    a, b = β[2], β[3]
    t = a .* grid.xy[1,:] .+ b .* grid.xy[2,:]             # projection
    # normalize to [0,1]
    t .-= minimum(t); t ./= (maximum(t) + eps())
    return t
end

"Low/high-tail boolean masks along the climate axis (default 20% each)."
function climate_tails(grid::Grid; q::Float64=0.20)
    t = _climate_axis_proj(grid)
    lo = quantile(t, q); hi = quantile(t, 1-q)
    low  = map(v->v ≤ lo, t)
    high = map(v->v ≥ hi, t)
    return (; low=BitVector(low), high=BitVector(high), t=t)
end

# ---- baseline indicators ----------------------------------------------------

"""
baseline_indicators(grid, S, basal_frac; A_level, B_level, M_level, τA, τocc, T_frac_on)
Return:
  D   :: Float64   # Climate divergence (tail-contrast of suitable density)
  R   :: NamedTuple(mean, p95)  # diet redundancy (#prey per consumer)
  Θ   :: Float64   # Movement stiffness = T / mean # climatically suitable cells
  LPRI:: NamedTuple(mean, min)  # Last-Prey Risk Index in threatened tail
Also returns structs you may want (pool, A, BAM & movement params).
"""
function baseline_indicators(grid::Grid, S::Int, basal_frac::Float64;
        A_level::Symbol, B_level::Symbol, M_level::Symbol,
        τA::Float64=0.5, τocc::Float64=0.35, T_frac_on::Float64=0.02,
        sim_seed::Int=1234, pool_seed::Int=1)

    rng_pool = MersenneTwister(hash((sim_seed, :pool, pool_seed)))
    pool     = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
    A        = abiotic_matrix(pool, grid)
    pars     = bam_from_axes(; B_level, M_level, τA, τocc)
    bam      = pars.bam
    T        = (M_level===:on) ? max(8, round(Int, T_frac_on*grid.C)) : 8
    mp       = MovementParams(; mode=(M_level===:on ? :component : :none), T=T)

    # Baseline assembly (keep everything)
    base_keep = trues(grid.C)
    P0, Bsup0, _, M0 = assemble_BAM(pool, grid, A, base_keep; bam=bam, mp=mp)

    # -- D: climate divergence (tail-contrast of suitable density)
    tails = climate_tails(grid)
    low, high = tails.low, tails.high
    # Suitability density across species (climate-only pass)
    Z = falses(pool.S, grid.C)
    @inbounds for s in 1:pool.S, c in 1:grid.C
        Z[s,c] = A[s,c] ≥ τA
    end
    dens = vec(mean(Z; dims=1))            # fraction of species suitable per cell
    sum_low  = sum(dens[low]);  sum_high = sum(dens[high])
    D = abs(sum_high - sum_low) / max(sum_low + sum_high, 1e-12)

    # -- R: diet redundancy (#prey per consumer) on the metaweb
    L = [length(pool.E[s]) for s in 1:pool.S if !pool.basal[s]]
    Rstats = (; mean=mean(L), p95=quantile(L, 0.95))

    # -- Θ: movement stiffness ~ how hard T is relative to suitable area
    # mean # of climatically suitable cells per consumer
    cons = .!pool.basal
    suitable_counts = [count(@view Z[s,:]) for s in 1:pool.S if cons[s]]
    Sbar = mean(suitable_counts)
    Θ = (M_level===:on) ? T / max(Sbar, 1e-9) : 0.0

    # -- LPRI: tail-specific prey support in threatened tail (baseline)
    # For each consumer: mean # prey present in low tail among cells where consumer is climatically suitable.
    function prey_count_in_tail(s::Int, tailmask::BitVector)
        cons_cells = findall(i -> Z[s,i] && tailmask[i], 1:grid.C)
        if isempty(cons_cells); return 0.0; end
        cnt = 0.0
        for i in cons_cells
            # count prey present in P0 at cell i
            pc = 0
            for q in pool.E[s]; pc += (P0[q,i] ? 1 : 0); end
            cnt += pc
        end
        return cnt / length(cons_cells)
    end
    prey_low  = [prey_count_in_tail(s, low)  for s in 1:pool.S if cons[s]]
    LPRI = (; mean=mean(prey_low), min=minimum(prey_low))

    return (; D, R=Rstats, Θ, LPRI, pool, A, bam, mp, tails, T)
end

"""
predict_geometry_ranking(D, R, Θ; movement_on::Bool)
Heuristic rule-of-thumb producing an ordering of expected damage:
returns a NamedTuple with symbols sorted from worst→best.
"""
function predict_geometry_ranking(D::Float64, R::NamedTuple, Θ::Float64; movement_on::Bool)
    # thresholds are conservative; tweak if you calibrate
    if movement_on && Θ ≥ 1.0
        # fragmentation dominates
        return (; worst=:random, middle=:clustered, best=:front, rationale="Θ≥1 ⇒ fragmentation → random worst")
    elseif D ≥ 0.2 && R.mean ≤ 3.0
        # climate tail + low redundancy ⇒ front risky
        return (; worst=:front, middle=:clustered, best=:random, rationale="large D & low R ⇒ last-prey in climate tail")
    else
        return (; worst=:clustered, middle=:front, best=:random, rationale="area loss dominates; minor geometry effect")
    end
end

# Mean ΔBSH (consumers) across seeds for a given geometry (threaded)
function sweep_dBSH_axes(; grid::Grid, S::Int=150, basal_frac::Float64=0.45,
        A_level::Symbol, B_level::Symbol, M_level::Symbol,
        hl_kind::Symbol, loss_fracs=0.2:0.1:0.8,
        seeds_pool=1:5, seeds_mask=1:5, sim_seed::Int=1234,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        τA::Float64=0.5, τocc::Float64=0.35, T_frac_on::Float64=0.02)

    pars = bam_from_axes(; B_level, M_level, τA, τocc)
    bam = pars.bam
    T   = (M_level===:on) ? max(8, round(Int, T_frac_on*grid.C)) : 8
    mp  = MovementParams(; mode=(M_level===:on ? :component : :none), T=T)

    xs = collect(loss_fracs); y = similar(xs, Float64)
    base_keep = trues(grid.C)

    Threads.@threads for k in eachindex(xs)
        f = xs[k]; keepfrac = 1 - f
        vals = Float64[]
        for ps in seeds_pool, ms in seeds_mask
            rng_pool = thread_rng(sim_seed, :pool, ps, k)
            rng_mask = thread_rng(sim_seed, :mask, ms, k, hl_kind)
            pool = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
            A    = abiotic_matrix(pool, grid)
            keep = hl_kind === :random    ? random_mask(rng_mask, grid.C, keepfrac) :
                   hl_kind === :clustered ? clustered_mask(rng_mask, grid, keepfrac; nseeds=nseeds_cluster) :
                   hl_kind === :front     ? frontlike_mask(rng_mask, grid, keepfrac; axis=front_axis, noise=front_noise) :
                   error("hl_kind")
            # baseline vs scenario BSH
            P0,B0,_,M0 = assemble_BAM(pool, grid, A, base_keep; bam=bam, mp=mp)
            S0 = species_stats(pool, grid, A, base_keep, P0,B0,M0, bam)
            P1,B1,_,M1 = assemble_BAM(pool, grid, A, keep;       bam=bam, mp=mp)
            S1 = species_stats(pool, grid, A, keep, P1,B1,M1, bam)
            cons = consumer_mask(pool)
            push!(vals, mean(S1.BSH[cons] .- S0.BSH[cons]))
        end
        y[k] = mean(vals)
    end
    return (; x=xs, y=y, meta=Dict(:hl=>hl_kind, :A=>A_level, :B=>B_level, :M=>M_level))
end

# Fragmentation diagnostic: fraction of suitable cells meeting T (M=1)
function frag_diagnostic_axes(; grid::Grid, S::Int=150, basal_frac::Float64=0.45,
        A_level::Symbol, B_level::Symbol, M_level::Symbol, hl_kind::Symbol,
        loss_fracs=0.2:0.1:0.8, seeds_pool=1:4, seeds_mask=1:4, sim_seed::Int=1234,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        τA::Float64=0.5, τocc::Float64=0.35, T_frac_on::Float64=0.02)

    pars = bam_from_axes(; B_level, M_level, τA, τocc)
    bam = pars.bam
    T   = (M_level===:on) ? max(8, round(Int, T_frac_on*grid.C)) : 8
    mp  = MovementParams(; mode=(M_level===:on ? :component : :none), T=T)

    xs = collect(loss_fracs); y = similar(xs, Float64)
    base_keep = trues(grid.C)

    Threads.@threads for k in eachindex(xs)
        f = xs[k]; keepfrac = 1 - f
        vals = Float64[]
        for ps in seeds_pool, ms in seeds_mask
            rng_pool = thread_rng(sim_seed, :pool, ps, k)
            rng_mask = thread_rng(sim_seed, :mask, ms, k, hl_kind)
            pool = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
            A    = abiotic_matrix(pool, grid)
            keep = hl_kind === :random    ? random_mask(rng_mask, grid.C, keepfrac) :
                   hl_kind === :clustered ? clustered_mask(rng_mask, grid, keepfrac; nseeds=nseeds_cluster) :
                   hl_kind === :front     ? frontlike_mask(rng_mask, grid, keepfrac; axis=front_axis, noise=front_noise) :
                   error("hl_kind")
            # compute M on the mask only (conditional on suitability)
            P1,B1,_,M1 = assemble_BAM(pool, grid, A, keep; bam=bam, mp=mp)
            cons = consumer_mask(pool)
            # share of suitable cells with M=1
            share = Float64[]
            for s in 1:pool.S
                cons[s] || continue
                idx = findall(i -> keep[i] && A[s,i] ≥ τA, 1:grid.C)
                if isempty(idx); continue; end
                push!(share, mean(@view M1[s, idx]))
            end
            if isempty(share); push!(vals, NaN); else push!(vals, mean(share)); end
        end
        y[k] = mean(skipmissing(share for share in [vals...]))  # quick safe mean
    end
    return (; x=xs, y=y, T=T)
end

function plot_dashboard(; nx::Int=60, ny::Int=60, S::Int=150, basal_frac::Float64=0.45,
        A_level::Symbol=:intermediate, B_level::Symbol=:soft, M_level::Symbol=:on,
        loss_fracs=0.2:0.1:0.8, show_loss_for_fingerprints::Float64=0.6,
        seeds_pool=1:6, seeds_mask=1:6, sim_seed::Int=1234,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        τA::Float64=0.5, τocc::Float64=0.35, T_frac_on::Float64=0.02)

    grid = make_grid(nx, ny; seed=42)

    # --- ΔBSH curves (impact) ---
    dR = sweep_dBSH_axes(; grid,S,basal_frac,A_level,B_level,M_level,hl_kind=:random,
                         loss_fracs,seeds_pool,seeds_mask,sim_seed,nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)
    dC = sweep_dBSH_axes(; grid,S,basal_frac,A_level,B_level,M_level,hl_kind=:clustered,
                         loss_fracs,seeds_pool,seeds_mask,sim_seed,nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)
    dF = sweep_dBSH_axes(; grid,S,basal_frac,A_level,B_level,M_level,hl_kind=:front,
                         loss_fracs,seeds_pool,seeds_mask,sim_seed,nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)

    # --- fingerprints: show in a separate window (unchanged) ---
    raw  = summarize_all_combos(; nx, ny, S, basal_frac, loss_pick=show_loss_for_fingerprints,
                                          seeds_pool, seeds_mask, sim_seed)
    elas = elasticity_summarize_all(; nx, ny, S, basal_frac, loss_pick=show_loss_for_fingerprints,
                                              seeds_pool, seeds_mask, sim_seed)
    subfig = compare_raw_vs_elasticity(raw; sum_elast=elas,
        title=@sprintf("Raw ΔF vs Elasticity shares (loss=%.2f)", show_loss_for_fingerprints))
    display(subfig)

    # --- fragmentation diagnostic ---
    fR = frag_diagnostic_axes(; grid,S,basal_frac,A_level,B_level,M_level,hl_kind=:random,
                               loss_fracs,seeds_pool,seeds_mask,sim_seed,nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)
    fC = frag_diagnostic_axes(; grid,S,basal_frac,A_level,B_level,M_level,hl_kind=:clustered,
                               loss_fracs,seeds_pool,seeds_mask,sim_seed,nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)
    fF = frag_diagnostic_axes(; grid,S,basal_frac,A_level,B_level,M_level,hl_kind=:front,
                               loss_fracs,seeds_pool,seeds_mask,sim_seed,nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)

    # --- baseline indicators + prediction ---
    base = baseline_indicators(grid, S, basal_frac; A_level,B_level,M_level,τA,τocc,T_frac_on)
    pred = predict_geometry_ranking(base.D, base.R, base.Θ; movement_on=(M_level===:on))

    # consistent colors (Makie wants RGBf / RGBAf)
    colors = Dict(:random   => RGBf(0.25,0.45,0.85),
                  :clustered=> RGBf(0.95,0.60,0.15),
                  :front    => RGBf(0.20,0.65,0.35))

    # ─────────────────────────────  FIGURE  ────────────────────────────────
    fig = Figure(; size=(1400, 1000))
    Label(fig[0, :], "BAM Dashboard — A=$(A_level), B=$(B_level), M=$(M_level)";
          fontsize=18, padding=(0,0,8,0))

    # TOP-LEFT: ΔBSH panel with legend in its own column
    gl_top_left = fig[1, 1:2] = GridLayout()
    ax1  = Axis(gl_top_left[1,1], xlabel="Area lost (fraction)",
                ylabel="ΔBSH (mean consumers)", title="ΔBSH vs loss")
    hR = lines!(ax1, dR.x, dR.y; color=colors[:random],    linewidth=2)
    hC = lines!(ax1, dC.x, dC.y; color=colors[:clustered], linewidth=2)
    hF = lines!(ax1, dF.x, dF.y; color=colors[:front],     linewidth=2)
    leg1 = Legend(gl_top_left[1,2], [hR,hC,hF], ["Random","Clustered","Front-like"];
                  framevisible=false, padding=(2,2,2,2))
    rowgap!(gl_top_left, 10); colgap!(gl_top_left, 10)
    gl_top_left[1,2].width = 120  # reserve space for legend

    # TOP-RIGHT: Fragmentation panel with its own legend column
    gl_top_right = fig[1, 3] = GridLayout()
    ax2 = Axis(gl_top_right[1,1], xlabel="Area lost (fraction)",
               ylabel="Share of suitable cells with M=1", title="Fragmentation diagnostic")
    gR = lines!(ax2, fR.x, fR.y; color=colors[:random],    linewidth=2)
    gC = lines!(ax2, fC.x, fC.y; color=colors[:clustered], linewidth=2)
    gF = lines!(ax2, fF.x, fF.y; color=colors[:front],     linewidth=2)
    hlines!(ax2, [0.5]; linestyle=:dash, color=RGBAf(0,0,0,0.4))
    leg2 = Legend(gl_top_right[1,2], [gR,gC,gF], ["Random","Clustered","Front-like"];
                  framevisible=false, padding=(2,2,2,2))
    gl_top_right[1,2].width = 120

    # BOTTOM-RIGHT: Baseline indices box (bigger so text never clips)
    ax3 = Axis(fig[2,3], title="Baseline indices & prediction"); hidedecorations!(ax3)
    y = 1.0; dy = 0.09
    text!(ax3, 0.0, y, text=@sprintf("D (climate divergence): %.2f", base.D), align=(:left,:top)); y -= dy
    text!(ax3, 0.0, y, text=@sprintf("R (diet redundancy): mean=%.2f, p95=%.2f", base.R.mean, base.R.p95), align=(:left,:top)); y -= dy
    text!(ax3, 0.0, y, text=@sprintf("Θ (movement stiffness): %.2f   (T=%d)", base.Θ, base.T), align=(:left,:top)); y -= dy
    text!(ax3, 0.0, y, text=@sprintf("LPRI (tail prey): mean=%.2f, min=%.2f", base.LPRI.mean, base.LPRI.min), align=(:left,:top)); y -= 1.2dy
    text!(ax3, 0.0, y, text="Predicted ranking (worst → best):", align=(:left,:top)); y -= dy
    text!(ax3, 0.0, y, text=string(pred), align=(:left,:top))

    display(fig)
    return fig
end

# Example dashboard for a “high-contrast” corner
# _ = plot_dashboard(
#     ; nx=60, ny=60, S=150,
#     A_level=:divergent, B_level=:strong, M_level=:on,
#     loss_fracs=0.2:0.1:0.8, show_loss_for_fingerprints=0.6,
#     T_frac_on=0.03
# )

"""
audit_phase_rule(loss_pick=0.6) — compute baseline indices for each combo,
observe geometry ranking from ΔBSH at loss_pick, and compare to prediction.
Returns a vector of NamedTuples (one per combo).
"""
function audit_phase_rule(; nx::Int=60, ny::Int=60, S::Int=150, basal_frac::Float64=0.45,
        loss_pick::Float64=0.6, seeds_pool=1:5, seeds_mask=1:5, sim_seed::Int=1234,
        τA::Float64=0.5, τocc::Float64=0.35, T_frac_on::Float64=0.02)

    grid = make_grid(nx, ny; seed=42)
    names = collect(keys(COMBOS))
    out = NamedTuple[]
    for nm in names
        axes = COMBOS[nm]
        # Baseline indicators
        base = baseline_indicators(grid, S, basal_frac; A_level=axes.A_level, B_level=axes.B_level,
                                   M_level=axes.M_level, τA, τocc, T_frac_on)
        pred = predict_geometry_ranking(base.D, base.R, base.Θ; movement_on=(axes.M_level===:on))
        # Observed ΔBSH at loss_pick
        for g in (:random, :clustered, :front)
            nothing
        end
        f = loss_pick
        d = Dict{Symbol,Float64}()
        for g in (:random, :clustered, :front)
            res = sweep_dBSH_axes(; grid,S,basal_frac,A_level=axes.A_level,B_level=axes.B_level,M_level=axes.M_level,
                                  hl_kind=g, loss_fracs=[f], seeds_pool, seeds_mask, sim_seed)
            d[g] = res.y[1]
        end
        # worst→best (most negative ΔBSH is “worst”)
        obs_order = sort(collect(keys(d)); by=x->d[x])
        push!(out, (; name=nm, A=axes.A_level, B=axes.B_level, M=axes.M_level,
                    D=base.D, R_mean=base.R.mean, R_p95=base.R.p95, Θ=base.Θ,
                    LPRI_mean=base.LPRI.mean, LPRI_min=base.LPRI.min,
                    pred_worst=pred.worst, pred_best=pred.best, obs_order=obs_order,
                    d_random=d[:random], d_clustered=d[:clustered], d_front=d[:front]))
    end
    return out
end

# A = audit_phase_rule(; loss_pick=0.6)
