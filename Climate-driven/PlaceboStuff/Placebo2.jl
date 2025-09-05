# --- 0.1 Climate + prey-support scores at the baseline (before loss) ---
function cell_scores(Z0::BitMatrix, P0::BitMatrix, pool::SpeciesPool)
    cons = .!pool.basal
    C    = size(Z0,2)

    # Climate score per cell: fraction of consumers climate-suitable
    A_cell = vec(mean(view(Z0, cons, :); dims=1))

    # Prey-support per cell: for each consumer, fraction of its prey present; then mean across consumers
    preyfrac = zeros(sum(cons), C)
    row = 1
    for s in findall(cons)
        prey = pool.E[s]; diet = max(1, length(prey))
        if diet == 0
            row += 1; continue
        end
        @inbounds for c in 1:C
            k = 0
            @inbounds for q in prey
                k += P0[q,c] ? 1 : 0
            end
            preyfrac[row, c] = k / diet
        end
        row += 1
    end
    prey_cell = vec(mean(preyfrac; dims=1))
    return A_cell, prey_cell
end

# --- 0.2 A-matched hotspot: within climate deciles, keep highest-prey cells ---
function a_matched_hotspot_mask(A_cell::Vector{Float64}, prey_cell::Vector{Float64};
                                keep_frac::Float64, nbins::Int=10, seed::Int=1)
    rng = MersenneTwister(seed)
    C = length(A_cell); keep = falses(C)
    qs = [quantile(A_cell, i/nbins) for i in 0:nbins]
    for b in 1:nbins
        idx = findall(c -> (A_cell[c] ≥ qs[b]) && (A_cell[c] ≤ qs[b+1]), 1:C)
        isempty(idx) && continue
        kkeep = round(Int, keep_frac * length(idx))
        ord   = sortperm(prey_cell[idx]; rev=true)  # highest-prey first
        keep[idx[ord[1:kkeep]]] .= true
    end
    return keep
end

# --- 0.3 Degree-preserving rewiring placebo (consumers keep diet size; prey are random basal) ---
function rewire_metaweb(pool::SpeciesPool; seed::Int=42)
    rng = MersenneTwister(seed)
    nb = count(pool.basal); S = pool.S
    newE = [Int[] for _ in 1:S]
    for s in 1:nb
        newE[s] = Int[]
    end
    for s in nb+1:S
        k = length(pool.E[s])
        newE[s] = rand(rng, 1:nb, k)
    end
    SpeciesPool(pool.S, copy(pool.masses), copy(pool.basal),
                copy(pool.mu), copy(pool.b), newE)
end

"""
run_excess_with_controls(pool, grid; τ, losses, hotspot=:A_matched, seed=101)

Returns (data, losses) where data is a Dict of scenario ⇒ Vector{NamedTuple}
with fields (dB, dA_only, dInt_only, synergy) from decomp_at_mask_fraction.
Scenarios: :random, :cluster, :hotspot and their *_rw rewired placebos.
"""
function run_excess_with_controls(pool::SpeciesPool, grid::Grid;
    τ::Float64, losses::Vector{Float64}, hotspot::Symbol=:A_matched, seed::Int=101)

    # Baseline for ORIGINAL pool
    Z0 = climate_pass(pool, grid; τ=τ)
    P0 = assemble(Z0, pool)
    A_cell, prey_cell = cell_scores(Z0, P0, pool)

    # Baseline for REWIRED placebo (build its own hotspot fairly)
    pool_rw = rewire_metaweb(pool; seed=seed)
    Z0_rw   = climate_pass(pool_rw, grid; τ=τ)
    P0_rw   = assemble(Z0_rw, pool_rw)
    A_cell_rw, prey_cell_rw = cell_scores(Z0_rw, P0_rw, pool_rw)

    # helpers to build keepmasks
    function hotspot_keep(keep_frac; use_rw::Bool=false)
        if hotspot == :A_matched
            return use_rw ?
                a_matched_hotspot_mask(A_cell_rw, prey_cell_rw; keep_frac=keep_frac, nbins=10, seed=seed) :
                a_matched_hotspot_mask(A_cell,    prey_cell;    keep_frac=keep_frac, nbins=10, seed=seed)
        elseif hotspot == :prey_only
            C = length(A_cell); kkeep = round(Int, keep_frac*C)
            if use_rw
                ord = sortperm(prey_cell_rw; rev=true)
                keep = falses(C); keep[ord[1:kkeep]] .= true; keep
            else
                ord = sortperm(prey_cell; rev=true)
                keep = falses(C); keep[ord[1:kkeep]] .= true; keep
            end
        else
            error("hotspot must be :A_matched or :prey_only")
        end
    end

    data = Dict{Symbol,Vector{NamedTuple}}(
        :random      => NamedTuple[],
        :cluster     => NamedTuple[],
        :hotspot     => NamedTuple[],
        :random_rw   => NamedTuple[],
        :cluster_rw  => NamedTuple[],
        :hotspot_rw  => NamedTuple[]
    )

    for loss in losses
        keep = 1 - loss
        kmR  = random_mask(grid.C, keep; seed=seed)
        kmC  = clustered_mask(grid, keep; nseeds=1, seed=seed)
        kmH  = hotspot_keep(keep)

        # original
        push!(data[:random],  decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmR))
        push!(data[:cluster], decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmC))
        push!(data[:hotspot], decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmH))

        # rewired placebo (re-use kmR/kmC; build fair kmH_rw)
        push!(data[:random_rw],  decomp_at_mask_fraction(pool_rw, grid; τ=τ, keepmask=kmR))
        push!(data[:cluster_rw], decomp_at_mask_fraction(pool_rw, grid; τ=τ, keepmask=kmC))
        kmH_rw = hotspot_keep(keep; use_rw=true)
        push!(data[:hotspot_rw], decomp_at_mask_fraction(pool_rw, grid; τ=τ, keepmask=kmH_rw))
    end
    return (; data, losses)
end

# --- Excess decomposition & interaction share for a scenario (vs random) ---
function plot_excess_and_share(results; scen=:hotspot, label="A-matched hotspot")
    losses = results.losses
    R = results.data[:random];  S = results.data[scen]

    ex    = [S[i].dB        - R[i].dB        for i in eachindex(losses)]
    exA   = [S[i].dA_only   - R[i].dA_only   for i in eachindex(losses)]
    exInt = [S[i].dInt_only - R[i].dInt_only for i in eachindex(losses)]
    exSyn = [S[i].synergy   - R[i].synergy   for i in eachindex(losses)]
    share = [abs(ex[i]) < 1e-12 ? NaN : exInt[i] / ex[i] for i in eachindex(losses)]

    fig = Figure(; size=(980,360))

    ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Excess ΔBSH (scenario − random)",
               title=label)
    lines!(ax1, losses, exA,   label="Climate-only")
    lines!(ax1, losses, exInt, label="Interaction-only")
    lines!(ax1, losses, exSyn, label="Synergy")
    hlines!(ax1, [0], color=(:gray,0.5), linestyle=:dash); axislegend(ax1, position=:lt)

    ax2 = Axis(fig[1,2], xlabel="Area lost (fraction)", ylabel="Interaction share of excess")
    lines!(ax2, losses, share)
    hlines!(ax2, [0,1], color=(:gray,0.4), linestyle=:dash)
    display(fig)

    return (; excess=ex, exA, exInt, exSyn, share)
end

# --- Placebo check (rewiring) ---
function plot_placebo(results; scen=:hotspot, label="Placebo (rewired)")
    losses = results.losses
    R  = results.data[:random];    S  = results.data[scen]
    Rw = results.data[:random_rw]; Sw = results.data[Symbol(scen, "_rw")]

    ex   = [S[i].dB  - R[i].dB  for i in eachindex(losses)]
    ex_w = [Sw[i].dB - Rw[i].dB for i in eachindex(losses)]

    fig = Figure(; size=(620,320))
    ax  = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Excess ΔBSH vs random",
               title=label)
    lines!(ax, losses, ex,  label="with interactions")
    lines!(ax, losses, ex_w, label="rewired placebo", linestyle=:dash)
    hlines!(ax, [0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax, position=:lt)
    display(fig)
end

# Four regimes (low/high redundancy × low/high synchrony)
function RS_regimes_for_max()
    Dict(
        :lowR_lowS  => (sigma=0.20, density=0.12, pmax=0.70, niche_mode=:uniform, mu_basal_sd=0.08),
        :lowR_highS => (sigma=0.20, density=0.12, pmax=0.70, niche_mode=:bimodal, mu_basal_sd=0.03),  # ★ strongest
        :highR_lowS => (sigma=0.45, density=0.35, pmax=0.90, niche_mode=:uniform, mu_basal_sd=0.08),
        :highR_highS=> (sigma=0.45, density=0.35, pmax=0.90, niche_mode=:bimodal, mu_basal_sd=0.03)
    )
end

# ----- Fraction-basis ΔBSH decomposition (mean over consumers) -----

# tiny helper: turn NaNs into 0 so species with no suitable cells don't explode
nan_to_zero!(v) = (for i in eachindex(v); if isnan(v[i]); v[i] = 0.0; end; end; v)

"""
decomp_at_mask_fraction(pool, grid; τ, keepmask)

Returns a NamedTuple with fields:
  dB, dA_only, dInt_only, synergy
Each is the mean across **consumers**.
"""
function decomp_at_mask_fraction(pool::SpeciesPool, grid::Grid; τ::Float64, keepmask::BitVector)
    # before / after climate masks and assembled presence
    Z0 = climate_pass(pool, grid; τ=τ);   P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask);        P1 = assemble(Z1, pool)

    cons = .!pool.basal

    # A = climate-suitable fraction per sp (fraction of all cells)
    A0 = climate_area_per_species(Z0)[cons]
    A1 = climate_area_per_species(Z1)[cons]

    # p = conditional occupancy given climate suitability (0..1; can be NaN if A==0)
    p0 = bsh1_cond_per_species(P0, Z0, pool)[cons]
    p1 = bsh1_cond_per_species(P1, Z1, pool)[cons]
    nan_to_zero!(p0); nan_to_zero!(p1)

    # B = A * p
    B0 = A0 .* p0
    B1 = A1 .* p1

    # midpoint (Shapley-style) split
    dA    = A1 .- A0
    dp    = p1 .- p0
    Abar  = 0.5 .* (A0 .+ A1)
    pbar  = 0.5 .* (p0 .+ p1)

    dB        = B1 .- B0
    dA_only   = pbar .* dA
    dInt_only = Abar .* dp
    synergy   = dB   .- dA_only .- dInt_only   # = 0 numerically (up to fp rounding)

    mean_cons(x) = mean(x)   # already consumers only

    return (
        dB        = mean_cons(dB),
        dA_only   = mean_cons(dA_only),
        dInt_only = mean_cons(dInt_only),
        synergy   = mean_cons(synergy)
    )
end

# --- One-shot driver: scans regimes, selects the strongest, then plots everything ---
begin
    grid    = make_grid(60,60; seed=11)
    S       = 220
    τ       = 0.58                 # tighter climate threshold accentuates p-effects
    losses  = collect(0.15:0.05:0.80)
    regimes = RS_regimes_for_max()

    # Build/score each regime and keep the best (max interaction share at 50% loss)
    best = nothing; best_share = -Inf; best_name = ""
    results_by_regime = Dict{Symbol,Any}()

    for (name, kw) in regimes
        pool = build_pool(S;
            basal_frac=0.35, seed=1,
            sigma=kw[:sigma], density=kw[:density], pmax=kw[:pmax],
            niche_mode=kw[:niche_mode], mu_basal_centers=(0.25,0.75), mu_basal_sd=kw[:mu_basal_sd],
            b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04)

        res = run_excess_with_controls(pool, grid; τ=τ, losses=losses, hotspot=:A_matched, seed=101)
        results_by_regime[name] = res

        # interaction share at ~50% loss for A-matched hotspot
        i50 = argmin(abs.(losses .- 0.50))
        R = res.data[:random]; H = res.data[:hotspot]
        ex   = H[i50].dB        - R[i50].dB
        exInt= H[i50].dInt_only - R[i50].dInt_only
        share = abs(ex) < 1e-12 ? -Inf : exInt / ex
        if share > best_share
            best_share = share
            best = res; best_name = String(name)
        end
    end

    println("Best regime at 50% loss (interaction share of excess) → ",
            best_name, " ; share = ", round(best_share, digits=3))

    # Plot decomposition & share for the best regime
    sums = plot_excess_and_share(best; scen=:hotspot,
        label="A-matched hotspot (regime: $(best_name))")
    # plot_total_and_relative(losses, ex; title="A-matched hotspot (lowR–highS)")
    # Placebo: rewiring should suppress excess if interactions drive it
    plot_placebo(best; scen=:hotspot, label="Placebo (rewired) — $(best_name)")

    # (Optional) also check 'clustered' quickly:
    plot_excess_and_share(best; scen=:cluster, label="Clustered (same regime)")
    # plot_total_and_relative(losses, ex; title="Clustered (lowR–highS)")
    plot_placebo(best; scen=:cluster, label="Placebo (clustered) — $(best_name)")
end

"Left: excess components; Right: interaction share."
function plot_excess_and_share(losses, ex; title="")
    fig = Figure(; size=(950,380))

    ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)",
               ylabel="Excess damage (scenario − random)", title=title)
    lines!(ax1, losses, ex.ED_clim, label="Climate-only")
    lines!(ax1, losses, ex.ED_int,  label="Interaction-only")
    lines!(ax1, losses, ex.ED_syn,  label="Synergy")
    hlines!(ax1, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax1, position=:lt)

    ax2 = Axis(fig[1,2], xlabel="Area lost (fraction)",
               ylabel="Interaction share of excess")
    lines!(ax2, losses, ex.share)
    hlines!(ax2, [0.0, 1.0], color=(:gray,0.4), linestyle=:dash)

    display(fig); fig
end

"Left: total excess; Right: relative excess vs random (unit-free)."
function plot_total_and_relative(losses, ex; title="")
    fig = Figure(; size=(900,360))

    ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)",
               ylabel="Excess damage (total)", title=title)
    lines!(ax1, losses, ex.ED_tot)
    hlines!(ax1, [0.0], color=(:gray,0.5), linestyle=:dash)

    ax2 = Axis(fig[1,2], xlabel="Area lost (fraction)",
               ylabel="Relative excess (vs random)")
    lines!(ax2, losses, ex.RED_tot)
    hlines!(ax2, [0.0], color=(:gray,0.5), linestyle=:dash)

    display(fig); fig
end

"Plot excess total damage of scenario vs random, and its rewired placebo."
function plot_placebo_excess(res; scen::Symbol=:hotspot, title="")
    dR = collect_decomp_series(res.data[:random])
    dS = collect_decomp_series(res.data[scen])
    exS = excess_from_decomps(dR, dS)
    ED_s = exS.ED_tot

    # allow either :placebo or :rewired container
    plac = haskey(res, :placebo) ? res.placebo :
           (haskey(res, :rewired) ? res.rewired :
            error("This result does not contain placebo/rewired outputs."))

    dP  = collect_decomp_series(plac[scen])
    ED_p = excess_from_decomps(dR, dP).ED_tot

    L = hasproperty(res, :losses) ? res.losses : error("res.losses not found")

    fig = Figure(; size=(760,320))
    ax  = Axis(fig[1,1], xlabel="Area lost (fraction)",
               ylabel="Excess damage vs random", title=title)
    lines!(ax, L, ED_s, label="with interactions")
    lines!(ax, L, ED_p, linestyle=:dash, label="rewired placebo")
    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax, position=:lt)
    display(fig); fig
end

begin
    grid    = make_grid(60,60; seed=11)
    S       = 220
    τ       = 0.58
    losses  = collect(0.15:0.05:0.80)
    regimes = RS_regimes_for_max()

    best_res = nothing
    best_name = ""
    best_share = -Inf

    for (name, kw) in regimes
        pool = build_pool(S;
            basal_frac=0.35, seed=1,
            sigma=kw[:sigma], density=kw[:density], pmax=kw[:pmax],
            niche_mode=kw[:niche_mode], mu_basal_centers=(0.25,0.75),
            mu_basal_sd=kw[:mu_basal_sd],
            b0_basal=0.08, bspread_basal=0.02,
            b0_cons=0.12,  bspread_cons=0.04)

        res = run_excess_with_controls(pool, grid; τ=τ, losses=losses, hotspot=:A_matched, seed=101)

        # compute interaction-share of excess for HOTSPOT at ~50% loss
        dR = collect_decomp_series(res.data[:random])
        dH = collect_decomp_series(res.data[:hotspot])
        exH = excess_from_decomps(dR, dH)
        i50 = argmin(abs.(losses .- 0.50))
        share50 = exH.share[i50]

        if share50 > best_share
            best_share = share50
            best_name  = String(name)
            best_res   = res
        end
    end

    println("Best regime at 50% loss (interaction share of excess) → ",
            best_name, " ; share = ", round(best_share, digits=3))

    # ----- Plot A-matched hotspot for BEST regime -----
    dR  = collect_decomp_series(best_res.data[:random])
    dH  = collect_decomp_series(best_res.data[:hotspot])
    exH = excess_from_decomps(dR, dH)

    plot_excess_and_share(best_res.losses, exH;
        title="A-matched hotspot (regime: $best_name)")
    plot_total_and_relative(best_res.losses, exH;
        title="A-matched hotspot (regime: $best_name)")
    plot_placebo_excess(best_res; scen=:hotspot,
        title="Placebo (rewired) — $best_name")

    # ----- Also show CLUSTERED for the SAME regime -----
    dC  = collect_decomp_series(best_res.data[:cluster])
    exC = excess_from_decomps(dR, dC)

    plot_excess_and_share(best_res.losses, exC;
        title="Clustered (regime: $best_name)")
    plot_total_and_relative(best_res.losses, exC;
        title="Clustered (regime: $best_name)")
    plot_placebo_excess(best_res; scen=:cluster,
        title="Placebo (clustered) — $best_name")
end

