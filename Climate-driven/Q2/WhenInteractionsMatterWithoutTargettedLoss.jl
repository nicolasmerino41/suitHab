# --------------------------
# Viability metric (Q2 version)
# --------------------------
"""
boolean viability per species: 1{A ≥ A_min} * 1{p ≥ p_min}
A, p are per-species (consumers); returns BitVector
"""
function viability_bits(A::AbstractVector, p::AbstractVector; A_min=0.15, p_min=0.15)
    @assert length(A) == length(p)
    v = falses(length(A))
    @inbounds for i in eachindex(A)
        if !isnan(A[i]) && !isnan(p[i]) && A[i] ≥ A_min && p[i] ≥ p_min
            v[i] = true
        end
    end
    return v
end

# mean viability over consumers
mean_viability(A, p; A_min=0.15, p_min=0.15) = mean(viability_bits(A, p; A_min=A_min, p_min=p_min))

"""
viability-based decomposition at a given keepmask.
Returns (clim, inter, syn) where: ΔU = clim + inter + syn
"""
function decomp_viability_at_mask(pool::SpeciesPool, grid::Grid;
                                  τ::Float64, keepmask::BitVector,
                                  A_min::Float64=0.15, p_min::Float64=0.15)

    # assemble before/after
    Z0 = climate_pass(pool, grid; τ=τ);  P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask);       P1 = assemble(Z1, pool)

    cons = .!pool.basal

    # per-consumer climate area fraction
    A0 = climate_area_per_species(Z0)[cons]
    A1 = climate_area_per_species(Z1)[cons]

    # per-consumer conditional biotic support (given climate)
    p0 = bsh1_cond_per_species(P0, Z0, pool)[cons]
    p1 = bsh1_cond_per_species(P1, Z1, pool)[cons]

    # guard NaNs → 0 (no climate or undefined assembly ⇒ not viable)
    fix!(x) = (for i in eachindex(x); if isnan(x[i]); x[i]=0.0; end; end; x)
    fix!(p0); fix!(p1)

    # four “states”: (A0,p0), (A1,p0), (A0,p1), (A1,p1)
    V00 = mean(viability_bits(A0, p0; A_min=A_min, p_min=p_min))
    V10 = mean(viability_bits(A1, p0; A_min=A_min, p_min=p_min))
    V01 = mean(viability_bits(A0, p1; A_min=A_min, p_min=p_min))
    V11 = mean(viability_bits(A1, p1; A_min=A_min, p_min=p_min))

    clim = V10 - V00                 # climate-only
    inter = V01 - V00                # interaction-only
    syn  = V11 - V10 - V01 + V00     # genuine synergy

    return (clim=clim, inter=inter, syn=syn, U0=V00, U1=V11)
end

"Extreme regimes to elicit interaction effects without prey-targeting."
function VI_regimes()
    Dict(
        :lowR_lowS  => (sigma=0.20, density=0.10, pmax=0.65, niche_mode=:uniform, mu_basal_sd=0.09),
        :lowR_highS => (sigma=0.20, density=0.10, pmax=0.65, niche_mode=:bimodal, mu_basal_sd=0.03),
        :highR_lowS => (sigma=0.45, density=0.35, pmax=0.90, niche_mode=:uniform, mu_basal_sd=0.09),
        :highR_highS=> (sigma=0.45, density=0.35, pmax=0.90, niche_mode=:bimodal, mu_basal_sd=0.03),
    )
end

"""
Run viability decomposition across regimes & loss fractions.
Returns a Vector of NamedTuples (easy to DataFrame) with:
(:regime,:scenario,:loss,:clim,:inter,:syn,:U0,:U1)
"""
function sweep_viability(
    ; grid::Grid, τ::Float64=0.55, S::Int=200, basal_frac=0.35,
    keep_fracs=collect(0.20:0.05:0.80), regimes=VI_regimes(),
    include_hotspot_cons::Bool=true, include_hotspot_prey::Bool=true,
    A_min=0.15, p_min=0.15,
    seed_pool=1, seed_masks=101
)

    rows = NamedTuple[]
    for (rname, kw) in regimes
        pool = build_pool(S;
            basal_frac=basal_frac, seed=seed_pool,
            sigma=kw[:sigma], density=kw[:density], pmax=kw[:pmax],
            niche_mode=kw[:niche_mode], mu_basal_centers=(0.25,0.75), mu_basal_sd=kw[:mu_basal_sd],
            b0_basal=0.08, bspread_basal=0.02, b0_cons=0.12, bspread_cons=0.04
        )
        for keep in keep_fracs
            kmR = random_mask(grid.C, keep; seed=seed_masks)
            kmC = clustered_mask(grid, keep; nseeds=1, seed=seed_masks)
            outR = decomp_viability_at_mask(pool, grid; τ=τ, keepmask=kmR, A_min=A_min, p_min=p_min)
            outC = decomp_viability_at_mask(pool, grid; τ=τ, keepmask=kmC, A_min=A_min, p_min=p_min)
            push!(rows, (; regime=String(rname), scenario="random",    loss=1-keep, outR...))
            push!(rows, (; regime=String(rname), scenario="clustered", loss=1-keep, outC...))

            if include_hotspot_cons
                kmH = consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep, power=2.5, seed=seed_masks)
                outH = decomp_viability_at_mask(pool, grid; τ=τ, keepmask=kmH, A_min=A_min, p_min=p_min)
                push!(rows, (; regime=String(rname), scenario="hotspot_cons", loss=1-keep, outH...))
            end
            if include_hotspot_prey
                kmH = prey_hotspot_mask(grid, pool; τ=τ, keep_frac=keep, power=2.5, seed=seed_masks)
                outH = decomp_viability_at_mask(pool, grid; τ=τ, keepmask=kmH, A_min=A_min, p_min=p_min)
                push!(rows, (; regime=String(rname), scenario="hotspot_prey", loss=1-keep, outH...))
            end
        end
    end
    return rows
end

# ==== RUN ====
grid = make_grid(60,60; seed=11)
rows = sweep_viability(
    ; grid=grid, τ=0.58, keep_fracs=collect(0.20:0.05:0.80),
    include_hotspot_cons=true, include_hotspot_prey=true,
    A_min=0.15, p_min=0.15
)

# to DataFrame (if you use DataFrames)
# using DataFrames
# df = DataFrame(rows)
# --- compute EXCESS (scenario − random) by regime & loss ---
function excess_by(df, regime::String, loss::Float64, scenario::String)
    R  = first(filter(r -> r.regime==regime && r.scenario=="random"    && isapprox(r.loss, loss; atol=1e-8), rows))
    Sc = first(filter(r -> r.regime==regime && r.scenario==scenario   && isapprox(r.loss, loss; atol=1e-8), rows))
    (; dclim = Sc.clim - R.clim, dinter = Sc.inter - R.inter, dsyn = Sc.syn - R.syn)
end

begin
    scen = "clustered" # Choose between "clustered", "hotspot_cons", "hotspot_prey"
    fig = Figure(; size=(980, 620))
    Label(
        fig[0, 1:2], "Excess ΔU across regimes — Scenario: $(scen)",
        fontsize = 22, tellwidth = false
    )
    regnames = collect(keys(VI_regimes()))
    losses = unique(sort!(map(r->r.loss, rows)))

    col = (clim=:dodgerblue, inter=:orange, syn=:seagreen)

    for (k, rname) in enumerate(regnames)
        ax = Axis(
            fig[div(k-1,2)+1, mod(k-1,2)+1],
            title=String(rname), xlabel=(k>=3 ? "Area lost (fraction)" : ""),
            ylabel=(mod(k-1,2)==0 ? "Excess ΔU (clustered − random)" : "")
        )
        dC = [excess_by(rows, String(rname), ℓ, scen) for ℓ in losses]
        lines!(ax, losses, [d.dclim for d in dC], color=col.clim, label="Climate")
        lines!(ax, losses, [d.dinter for d in dC], color=col.inter, label="Interaction")
        lines!(ax, losses, [d.dsyn  for d in dC], color=col.syn,  label="Synergy")
        hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
        if k==3; axislegend(ax, position=:lb); end
    end
    display(fig)
end
