# --- species-level trophic level from the metaweb (acyclic via body size) ---
# TL=1 for basal; otherwise TL = 1 + max(TL of prey)
function trophic_level_metaweb(pool::SpeciesPool)
    S = pool.S
    TL = zeros(Int, S)
    ord = sortperm(pool.masses)  # increasing mass (prey before predators)
    for s in ord
        if pool.basal[s] || isempty(pool.E[s])
            TL[s] = 1
        else
            TL[s] = 1 + maximum(TL[pool.E[s]])
        end
    end
    return TL
end

# --- group TL into bins (edit if you want different bins) ---
# returns a Symbol per consumer species: :TL2, :TL3, :TL4p
function tl_group_labels(TL::Vector{Int}, basal::BitVector)
    S = length(TL)
    labs = Vector{Symbol}(undef, S)
    @inbounds for s in 1:S
        if basal[s]
            labs[s] = :basal
        else
            t = TL[s]
            labs[s] = t <= 2 ? :TL2 : (t == 3 ? :TL3 : :TL4p)
        end
    end
    return labs
end

# --- keepmask dispatcher (uses your existing masks incl. hotspot) ---
function keepmask_for(kind::Symbol, grid::Grid, pool::SpeciesPool;
                      τ::Float64, keep::Float64, nseeds_cluster::Int, seed::Int,
                      hotspot_power::Float64=2.5, hotspot_what::Symbol=:cons)
    if kind === :random
        return random_mask(grid.C, keep; seed=seed)
    elseif kind === :clustered
        return clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=seed)
    elseif kind === :hotspot
        return consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep,
                                     power=hotspot_power, what=hotspot_what, seed=seed)
    else
        error("unknown kind = $kind")
    end
end

# --- FRACTION-basis per-TL decomposition at one mask ---
# B_s = A_s * p_s, and ΔB_s = p̄_s ΔA_s + Ā_s Δp_s + ΔA_s Δp_s (midpoint)
function decomp_by_TL_fraction(pool::SpeciesPool, grid::Grid;
                               τ::Float64, keepmask::BitVector)
    Z0 = climate_pass(pool, grid; τ=τ)
    P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask)
    P1 = assemble(Z1, pool)

    cons = .!pool.basal
    TL_all = trophic_level_metaweb(pool)
    labs   = tl_group_labels(TL_all, pool.basal)

    A0 = climate_area_per_species(Z0)
    p0 = bsh1_cond_per_species(P0, Z0, pool)
    B0 = A0 .* p0

    A1 = climate_area_per_species(Z1)
    p1 = bsh1_cond_per_species(P1, Z1, pool)
    B1 = A1 .* p1

    dA = A1 .- A0
    dp = p1 .- p0
    dB = B1 .- B0
    Abar = 0.5 .* (A0 .+ A1)
    pbar = 0.5 .* (p0 .+ p1)

    parts = Dict(:TL2 => (Float64[], Float64[], Float64[]),
                 :TL3 => (Float64[], Float64[], Float64[]),
                 :TL4p=> (Float64[], Float64[], Float64[]))

    @inbounds for s in 1:pool.S
        (cons[s] && isfinite(p0[s]) && isfinite(p1[s])) || continue
        g = labs[s]; g === :basal && continue
        clim = pbar[s] * dA[s]
        inter= Abar[s] * dp[s]
        syn  = dA[s] * dp[s]
        push!(parts[g][1], clim)
        push!(parts[g][2], inter)
        push!(parts[g][3], syn)
    end

    # return means per group; missing groups get zeros
    out = Dict{Symbol,NamedTuple}()
    for g in (:TL2, :TL3, :TL4p)
        if isempty(parts[g][1])
            out[g] = (; climate=0.0, interaction=0.0, synergy=0.0)
        else
            out[g] = (; climate=mean(parts[g][1]),
                       interaction=mean(parts[g][2]),
                       synergy=mean(parts[g][3]))
        end
    end
    return out
end

# --- curves across loss for one scenario ---
function decomp_TL_curves(kind::Symbol; grid::Grid, pool::SpeciesPool,
                          τ::Float64=0.55, losses=0.0:0.05:0.8,
                          nseeds_cluster::Int=1, seed::Int=123,
                          hotspot_power::Float64=2.5, hotspot_what::Symbol=:cons)
    curves = Dict(
        :TL2 => (loss=Float64[], climate=Float64[], interaction=Float64[], synergy=Float64[]),
        :TL3 => (loss=Float64[], climate=Float64[], interaction=Float64[], synergy=Float64[]),
        :TL4p=> (loss=Float64[], climate=Float64[], interaction=Float64[], synergy=Float64[])
    )
    for f in losses
        keep = 1.0 - f
        km = keepmask_for(kind, grid, pool; τ=τ, keep=keep,
                          nseeds_cluster=nseeds_cluster, seed=seed,
                          hotspot_power=hotspot_power, hotspot_what=hotspot_what)
        d = decomp_by_TL_fraction(pool, grid; τ=τ, keepmask=km)
        for g in (:TL2, :TL3, :TL4p)
            push!(curves[g].loss, f)
            push!(curves[g].climate,     d[g].climate)
            push!(curves[g].interaction, d[g].interaction)
            push!(curves[g].synergy,     d[g].synergy)
        end
    end
    return curves
end

# --- build a pool that shows interaction differences clearly (adjust if you like) ---
grid = make_grid(60, 60; seed=11)
τ = 0.55
pool = build_pool(200;
    basal_frac=0.35, seed=1,
    sigma=0.22, density=0.12, pmax=0.70,
    R0_mean=10.0, R0_sd=0.25,
    niche_mode=:bimodal, mu_basal_centers=(0.25, 0.75), mu_basal_sd=0.04,
    b0_basal=0.08, bspread_basal=0.02,
    b0_cons=0.12, bspread_cons=0.04
)
losses = 0.0:0.05:0.8

curR = decomp_TL_curves(:random;   grid=grid, pool=pool, τ=τ, losses=losses)
curC = decomp_TL_curves(:clustered;grid=grid, pool=pool, τ=τ, losses=losses, nseeds_cluster=1)
curH = decomp_TL_curves(:hotspot;  grid=grid, pool=pool, τ=τ, losses=losses,
                         hotspot_power=2.5, hotspot_what=:cons)

# --- plot grid: rows = components, cols = scenarios; lines = TL groups ---
begin
    comps  = (:climate, :interaction, :synergy)
    titles = Dict(:random=>"Random", :clustered=>"Clustered", :hotspot=>"Hotspot")
    colors = Dict(:TL2=>:dodgerblue, :TL3=>:orange, :TL4p=>:forestgreen)

    fig = Figure(; size=(1180, 900))
    sets = [curR, curC, curH]

    for (j, curves) in enumerate(sets)             # columns (scenarios)
        for (i, comp) in enumerate(comps)          # rows (components)
            ax = Axis(fig[i, j],
                      title = (i==1 ? titles[[ :random, :clustered, :hotspot ][j]] : ""),
                      xlabel = (i==length(comps) ? "Area lost (fraction)" : ""),
                      ylabel = (j==1 ? string(comp, " ΔBSH") : ""))
            for g in (:TL2, :TL3, :TL4p)
                lines!(ax, curves[g].loss, getfield(curves[g], comp);
                       color=colors[g], label=(j==3 && i==1 ? string(g) : nothing), linewidth=2)
            end
            hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
            if j==3 && i==1
                axislegend(ax, position=:lt)
            end
        end
    end
    display(fig)
end

begin
    # data you already have:
    # curR, curC, curH :: Dict(:TL2,:TL3,:TL4p => (loss, climate, interaction, synergy))
    # losses :: 0.0:0.05:0.8 (or whatever you used)

    scenarios = [(curR, "Random"), (curC, "Clustered"), (curH, "Hotspot")]
    comps     = (:climate, :interaction, :synergy)
    compnames = Dict(:climate=>"Climate ΔBSH", :interaction=>"Interaction ΔBSH", :synergy=>"Synergy ΔBSH")
    tl_order  = (:TL2, :TL3, :TL4p)
    tl_names  = Dict(:TL2=>"TL2", :TL3=>"TL3", :TL4p=>"TL4+")
    tl_cols   = Dict(:TL2=>:dodgerblue, :TL3=>:orange, :TL4p=>:forestgreen)

    fig = Figure(; size=(1180, 900))

    # a small legend at the top
    legend_elems = [PolyElement(color=tl_cols[g]) for g in tl_order]
    legend_labels = [tl_names[g] for g in tl_order]
    # Legend(fig, legend_elems, legend_labels; orientation=:horizontal)
    fig[0, 1:3] = GridLayout()

    for (j, (curves, scename)) in enumerate(scenarios)             # columns
        for (i, comp) in enumerate(comps)                          # rows
            ax = Axis(fig[i, j],
                      title = (i == 1 ? scename : ""),
                      xlabel = (i == length(comps) ? "Area lost (fraction)" : ""),
                      ylabel = (j == 1 ? compnames[comp] : ""))

            n = length(curves[:TL2].loss)
            # X positions per loss (1..n), repeated for each TL (inner=3)
            x = repeat(1:n, inner=3)
            # group ids for dodging (1=TL2, 2=TL3, 3=TL4+), repeated per loss
            g = repeat(1:3, outer=n)

            # pull values per TL for the selected component
            v2 = getfield(curves[:TL2], comp)
            v3 = getfield(curves[:TL3], comp)
            v4 = getfield(curves[:TL4p], comp)

            # interleave by loss: [v2[1], v3[1], v4[1], v2[2], v3[2], v4[2], ...]
            y = Float64[]
            for k in 1:n
                append!(y, (v2[k], v3[k], v4[k]))
            end

            # color by group
            cols = map(gi -> (gi == 1 ? tl_cols[:TL2] : gi == 2 ? tl_cols[:TL3] : tl_cols[:TL4p]), g)

            # draw bars (dodged)
            barplot!(ax, x, y; dodge = g, color = cols, width = 0.8/3)

            # zero line & tidy ticks (show every ~second loss to reduce clutter)
            hlines!(ax, [0.0], color=(:gray, 0.5), linestyle=:dash)
            tick_idx = collect(1:2:n)
            # xticks!(ax, tick_idx, string.(round.(losses[tick_idx]; digits=2)))
        end
    end

    display(fig)
end

begin

    # expect these to be in scope already:
    # curR, curC, curH :: Dict(:TL2,:TL3,:TL4p => (loss, climate, interaction, synergy))

    scenarios  = [(curR, "Random"), (curC, "Clustered"), (curH, "Hotspot")]
    comps      = (:climate, :interaction, :synergy)
    comp_names = ["Climate", "Interaction", "Synergy"]
    tl_order   = (:TL2, :TL3, :TL4p)
    tl_names   = Dict(:TL2=>"TL2", :TL3=>"TL3", :TL4p=>"TL4+")
    tl_cols    = Dict(:TL2=>:dodgerblue, :TL3=>:orange, :TL4p=>:forestgreen)

    loss_target = 0.50

    fig = Figure(; size=(1000, 360))

    for (j, (curves, ttl)) in enumerate(scenarios)
        ax = Axis(fig[1, j],
                  title = ttl,
                  xlabel = (j == 2 ? "Component" : ""),
                  ylabel = (j == 1 ? "ΔBSH at 50% loss" : ""))

        # locate index closest to the requested loss in this curves set
        L = curves[:TL2].loss
        idx = argmin(abs.(L .- loss_target))

        # x positions are the three components; each has 3 dodged bars (one per TL)
        x = repeat(1:3, inner=length(tl_order))                 # 1..3 repeated for TLs
        g = repeat(1:length(tl_order), outer=3)                 # group ids for dodge

        # interleave values as [TL2(clim), TL3(clim), TL4p(clim), TL2(int), ...]
        y = Float64[]
        for comp in comps
            for tl in tl_order
                vals = getfield(curves[tl], comp)
                push!(y, vals[idx])
            end
        end

        cols = map(gi -> gi==1 ? tl_cols[:TL2] : gi==2 ? tl_cols[:TL3] : tl_cols[:TL4p], g)

        barplot!(ax, x, y; dodge=g, color=cols, width=0.8/length(tl_order))
        hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
        ax.xticks = (1:3, comp_names)

        if j == 2
            elems  = [PolyElement(color=tl_cols[k]) for k in tl_order]
            labels = [tl_names[k] for k in tl_order]
            axislegend(ax, elems, labels; position=:ct, orientation=:vertical)
        end
    end

    display(fig)
end
