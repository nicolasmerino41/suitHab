# --- cell-level climate & prey-support scores from the *baseline* state ---
function cell_scores(Z0::BitMatrix, P0::BitMatrix, pool::SpeciesPool)
    cons = .!pool.basal; nb = count(pool.basal); C = size(Z0,2)
    # Climate score: fraction of consumers that are climate-suitable in the cell
    A_cell = vec(mean(view(Z0, cons, :); dims=1))
    # Prey-support score: for each consumer at each cell, fraction of its prey present
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
    # Aggregate prey support across consumers
    prey_cell = vec(mean(preyfrac; dims=1))
    return A_cell, prey_cell
end

########## DAMAGE / EXCESS HELPERS ##########
# meaning: Δ = B_after - B_before   (usually ≤ 0)
damage(Δ) = -Δ                                     # positive = worse
excess_damage(Δ_s, Δ_r) = -(Δ_s - Δ_r)             # >0 scenario worse than random
rel_excess(Δ_s, Δ_r; ε=1e-12) = excess_damage(Δ_s,Δ_r) / max(damage(Δ_r), ε)

interaction_share(ED_clim, ED_int, ED_syn; ε=1e-12) =
    abs(ED_int) / (abs(ED_clim) + abs(ED_int) + abs(ED_syn) + ε)

########## COLLECT & COMBINE ##########
"Turn a vector of NamedTuples with fields :clim,:inter,:syn into arrays."
function collect_decomp_series(vec_of_namedtuples)
    clim  = Float64[ v.clim  for v in vec_of_namedtuples ]
    inter = Float64[ v.inter for v in vec_of_namedtuples ]
    syn   = Float64[ v.syn   for v in vec_of_namedtuples ]
    (; clim, inter, syn)
end

"Given series for random (dR) and a scenario (dS), compute excess curves."
function excess_from_decomps(dR, dS)
    ED_clim = excess_damage.(dS.clim,  dR.clim)
    ED_int  = excess_damage.(dS.inter, dR.inter)
    ED_syn  = excess_damage.(dS.syn,   dR.syn)
    ED_tot  = ED_clim .+ ED_int .+ ED_syn

    ΔR_tot = dR.clim .+ dR.inter .+ dR.syn
    ΔS_tot = dS.clim .+ dS.inter .+ dS.syn
    RED_tot = rel_excess.(ΔS_tot, ΔR_tot)
    share   = interaction_share.(ED_clim, ED_int, ED_syn)
    (; ED_clim, ED_int, ED_syn, ED_tot, RED_tot, share)
end


# --- "A-matched" hotspot: within climate deciles, remove lowest-prey cells (or keep highest) ---
function a_matched_hotspot_mask(A_cell::Vector{Float64}, prey_cell::Vector{Float64};
                                keep_frac::Float64, nbins::Int=10, seed::Int=1)
    C = length(A_cell); keep = falses(C)
    rng = MersenneTwister(seed)
    # bin edges
    q = [quantile(A_cell, i/nbins) for i in 0:nbins]
    for b in 1:nbins
        idx = findall(c -> (A_cell[c] ≥ q[b]) && (A_cell[c] ≤ q[b+1]), 1:C)
        if isempty(idx); continue; end
        kkeep = round(Int, keep_frac * length(idx))
        # keep the k cells with *largest* prey support within this A bin
        ord = sortperm(prey_cell[idx]; rev=true)
        keep[idx[ord[1:kkeep]]] .= true
    end
    return keep
end

function rewire_metaweb(pool::SpeciesPool; seed::Int=42)
    rng = MersenneTwister(seed)
    nb = count(pool.basal); S = pool.S
    newE = [Int[] for _ in 1:S]
    for s in 1:nb
        newE[s] = Int[]   # basal
    end
    for s in nb+1:S
        k = length(pool.E[s])
        newE[s] = rand(rng, 1:nb, k)  # same diet size, random basal prey
    end
    SpeciesPool(pool.S, copy(pool.masses), copy(pool.basal),
                copy(pool.mu), copy(pool.b), newE)
end

"""
run_excess_with_controls(pool, grid; τ, keep_grid, losses, hotspot_mode=:A_matched)
Returns a NamedTuple with time series for each scenario and the rewired placebo.
"""
function run_excess_with_controls(pool::SpeciesPool, grid::Grid;
        τ::Float64, losses::Vector{Float64}, hotspot_mode::Symbol = :A_matched,
        seed::Int=101)

    # Baseline before any loss
    Z0 = climate_pass(pool, grid; τ=τ)
    P0 = assemble(Z0, pool)
    A_cell, prey_cell = cell_scores(Z0, P0, pool)

    # helper to build hotspot keepmask
    function hotspot_keep(keep_frac)
        if hotspot_mode == :prey_only
            # simple prey hotspot (no matching)
            C = length(A_cell); kkeep = round(Int, keep_frac*C)
            ord = sortperm(prey_cell; rev=true)
            keep = falses(C); keep[ord[1:kkeep]] .= true; return keep
        else
            return a_matched_hotspot_mask(A_cell, prey_cell; keep_frac=keep_frac, nbins=10, seed=seed)
        end
    end

    # rewired placebo
    pool_rw = rewire_metaweb(pool; seed=seed)
    Z0_rw   = climate_pass(pool_rw, grid; τ=τ)
    P0_rw   = assemble(Z0_rw, pool_rw)
    A_cell_rw, prey_cell_rw = cell_scores(Z0_rw, P0_rw, pool_rw)

    data = Dict{Symbol,Vector{NamedTuple}}(
        :random   => NamedTuple[],
        :cluster  => NamedTuple[],
        :hotspot  => NamedTuple[],
        :random_rw=> NamedTuple[],
        :cluster_rw=>NamedTuple[],
        :hotspot_rw=>NamedTuple[]
    )

    for loss in losses
        keep = 1 - loss
        kmR = random_mask(grid.C, keep; seed=seed)
        kmC = clustered_mask(grid, keep; nseeds=1, seed=seed)
        kmH = hotspot_keep(keep)

        # ORIGINAL pool
        push!(data[:random],  decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmR))
        push!(data[:cluster], decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmC))
        push!(data[:hotspot], decomp_at_mask_fraction(pool, grid; τ=τ, keepmask=kmH))

        # REWIRED placebo (reuse the *same masks* so climate piece is identical)
        push!(data[:random_rw],  decomp_at_mask_fraction(pool_rw, grid; τ=τ, keepmask=kmR))
        push!(data[:cluster_rw], decomp_at_mask_fraction(pool_rw, grid; τ=τ, keepmask=kmC))
        # hotspot defined from rewired baseline too (fair)
        kmH_rw = a_matched_hotspot_mask(A_cell_rw, prey_cell_rw; keep_frac=keep, nbins=10, seed=seed)
        push!(data[:hotspot_rw], decomp_at_mask_fraction(pool_rw, grid; τ=τ, keepmask=kmH_rw))
    end

    return (; data, losses)
end

function plot_excess_and_share(results; scen=:hotspot, label="A-matched hotspot")
    losses = results.losses
    R = results.data[:random];  S = results.data[scen]
    excess   = [S[i].dB - R[i].dB for i in eachindex(losses)]
    ex_clim  = [S[i].dA_only   - R[i].dA_only   for i in eachindex(losses)]
    ex_int   = [S[i].dInt_only - R[i].dInt_only for i in eachindex(losses)]
    ex_syn   = [S[i].synergy   - R[i].synergy   for i in eachindex(losses)]
    share    = [abs(excess[i]) < 1e-12 ? NaN : ex_int[i] / excess[i] for i in eachindex(losses)]

    fig = Figure(; size=(980,360))
    ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Excess ΔBSH (scenario − random)",
               title=label)
    lines!(ax1, losses, ex_clim, label="Climate-only")
    lines!(ax1, losses, ex_int,  label="Interaction-only")
    lines!(ax1, losses, ex_syn,  label="Synergy")
    hlines!(ax1, [0], color=(:gray,0.5), linestyle=:dash); axislegend(ax1, position=:lt)

    ax2 = Axis(fig[1,2], xlabel="Area lost (fraction)", ylabel="Interaction share of excess")
    lines!(ax2, losses, share)
    hlines!(ax2, [0,1], color=(:gray,0.4), linestyle=:dash)
    display(fig)
end

function plot_placebo(results; scen=:hotspot, label="A-matched hotspot")
    losses = results.losses
    R  = results.data[:random];    S  = results.data[scen]
    Rw = results.data[:random_rw]; Sw = results.data[Symbol(scen, "_rw")]

    ex   = [S[i].dB  - R[i].dB  for i in eachindex(losses)]
    ex_w = [Sw[i].dB - Rw[i].dB for i in eachindex(losses)]

    fig = Figure(; size=(620,320))
    ax  = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Excess ΔBSH vs random",
               title="Placebo (degree-preserving rewiring)")
    lines!(ax, losses, ex,  label="with interactions")
    lines!(ax, losses, ex_w, label="rewired placebo", linestyle=:dash)
    hlines!(ax, [0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax, position=:lt)
    display(fig)
end
