# Point 2 — Regional viability curves (share of consumers that “fail”)
function viability_curve(pool::SpeciesPool, grid::Grid;
                         τ::Float64=0.64,
                         losses=collect(0.0:0.05:0.8),
                         ρ::Float64=0.10,
                         strategy::Symbol=:random,
                         nseeds_cluster::Int=1,
                         score_type::Symbol=:support)
    cons = .!pool.basal
    Z0 = climate_pass(pool, grid; τ=τ)
    P0 = assemble(Z0, pool)
    A0 = sum(Z0; dims=2)[:]                      # baseline climate-suitable cells

    # hotspot score (if needed)
    score = score_type === :support ?
        prey_support_score(pool, grid; τ=τ) :
        prey_vulnerability_score(pool, grid; τ=τ)

    out = Float64[]
    for f in losses
        keep = 1.0 - f
        km = strategy === :random    ? random_mask(grid.C, keep; seed=101) :
             strategy === :clustered ? clustered_mask(grid, keep; nseeds=nseeds_cluster, seed=202) :
             hotspot_clustered_mask_bestfirst(grid, keep; score=score, nseeds=nseeds_cluster, seed=303)

        Z1 = apply_mask(Z0, km)
        P1 = assemble(Z1, pool)
        C1 = size(Z1, 2)                               # << correct number of columns

        B1 = zeros(Int, pool.S)                        # effective area in cells
        @inbounds for s in 1:pool.S
            if cons[s]
                cnt = 0
                for c in 1:C1
                    if Z1[s,c]
                        has_pre = false
                        for q in pool.E[s]; if P1[q,c]; has_pre = true; break; end; end
                        cnt += has_pre ? 1 : 0
                    end
                end
                B1[s] = cnt
            else
                B1[s] = sum(@view Z1[s, :])           # basal: climate-only
            end
        end

        thr = @. max(1, ceil(Int, ρ * A0))
        use = cons .& (A0 .> 0)
        frac_lost = sum(B1[use] .< thr[use]) / sum(use)
        push!(out, frac_lost)
    end
    return (losses=losses, frac=out)
end

function plot_viability_three(pool::SpeciesPool, grid::Grid; τ=0.64, losses=collect(0.0:0.05:0.8), ρ=0.10)
    r = viability_curve(pool, grid; τ=τ, losses=losses, ρ=ρ, strategy=:random)
    c = viability_curve(pool, grid; τ=τ, losses=losses, ρ=ρ, strategy=:clustered)
    h = viability_curve(pool, grid; τ=τ, losses=losses, ρ=ρ, strategy=:hotspot, score_type=:vulnerable)
    begin
        fig = Figure(resolution=(760,360))
        ax  = Axis(fig[1,1], title="Consumers losing viability (ρ=$(ρ))",
                   xlabel="Area lost (fraction)", ylabel="Fraction")
        lines!(ax, r.losses, r.frac, label="Random")
        lines!(ax, c.losses, c.frac, label="Clustered")
        lines!(ax, h.losses, h.frac, label="Hotspot (vulnerable)")
        axislegend(ax, position=:lt)
        display(fig)
    end
end

plot_viability_three(pool, grid; τ=τ, losses=losses, ρ=0.10)
