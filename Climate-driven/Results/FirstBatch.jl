# --- Setup ---
grid = make_grid(60, 60; seed=11)
τ = 0.55
# Low redundancy + high prey synchrony (what gave us the clearest effects)
pool = build_pool(200;
    basal_frac=0.35, seed=1,
    sigma=0.22, density=0.12, pmax=0.70,
    R0_mean=10.0, R0_sd=0.25,
    niche_mode=:bimodal, mu_basal_centers=(0.25, 0.75), mu_basal_sd=0.04,
    b0_basal=0.08, bspread_basal=0.02,
    b0_cons=0.12, bspread_cons=0.04
)

losses = 0.0:0.05:0.8  # fraction of area lost

#### FIGURE A ####
# Helper: one curve for a given mask type
function decomp_curve(kind; hotspot_power=2.0, nseeds=6)
    xs = Float64[]; c_only = Float64[]; i_only = Float64[]; syn = Float64[]
    for f in losses
        keep = 1.0 - f
        score = prey_vulnerability_score(pool, grid; τ=τ)
        # score = Float64.(score)
        # score .= (score .+ 1e-12) .^ power
        keepmask = kind === :random    ? random_mask(grid.C, keep; seed=123) :
                   kind === :clustered ? clustered_mask(grid, keep; nseeds=nseeds, seed=123) :
                #    hotspot_clustered_mask_bestfirst(grid, keep; score=score, nseeds=1, seed=303)
                   consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep, power=hotspot_power)
        d = decomp_at_mask(pool, grid; τ=τ, keepmask=keepmask)  # returns .agg with fields
        d = decomp_at_mask(pool, grid; τ=τ, keepmask=keepmask)
        push!(xs, f)
        push!(c_only, d.agg.dAcl); push!(i_only, d.agg.dInt); push!(syn, d.agg.dSyn)
    end
    (; loss=xs, climate=c_only, interact=i_only, synergy=syn)
end

curR = decomp_curve(:random)
curC = decomp_curve(:clustered; nseeds=1)           # 1 seed makes clusters more coherent
curH = decomp_curve(:hotspot; hotspot_power=2.5)    # stronger targeting, bigger effect

# Plot three panels (Random / Clustered / Hotspot)
begin
    fig = Figure(; size=(1180,380))
    curves = [curR, curC, curH]; titles = ["Random", "Clustered", "Hotspot"]
    for j in 1:3
        ax = Axis(fig[1,j], title=titles[j], xlabel="Area lost (fraction)",
                  ylabel = j==1 ? "ΔBSH decomposition" : "")
        lines!(ax, curves[j].loss, curves[j].climate,   label=j==3 ? "Climate-only" : nothing)
        lines!(ax, curves[j].loss, curves[j].interact,  label=j==3 ? "Interaction-only" : nothing)
        lines!(ax, curves[j].loss, curves[j].synergy,   label=j==3 ? "Synergy" : nothing)
        hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
        if j==3; axislegend(ax, position=:lt); end
    end
    display(fig)
end

#### FIGURE B ####
f = 0.5; keep = 1.0 - f
km_hot = consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep, power=2.5)
# km_hot = hotspot_clustered_mask_bestfirst(grid, keep; score=score, nseeds=1, seed=303)
km_rnd = random_mask(grid.C, keep; seed=99)

function ΔA_Δp(pool, grid; τ, keepmask)
    Z0 = climate_pass(pool, grid; τ=τ)
    P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask)
    P1 = assemble(Z1, pool)
    cons = .!pool.basal
    A0 = climate_area_per_species(Z0)[cons]
    A1 = climate_area_per_species(Z1)[cons]
    p0 = bsh1_cond_per_species(P0, Z0, pool)[cons]
    p1 = bsh1_cond_per_species(P1, Z1, pool)[cons]
    (ΔA=A1.-A0, Δp=p1.-p0)
end

ΔH = ΔA_Δp(pool, grid; τ=τ, keepmask=km_hot)
ΔR = ΔA_Δp(pool, grid; τ=τ, keepmask=km_rnd)

# Quick scatter (hotspot)
begin
    fig = Figure(; size=(560,380))
    ax = Axis(fig[1,1], xlabel="ΔA", ylabel="Δp", title="Hotspot, loss=0.5")
    scatter!(ax, ΔH.ΔA, ΔH.Δp, markersize=5)
    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    vlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    display(fig)
end

#### FIGURE C ####
function ΔB_mean(pool, grid; τ, keepmask)
    Z0 = climate_pass(pool, grid; τ=τ); P0 = assemble(Z0, pool)
    Z1 = apply_mask(Z0, keepmask);      P1 = assemble(Z1, pool)
    cons = .!pool.basal
    B0 = bsh1_per_species(P0, Z0, pool)[cons]
    B1 = bsh1_per_species(P1, Z1, pool)[cons]
    mean(B1 .- B0)
end

xs = Float64[]; excessH = Float64[]; excessC = Float64[]
for f in losses
    keep = 1.0 - f
    kmR = random_mask(grid.C, keep; seed=10)
    kmC = clustered_mask(grid, keep; nseeds=1, seed=10)
    kmH = consumer_hotspot_mask(grid, pool; τ=τ, keep_frac=keep, power=2.5)
    # kmH = hotspot_clustered_mask_bestfirst(grid, keep; score=score, nseeds=1, seed=303)
    dR = ΔB_mean(pool, grid; τ=τ, keepmask=kmR)
    dC = ΔB_mean(pool, grid; τ=τ, keepmask=kmC)
    dH = ΔB_mean(pool, grid; τ=τ, keepmask=kmH)
    push!(xs, f); push!(excessC, dC - dR); push!(excessH, dH - dR)
end

begin
    fig = Figure(; size=(700,360))
    ax  = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="Excess ΔB vs random")
    lines!(ax, xs, excessC, label="Clustered − Random")
    lines!(ax, xs, excessH, label="Hotspot − Random")
    hlines!(ax, [0.0], color=(:gray,0.5), linestyle=:dash)
    axislegend(ax, position=:lt)
    display(fig)
end
