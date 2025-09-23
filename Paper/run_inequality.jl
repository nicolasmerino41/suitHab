# ----------------------------- run_inequality.jl ------------------------------
# Narrative 1: Biotic interactions increase inequality of species-level losses
# Figures:
#   F1_Δgini_bars.png, F2_Δtail_bars.png, F3_regime_map_DeltaGini.png, F4_cdf.png
# ------------------------------------------------------------------------------

# --- modules
include("../SetUp.jl")
include("src/grids.jl");      using .Grids
include("src/grids_init.jl"); using .GridsInit
include("src/metawebs.jl");   using .Metawebs
include("src/hl.jl");         using .HL
include("src/bsh.jl");        using .BSH
include("src/metrics.jl");    using .Metrics


# --- settings
rng        = MersenneTwister(123)
nx, ny     = 120, 120
S          = 200
basal_frac = 0.25
fstar      = 0.60
geoms      = (:random,:clustered,:front)
outdir     = "Paper/narratives/inequality"; isdir(outdir) || mkpath(outdir)

# --- data
G     = GridsInit.build_default_grids(nx, ny; seeds=(11,12,13,14))
grid  = G.grad                       # gradient grid for the headline bars
p_mid = Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:mid)
pars  = BSH.BAMParams(; τA=0.5, τB=0.5, τocc=0.2, γ=3.0, movement=:component, T=8)

# --- helpers
gini(x) = begin
    y=sort(clamp.(x,0,1)); n=length(y)
    n==0 && return 0.0
    2*sum((1:n).*y)/(n*sum(y)+eps()) - (n+1)/n
end

function per_species_losses(; rng, pool, grid, pars, fstar, geometry, seed_A=1)
    rA, rB = BSH.per_species_relative_loss(rng, pool, grid, pars; fstar, geometry, seed_A)
    LA  = clamp.(-rA, 0.0, 1.0)     # positive loss in [0,1]
    LB  = clamp.(-rB, 0.0, 1.0)
    return LA, LB
end

# --- Fig 1 & 2: ΔGini and ΔTail bars by geometry
tail_cut = 0.80
Δgini = Float64[]; Δtail = Float64[]
for g in geoms
    LA, LB = per_species_losses(; rng, pool=p_mid, grid, pars, fstar, geometry=g, seed_A=1)
    push!(Δgini, gini(LB) - gini(LA))
    push!(Δtail, mean(LB .> tail_cut) - mean(LA .> tail_cut))
end

labels = collect(String.(geoms))  # ["random", "clustered", "front"]

begin
    f1 = Figure(; size=(800,350))
    ax1 = Axis(f1[1,1];
        title="ΔGini (BAM − AM) at f*=$(fstar)",
        ylabel="ΔGini",
        xticks=(1:3, ["random","clustered","front"])
    )
    barplot!(ax1, 1:3, Δgini, color=:steelblue)
    hlines!(ax1, [0.0], color=:gray, linestyle=:dash)
    display(f1)
end

save(joinpath(outdir,"F1_Δgini_bars.png"), f1)

begin
    f2 = Figure(; size=(800,350))
    ax2 = Axis(
        f2[1,1];
        title="ΔTail (BAM − AM), tail>$(tail_cut) at f*=$(fstar)",
        ylabel="ΔTail",
        xticks=(1:3, ["random","clustered","front"])
    )
    barplot!(ax2, 1:3, Δtail, color=:steelblue)
    hlines!(ax2, [0.0], color=:gray, linestyle=:dash)
    display(f2)
end
save(joinpath(outdir,"F2_Δtail_bars.png"), f2)

# --- Fig 3: Regime map ΔGini over (D,R) across archetypes × grids
grids  = (G.grad, G.patch, G.mosaic, G.ridge)
gnames = ("gradient","patchy","mosaic","ridge")
pools  = (Metawebs.build_metaweb_archetype(MersenneTwister(7);  S, basal_frac, archetype=:low),
          Metawebs.build_metaweb_archetype(MersenneTwister(8);  S, basal_frac, archetype=:mid),
          Metawebs.build_metaweb_archetype(MersenneTwister(10); S, basal_frac, archetype=:high))
rnames = ("low","mid","high")

Dvals = [Grids.climate_tail_index(g) for g in grids]
Rvals = Float64[]
for p in pools
    diets = [length(p.prey[s]) for s in 1:p.S if !p.basal[s]]
    push!(Rvals, isempty(diets) ? 0.0 : quantile(diets, 0.95))
end

Z = zeros(length(Rvals), length(Dvals))
for (j,g) in enumerate(grids), (i,p) in enumerate(pools)
    LA, LB = per_species_losses(; rng, pool=p, grid=g, pars, fstar, geometry=:random, seed_A=1)
    Z[i,j] = gini(LB) - gini(LA)
end

begin
    f3 = Figure(; size=(700,420))
    ax = Axis(
        f3[1,1], title="ΔGini (random geometry) over (D,R)",
        xlabel="Climate tail index D", ylabel="Redundancy R"
    )

    hm = heatmap!(ax, 1:length(Dvals), 1:length(Rvals), Z; colormap=:viridis)

    # convert tuples → vectors for labels
    ax.xticks = (1:length(Dvals), collect(string.(gnames)))
    ax.yticks = (1:length(Rvals), collect(string.(rnames)))

    Colorbar(f3[1,2], hm, label="ΔGini")
    display(f3)

end
save(joinpath(outdir,"F3_regime_map_DeltaGini.png"), f3)

# --- Fig 4: A single CDF panel (front geometry) for storytelling
geom = :front
LA, LB = per_species_losses(; rng, pool=p_mid, grid, pars, fstar, geometry=geom, seed_A=1)
function ecdf_xy(x)
    s = sort(x); n=length(s); (s, (1:n)./n)
end
xA,yA = ecdf_xy(LA); xB,yB = ecdf_xy(LB)
begin
    f4 = Figure(; size=(800,450))
    a4 = Axis(
        f4[1,1], title="Per-species loss CDFs — $(String(geom)) at f*=$(fstar)",
        xlabel="loss (−Δ/BSH₀)", ylabel="fraction of consumers"
    )
    lines!(a4, xA, yA, color=:gray, linestyle=:dash, label="AM")
    lines!(a4, xB, yB, color=:darkgreen, label="BAM")
    axislegend(a4, position=:rb)
    display(f4)
end
save(joinpath(outdir,"F4_cdf.png"), f4)

@info "Inequality narrative figs saved in $outdir"
