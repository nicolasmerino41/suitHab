# ---------------------- run_kofn_transition.jl -------------------------------
# Narrative 3: With k-of-n prey rules, retained area shows a transition. Random
# shifts the critical loss left; front shifts it right; redundancy pushes right.
# Figures:
#   F1_y_vs_f_k123.png, F2_fc_vs_k.png
# ------------------------------------------------------------------------------

# --- modules
include("../SetUp.jl")
include("src/grids.jl");      using .Grids
include("src/grids_init.jl"); using .GridsInit
include("src/metawebs.jl");   using .Metawebs
include("src/hl.jl");         using .HL
include("src/bsh.jl");        using .BSH

# --- settings
rng        = MersenneTwister(123)
nx, ny     = 120, 120
S          = 200
basal_frac = 0.25
loss_fracs = 0.20:0.05:0.80
geoms      = (:random,:clustered,:front)
Klist      = [1,2,3]                 # kreq values
outdir     = "Paper/narratives/kofn"; isdir(outdir) || mkpath(outdir)

# --- build grid & pools
include("src/grids_init.jl"); using .GridsInit
G     = GridsInit.build_default_grids(nx, ny; seeds=(11,12,13,14))
grid  = G.grad
pool  = Metawebs.build_metaweb_archetype(rng; S, basal_frac, archetype=:mid)
pars  = BSH.BAMParams(; τA=0.5, τB=0.5, τocc=0.2, γ=3.0, movement=:component, T=8)

A = BSH.abiotic_matrix(pool, grid; seed=1)
cons_inds = findall(!, pool.basal)

# --- compute mean retained area vs f for AM and BAM(k)
function mean_area_vs_f(; pool, grid, pars, A, geom::Symbol, kreq::Int)
    xs = Float64[]; yAM = Float64[]; yBM = Float64[]
    for f in loss_fracs
        keepfrac = 1 - f
        keep = geom===:random    ? HL.random_mask(rng, grid.C, keepfrac) :
               geom===:clustered ? HL.clustered_mask(rng, grid.nx, grid.ny, keepfrac; nseeds=8) :
                                   HL.front_mask(rng, grid.xy, keepfrac; axis=:x, noise=0.05)
        _, bAM = BSH.assemble_AM(pool, grid, A, keep; pars=pars)
        bBM    = BSH.assemble_BAM(pool, grid, A, keep; pars=pars, agg=:kofn, kreq=kreq).bsh
        push!(xs, f); push!(yAM, mean(bAM[cons_inds])); push!(yBM, mean(bBM[cons_inds]))
    end
    (; f=xs, AM=yAM, BAM=yBM)
end

# --- curvature & critical f
function critical_f(f::Vector{Float64}, y::Vector{Float64})
    # discrete 2nd derivative
    n=length(f); n<3 && return f[clamp(argmin(y),1,n)]
    curv = similar(y); curv .= NaN
    for i in 2:n-1
        df1 = f[i]   - f[i-1]; df2 = f[i+1] - f[i]
        dy1 = y[i]   - y[i-1]; dy2 = y[i+1] - y[i]
        s1  = dy1/df1; s2 = dy2/df2
        curv[i] = (s2 - s1) / ((df1 + df2)/2)
    end
    # pick most negative curvature (steepening drop)
    idx = argmin(curv)         # (ties handled by first)
    return f[idx], curv
end

# --- Fig 1: y(f) curves for k=1,2,3 for each geometry
begin
    f1 = Figure(; size=(1100,350))
    for (col, g) in enumerate(geoms)
        ax = Axis(f1[1,col], title=String(g), xlabel="area lost f", ylabel="retained area")
        for kreq in Klist
            dat = mean_area_vs_f(; pool, grid, pars, A, geom=g, kreq=kreq)
            lines!(ax, dat.f, dat.BAM, label="BAM k=$kreq")
        end
        lines!(ax, collect(loss_fracs),
            mean_area_vs_f(; pool, grid, pars, A, geom=g, kreq=1).AM,
            color=:black, linestyle=:dash, label="AM")
        axislegend(ax, position=:rt)
    end
    display(f1)
end
save(joinpath(outdir,"F1_y_vs_f_k123.png"), f1)

# --- Fig 2: f_c vs k per geometry
fcs = Dict(g => Float64[] for g in geoms)
for g in geoms, kreq in Klist
    dat = mean_area_vs_f(; pool, grid, pars, A, geom=g, kreq=kreq)
    fc,_ = critical_f(dat.f, dat.BAM)
    push!(fcs[g], fc)
end

begin
    f2 = Figure(; size=(800,360))
    ax2 = Axis(
        f2[1,1], title="Critical loss f_c vs k (BAM)",
        xlabel="k (required prey)", ylabel="f_c"
    )
    cols = Dict(:random=>:steelblue, :clustered=>:orange, :front=>:green)
    for g in geoms
    lines!(ax2, Klist, fcs[g]; color=cols[g], label=String(g))
    end
    axislegend(ax2, position=:rb)
    display(f2)
end
save(joinpath(outdir,"F2_fc_vs_k.png"), f2)

@info "k-of-n narrative figs saved in $outdir"
