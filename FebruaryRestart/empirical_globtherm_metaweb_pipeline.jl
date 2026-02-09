#!/usr/bin/env julia
# ============================================================
# EMPIRICAL GLOBTHERM × EU META-WEB PIPELINE (CANONICAL)
#
# PURPOSE
# -------
# • Load GlobalTherm thermal limits
# • Canonicalize species names (canon_taxon)
# • Map EU predator–prey metaweb onto GlobalTherm traits
# • Compute:
#     - Edge-level predator–prey trait relationships
#     - Predator vs mean(prey) relationships
#     - Directed thermal overlap
# • Run degree-preserving null models for ALL metrics
# • Automatically generate CSV summaries + plots
#
# THIS SCRIPT REPLACES:
#   GlobTherm.jl
#   GlobTherm2.jl
#   GlobTherm3.jl
#   GlobTherm3_1.jl
# ============================================================

using CSV, DataFrames, Statistics, Random
using CairoMakie

# ============================================================
# 0) PATHS (OVERWRITE MODE)
# ============================================================

script_dir   = @__DIR__
project_root = dirname(script_dir)
data_dir     = joinpath(project_root, "data")

globtherm_csv = joinpath(project_root, "FebruaryRestart",
                         "GlobalTherm_upload_02_11_17.csv")
pairwise_csv  = joinpath(data_dir,
                         "TetraEU_pairwise_interactions.csv")

OUTDIR = joinpath(project_root, "outputs_empirical_GlobTherm")
isdir(OUTDIR) || mkpath(OUTDIR)

# ============================================================
# 1) CANONICAL TAXONOMY
# ============================================================
function canon_taxon(s::AbstractString)::String
    t = replace(String(s), r"[_(),\[\]]" => " ")
    t = replace(t, r"\s+" => " ")
    t = strip(t)
    parts = split(t)
    length(parts) < 2 && return ""

    genus   = String(filter(isletter, parts[1]))
    species = String(filter(isletter, parts[2]))

    isempty(genus) || isempty(species) && return ""
    lowercase(species) in ("sp", "spp") && return ""

    genus2   = uppercase(genus[1]) * lowercase(genus[2:end])
    species2 = lowercase(species)

    return "$genus2 $species2"
end

@inline function safe_get(v, i)
    (1 ≤ i ≤ length(v) && isassigned(v, i)) ? v[i] : missing
end

# ============================================================
# 2) LOAD DATA
# ============================================================
@info "Loading GlobalTherm"
gt = CSV.read(globtherm_csv, DataFrame; missingstring="", ntasks=1, pool=false)

@info "Loading EU metaweb"
raw_web = CSV.read(pairwise_csv, DataFrame)

web = DataFrame(
    pred = canon_taxon.(string.(raw_web.sourceTaxonName)),
    prey = canon_taxon.(string.(raw_web.targetTaxonName))
)

filter!(r -> !isempty(r.pred) && !isempty(r.prey), web)

# ============================================================
# 3) BUILD GLOBALTHERM TRAIT TABLE
# ============================================================
function pick_first_numeric(vals...)
    for v in vals
        v === missing && continue
        v isa Real && return Float64(v)
        v isa AbstractString && (x=tryparse(Float64,v); x!==nothing && return x)
    end
    return NaN
end

canon = String[]
Tmin  = Float64[]
Tmax  = Float64[]

for i in 1:nrow(gt)
    g = safe_get(gt[!, :Genus], i)
    s = safe_get(gt[!, :Species], i)
    g === missing || s === missing && continue

    cname = canon_taxon(string(g," ",s))
    isempty(cname) && continue

    tmin = pick_first_numeric(gt[i,:tmin], gt[i,:tmin_2])
    tmax = pick_first_numeric(gt[i,:Tmax], gt[i,:Tmax_2])

    isfinite(tmin) && isfinite(tmax) && tmax > tmin || continue

    push!(canon, cname)
    push!(Tmin, tmin)
    push!(Tmax, tmax)
end

gt_raw = DataFrame(
    canon = canon,
    Tmin  = Tmin,
    Tmax  = Tmax
)

gt_sum = combine(groupby(gt_raw, :canon),
    :Tmin => mean => :Tmin,
    :Tmax => mean => :Tmax
)

gt_sum.Tmid    = (gt_sum.Tmin .+ gt_sum.Tmax) ./ 2
gt_sum.breadth = gt_sum.Tmax .- gt_sum.Tmin

gt_lookup = Dict(r.canon => r for r in eachrow(gt_sum))

# ============================================================
# 4) MAP META-WEB TO TRAITS
# ============================================================
filter!(r -> haskey(gt_lookup,r.pred) && haskey(gt_lookup,r.prey), web)

@info "Matched interactions" nrow(web)

# ============================================================
# 5) EDGE-LEVEL TABLE
# ============================================================
function interval_overlap(a1,a2,b1,b2)
    max(0.0, min(a2,b2) - max(a1,b1))
end

edges = DataFrame(
    pred = web.pred,
    prey = web.prey
)

edges.pred_Tmin = [gt_lookup[p].Tmin for p in edges.pred]
edges.pred_Tmax = [gt_lookup[p].Tmax for p in edges.pred]
edges.pred_Tmid = [gt_lookup[p].Tmid for p in edges.pred]
edges.pred_breadth = [gt_lookup[p].breadth for p in edges.pred]

edges.prey_Tmin = [gt_lookup[p].Tmin for p in edges.prey]
edges.prey_Tmax = [gt_lookup[p].Tmax for p in edges.prey]
edges.prey_Tmid = [gt_lookup[p].Tmid for p in edges.prey]
edges.prey_breadth = [gt_lookup[p].breadth for p in edges.prey]

edges.overlap_pred = [
    interval_overlap(edges.pred_Tmin[i], edges.pred_Tmax[i],
                     edges.prey_Tmin[i], edges.prey_Tmax[i]) /
    max(edges.pred_breadth[i], eps())
    for i in 1:nrow(edges)
]

# ============================================================
# 6) DEGREE-PRESERVING NULL MODEL
# ============================================================
function edge_swap(preds, preys, rng; nsweeps=10)
    preds = copy(preds); preys = copy(preys)
    E = length(preds)
    for _ in 1:(nsweeps*E)
        i,j = rand(rng,1:E,2)
        i==j && continue
        preds[i]==preds[j] && continue
        preys[i]==preys[j] && continue
        preds[i],preds[j] = preds[j],preds[i]
    end
    return preds, preys
end

function null_stats(edges, statfn; nperm=1000, nsweeps=10)
    rng = MersenneTwister(1)
    vals = Float64[]
    for _ in 1:nperm
        p2,q2 = edge_swap(edges.pred, edges.prey, rng; nsweeps=nsweeps)
        df = copy(edges)
        df.pred .= p2; df.prey .= q2
        push!(vals, statfn(df))
    end
    return vals
end

# ============================================================
# 7) STATISTICS
# ============================================================
stats = Dict(
    :mean_overlap_pred => df -> mean(df.overlap_pred),
    :mean_absdiff_Tmid => df -> mean(abs.(df.pred_Tmid .- df.prey_Tmid)),
    :mean_absdiff_Tmax => df -> mean(abs.(df.pred_Tmax .- df.prey_Tmax)),
    :mean_absdiff_Tmin => df -> mean(abs.(df.pred_Tmin .- df.prey_Tmin)),
    :mean_absdiff_breadth => df -> mean(abs.(df.pred_breadth .- df.prey_breadth))
)

results = DataFrame(stat=String[], obs=Float64[],
                    null_mean=Float64[], null_sd=Float64[],
                    z=Float64[], p=Float64[])

for (k,f) in stats
    obs = f(edges)
    null = null_stats(edges,f)
    μ,σ = mean(null), std(null)
    z   = (obs-μ)/(σ>0 ? σ : NaN)
    p   = k==:mean_overlap_pred ? mean(null .>= obs) : mean(null .<= obs)
    push!(results,(string(k),obs,μ,σ,z,p))
end

CSV.write(joinpath(OUTDIR,"null_test_summary.csv"),results)

@info "DONE"
show(results, allrows=true, allcols=true)

# ============================================================
# 8) PLOTS — EDGE-LEVEL & NODE-LEVEL
# ============================================================

function fit_line(x::Vector{Float64}, y::Vector{Float64})
    x̄ = mean(x); ȳ = mean(y)
    denom = sum((x .- x̄).^2)
    denom == 0 && return (ȳ, 0.0)
    b = sum((x .- x̄) .* (y .- ȳ)) / denom
    a = ȳ - b*x̄
    return a, b
end

function scatter_with_fit(x, y; title="", xlabel="", ylabel="", outpath="")
    x = Float64.(x); y = Float64.(y)
    good = isfinite.(x) .& isfinite.(y)
    x = x[good]; y = y[good]
    n = length(x)

    r = (n ≥ 3) ? cor(x,y) : NaN
    has_fit = (n ≥ 3 && std(x) > 0)
    a,b = has_fit ? fit_line(x,y) : (NaN,NaN)

    fig = Figure(size=(900,650))
    ax = Axis(fig[1,1],
        title = string(title, " (r=", round(r,digits=2), ", n=", n, ")"),
        xlabel = xlabel,
        ylabel = ylabel
    )

    scatter!(ax, x, y)
    if has_fit
        xs = range(minimum(x), maximum(x), length=100)
        lines!(ax, xs, a .+ b .* xs, linewidth=2)
    end

    outpath != "" && save(outpath, fig)
    display(fig)
end

# ----------------------------
# Edge-level scatters
# ----------------------------
edge_pairs = [
    (:pred_Tmid, :prey_Tmid, "Edge-level: predator vs prey Tmid"),
    (:pred_Tmax, :prey_Tmax, "Edge-level: predator vs prey Tmax"),
    (:pred_Tmin, :prey_Tmin, "Edge-level: predator vs prey Tmin"),
    (:pred_breadth, :prey_breadth, "Edge-level: predator vs prey breadth")
]

for (x,y,ttl) in edge_pairs
    scatter_with_fit(edges[!,x], edges[!,y];
        title=ttl,
        xlabel=String(x),
        ylabel=String(y),
        outpath=joinpath(OUTDIR,"edge_$(x)_vs_$(y).png")
    )
end

# ----------------------------
# Predator vs mean(prey)
# ----------------------------
pred_level = combine(groupby(edges,:pred),
    :pred_Tmin => first => :pred_Tmin,
    :pred_Tmax => first => :pred_Tmax,
    :pred_Tmid => first => :pred_Tmid,
    :pred_breadth => first => :pred_breadth,
    :prey_Tmin => mean => :mean_prey_Tmin,
    :prey_Tmax => mean => :mean_prey_Tmax,
    :prey_Tmid => mean => :mean_prey_Tmid,
    :prey_breadth => mean => :mean_prey_breadth,
    nrow => :nprey
)

node_pairs = [
    (:pred_Tmid, :mean_prey_Tmid, "Predator vs mean(prey): Tmid"),
    (:pred_Tmax, :mean_prey_Tmax, "Predator vs mean(prey): Tmax"),
    (:pred_Tmin, :mean_prey_Tmin, "Predator vs mean(prey): Tmin"),
    (:pred_breadth, :mean_prey_breadth, "Predator vs mean(prey): breadth")
]

for (x,y,ttl) in node_pairs
    scatter_with_fit(pred_level[!,x], pred_level[!,y];
        title=ttl,
        xlabel=String(x),
        ylabel=String(y),
        outpath=joinpath(OUTDIR,"node_$(x)_vs_$(y).png")
    )
end

# ----------------------------
# Observed overlap distribution
# ----------------------------
begin
    fig = Figure(size=(800,600))
    ax = Axis(fig[1,1],
        title="Observed directed overlap (predator-side)",
        xlabel="Overlap",
        ylabel="Count"
    )
    hist!(ax, edges.overlap_pred, bins=30)
    save(joinpath(OUTDIR,"observed_overlap_hist.png"), fig)
    display(fig)
end
