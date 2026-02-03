# ============================================================
# Empirical upgrade script (FIXED)
# ============================================================

using CSV, DataFrames, Statistics, Random
using CairoMakie

# ---------------------------
# 0) PATHS
# ---------------------------

script_dir    = @__DIR__
project_root  = dirname(script_dir)
data_dir      = joinpath(project_root, "data")

globtherm_csv = joinpath(project_root, "FebruaryRestart", "GlobalTherm_upload_02_11_17.csv")
pairwise_file = joinpath(data_dir, "TetraEU_pairwise_interactions.csv")

outdir = joinpath(project_root, "outputs_empirical")
isdir(outdir) || mkpath(outdir)

# ---------------------------
# 1) LOAD DATA
# ---------------------------

gt = CSV.read(globtherm_csv, DataFrame; missingstring="", ntasks=1, pool=false)

raw_web = CSV.read(pairwise_file, DataFrame)
web = DataFrame(
    predator = string.(raw_web.sourceTaxonName),
    prey     = string.(raw_web.targetTaxonName)
)

# ---------------------------
# 2) CANONICAL TAXON
# ---------------------------

function canon_taxon(s::AbstractString)::String
    t = replace(String(s), r"[_(),\[\]]" => " ")
    t = replace(t, r"\s+" => " ")
    t = strip(t)
    parts = split(t)
    length(parts) < 2 && return ""

    genus   = String(filter(isletter, parts[1]))
    species = String(filter(isletter, parts[2]))

    isempty(genus) || isempty(species) && return ""
    lowercase(species) in ("sp","spp") && return ""

    genus2   = uppercase(genus[1]) * lowercase(genus[2:end])
    species2 = lowercase(species)

    return "$genus2 $species2"
end

@inline function safe_get(v, i)
    (1 ≤ i ≤ length(v) && isassigned(v, i)) ? v[i] : missing
end

function gt_canon_name(gt::DataFrame, i::Int)::String
    g = safe_get(gt[!, :Genus], i)
    s = safe_get(gt[!, :Species], i)
    g === missing || s === missing ? "" : canon_taxon(string(g," ",s))
end

# ---------------------------
# 3) BUILD THERMAL NICHE TABLE
# ---------------------------

function pick_first_numeric(vals...)
    for v in vals
        v === missing && continue
        v isa Real && return Float64(v)
        v isa AbstractString && (x=tryparse(Float64,v); x!==nothing && return x)
    end
    return NaN
end

gt2 = copy(gt)
gt2.canon = [gt_canon_name(gt2,i) for i in 1:nrow(gt2)]

rows = Tuple{String,Float64,Float64}[]
for i in 1:nrow(gt2)
    cname = gt2.canon[i]
    isempty(cname) && continue
    tmin = pick_first_numeric(gt2[i,:tmin], gt2[i,:tmin_2])
    tmax = pick_first_numeric(gt2[i,:Tmax], gt2[i,:Tmax_2])
    isfinite(tmin) && isfinite(tmax) && tmax>tmin && push!(rows,(cname,tmin,tmax))
end

tmp = DataFrame(rows, [:canon,:Tmin,:Tmax])

gt_sum = combine(groupby(tmp,:canon),
    :Tmin => mean => :Tmin,
    :Tmax => mean => :Tmax)

gt_sum.Tmid    = (gt_sum.Tmin .+ gt_sum.Tmax)./2
gt_sum.breadth = gt_sum.Tmax .- gt_sum.Tmin

gt_lookup = Dict(r.canon => (Tmin=r.Tmin,Tmax=r.Tmax,Tmid=r.Tmid,breadth=r.breadth)
                 for r in eachrow(gt_sum))

# ---------------------------
# 4) MAP META-WEB
# ---------------------------

web2 = DataFrame(
    pred = canon_taxon.(web.predator),
    prey = canon_taxon.(web.prey)
)

filter!(r -> haskey(gt_lookup,r.pred) && haskey(gt_lookup,r.prey), web2)

# ---------------------------
# 5) EDGE-LEVEL METRICS
# ---------------------------

function interval_overlap(a1,a2,b1,b2)
    max(0.0, min(a2,b2)-max(a1,b1))
end

E = nrow(web2)

ov_pred = Float64[]
sizehint!(ov_pred,E)

for r in eachrow(web2)
    p = gt_lookup[r.pred]
    q = gt_lookup[r.prey]
    ov = interval_overlap(p.Tmin,p.Tmax,q.Tmin,q.Tmax)
    push!(ov_pred, ov / max(p.breadth, eps()))
end

edges = DataFrame(pred=web2.pred, prey=web2.prey, overlap_pred=ov_pred)

obs_stat = mean(edges.overlap_pred)
@info "Observed mean overlap" obs_stat

# ---------------------------
# 6) DEGREE-PRESERVING NULL
# ---------------------------

function edge_swap(edges::DataFrame, rng; nsweeps=10)
    preds = copy(edges.pred)
    preys = copy(edges.prey)
    E = length(preds)

    for _ in 1:(nsweeps*E)
        i,j = rand(rng,1:E,2)
        i==j && continue
        preds[i]==preds[j] && continue
        preys[i]==preys[j] && continue
        preds[i],preds[j] = preds[j],preds[i]
    end
    return DataFrame(pred=preds, prey=preys)
end

function edge_swap_null_stats(edges::DataFrame;
                              nperm=1000, nsweeps=10, seed=1)

    rng = MersenneTwister(seed)
    vals = Float64[]
    sizehint!(vals,nperm)

    for _ in 1:nperm
        rw = edge_swap(edges,rng;nsweeps=nsweeps)
        ovs = Float64[]
        for r in eachrow(rw)
            p = gt_lookup[r.pred]
            q = gt_lookup[r.prey]
            ov = interval_overlap(p.Tmin,p.Tmax,q.Tmin,q.Tmax)
            push!(ovs, ov/max(p.breadth,eps()))
        end
        push!(vals, mean(ovs))
    end
    return vals
end

@info "Running null model"
null_vals = edge_swap_null_stats(edges; nperm=1000, nsweeps=10)

@info "Null summary" (mean=mean(null_vals), sd=std(null_vals))

# ------------------------------------------------------------
# 8) NULL TEST SUMMARY
# ------------------------------------------------------------
null_stats = null_vals   # rename for clarity

μ0 = mean(null_stats)
σ0 = std(null_stats)
z  = (obs_stat - μ0) / (σ0 > 0 ? σ0 : NaN)

# one-sided p-value: enrichment (observed > null)
pval = mean(null_stats .>= obs_stat)

@info "Null test (mean overlap_pred)" (
    obs = obs_stat,
    null_mean = μ0,
    null_sd = σ0,
    z = z,
    p_one_sided = pval
)

# ------------------------------------------------------------
# 9) PLOTS
# ------------------------------------------------------------

# 9.1 Observed directed overlap distribution
begin
    fig_obs = Figure(size = (1200, 720))
    ax_obs = Axis(
        fig_obs[1, 1],
        title  = "Observed predator–prey thermal alignment",
        xlabel = "Directed overlap (predator-side)",
        ylabel = "Count"
    )
    hist!(ax_obs, edges.overlap_pred, bins = 30)
    display(fig_obs)
end

save(joinpath(outdir, "observed_overlap_distribution.png"), fig_obs)

# 9.2 Null distribution of mean overlap
begin
    fig_null = Figure(size = (800, 600))
    ax_null = Axis(
        fig_null[1, 1],
        title  = "Null distribution: mean directed overlap",
        xlabel = "Mean overlap under degree-preserving null",
        ylabel = "Count"
    )
    hist!(ax_null, null_stats, bins = 30)
    vlines!(ax_null, [obs_stat], linewidth = 3, linestyle = :dash)

    Label(
        fig_null[2, 1],
        "Observed = $(round(obs_stat, digits=4)) | " *
        "Null mean = $(round(μ0, digits=4)) ± $(round(σ0, digits=4)) | " *
        "z = $(round(z, digits=2)) | p(one-sided) = $(round(pval, digits=4))",
        fontsize = 12
    )
    display(fig_null)
end
save(joinpath(outdir, "null_mean_overlap.png"), fig_null)

@info "Saved plots to" outdir

println("\n=== DONE ===")
println("Observed mean overlap_pred = ", obs_stat)
println("Null mean = ", μ0, "  sd = ", σ0)
println("z-score = ", z)
println("One-sided p-value = ", pval)
