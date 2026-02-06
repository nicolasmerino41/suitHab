# ============================================================
# Empirical upgrade script — MULTI-STAT NULL TEST + AUTO PLOTS
#   Runs null tests for:
#     - mean directed overlap_pred
#     - mean abs mismatch for Tmin/Tmax/Tmid/breadth
#   And plots ALL outputs automatically.
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
    isfinite(tmin) && isfinite(tmax) && tmax > tmin && push!(rows,(cname,tmin,tmax))
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

@info "Matched interactions" nrow(web2)

# ---------------------------
# 5) EDGE-LEVEL TABLE (traits + overlap)
# ---------------------------

function interval_overlap(a1,a2,b1,b2)
    max(0.0, min(a2,b2)-max(a1,b1))
end

E = nrow(web2)

ov_pred = Float64[]
pred_Tmin = Float64[]; pred_Tmax = Float64[]; pred_Tmid = Float64[]; pred_breadth = Float64[]
prey_Tmin = Float64[]; prey_Tmax = Float64[]; prey_Tmid = Float64[]; prey_breadth = Float64[]

sizehint!(ov_pred,E)
sizehint!(pred_Tmin,E); sizehint!(pred_Tmax,E); sizehint!(pred_Tmid,E); sizehint!(pred_breadth,E)
sizehint!(prey_Tmin,E); sizehint!(prey_Tmax,E); sizehint!(prey_Tmid,E); sizehint!(prey_breadth,E)

for r in eachrow(web2)
    p = gt_lookup[r.pred]
    q = gt_lookup[r.prey]

    push!(pred_Tmin, p.Tmin); push!(pred_Tmax, p.Tmax); push!(pred_Tmid, p.Tmid); push!(pred_breadth, p.breadth)
    push!(prey_Tmin, q.Tmin); push!(prey_Tmax, q.Tmax); push!(prey_Tmid, q.Tmid); push!(prey_breadth, q.breadth)

    ov = interval_overlap(p.Tmin,p.Tmax,q.Tmin,q.Tmax)
    push!(ov_pred, ov / max(p.breadth, eps()))
end

edges = DataFrame(
    pred=web2.pred, prey=web2.prey,
    overlap_pred=ov_pred,
    pred_Tmin=pred_Tmin, pred_Tmax=pred_Tmax, pred_Tmid=pred_Tmid, pred_breadth=pred_breadth,
    prey_Tmin=prey_Tmin, prey_Tmax=prey_Tmax, prey_Tmid=prey_Tmid, prey_breadth=prey_breadth
)

# ---------------------------
# 6) DEGREE-PRESERVING NULL (swap on predator identities)
#     IMPORTANT: because your swap shuffles predators among edges,
#     the edge table's pred_* columns must be recomputed per perm.
# ---------------------------

function edge_swap_predprey(preds::Vector{String}, preys::Vector{String}, rng; nsweeps=10)
    preds = copy(preds)
    preys = copy(preys)
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

# Build a lightweight permuted edge table (only what we need for stats)
function build_perm_edges(preds::Vector{String}, preys::Vector{String}, gt_lookup)::DataFrame
    E = length(preds)

    ov = Vector{Float64}(undef, E)
    pTmin = Vector{Float64}(undef, E)
    pTmax = Vector{Float64}(undef, E)
    pTmid = Vector{Float64}(undef, E)
    pB    = Vector{Float64}(undef, E)
    qTmin = Vector{Float64}(undef, E)
    qTmax = Vector{Float64}(undef, E)
    qTmid = Vector{Float64}(undef, E)
    qB    = Vector{Float64}(undef, E)

    for i in 1:E
        p = gt_lookup[preds[i]]
        q = gt_lookup[preys[i]]
        pTmin[i]=p.Tmin; pTmax[i]=p.Tmax; pTmid[i]=p.Tmid; pB[i]=p.breadth
        qTmin[i]=q.Tmin; qTmax[i]=q.Tmax; qTmid[i]=q.Tmid; qB[i]=q.breadth
        inter = interval_overlap(p.Tmin,p.Tmax,q.Tmin,q.Tmax)
        ov[i] = inter / max(p.breadth, eps())
    end

    return DataFrame(
        pred=preds, prey=preys,
        overlap_pred=ov,
        pred_Tmin=pTmin, pred_Tmax=pTmax, pred_Tmid=pTmid, pred_breadth=pB,
        prey_Tmin=qTmin, prey_Tmax=qTmax, prey_Tmid=qTmid, prey_breadth=qB
    )
end

function null_stats_for(edges::DataFrame, gt_lookup, stat_fn;
                        nperm=1000, nsweeps=10, seed=1)

    rng = MersenneTwister(seed)
    base_preds = collect(edges.pred)
    base_preys = collect(edges.prey)

    vals = Float64[]
    sizehint!(vals, nperm)

    for _ in 1:nperm
        p2, q2 = edge_swap_predprey(base_preds, base_preys, rng; nsweeps=nsweeps)
        df = build_perm_edges(p2, q2, gt_lookup)
        push!(vals, stat_fn(df))
    end

    return vals
end

# ---------------------------
# 7) DEFINE *ALL* STATS (observed + null)
# ---------------------------

stats = Dict{Symbol,NamedTuple{(:title,:xlabel,:obs_edge,:stat_fn),Tuple{String,String,Symbol,Function}}}()

# obs_edge tells us which edge-level column to histogram for "observed distribution"
# (for mismatch stats we’ll create the edge-level mismatch vector on the fly)

stats[:mean_overlap_pred] = (
    title="Mean directed overlap (predator-side)",
    xlabel="Mean directed overlap",
    obs_edge=:overlap_pred,
    stat_fn = df -> mean(df.overlap_pred)
)

stats[:mean_abs_diff_Tmin] = (
    title="Mean |Tmin(pred) - Tmin(prey)|",
    xlabel="Mean abs difference (°C)",
    obs_edge=Symbol("absdiff_Tmin"),
    stat_fn = df -> mean(abs.(df.pred_Tmin .- df.prey_Tmin))
)

stats[:mean_abs_diff_Tmax] = (
    title="Mean |Tmax(pred) - Tmax(prey)|",
    xlabel="Mean abs difference (°C)",
    obs_edge=Symbol("absdiff_Tmax"),
    stat_fn = df -> mean(abs.(df.pred_Tmax .- df.prey_Tmax))
)

stats[:mean_abs_diff_Tmid] = (
    title="Mean |Tmid(pred) - Tmid(prey)|",
    xlabel="Mean abs difference (°C)",
    obs_edge=Symbol("absdiff_Tmid"),
    stat_fn = df -> mean(abs.(df.pred_Tmid .- df.prey_Tmid))
)

stats[:mean_abs_diff_breadth] = (
    title="Mean |breadth(pred) - breadth(prey)|",
    xlabel="Mean abs difference (°C)",
    obs_edge=Symbol("absdiff_breadth"),
    stat_fn = df -> mean(abs.(df.pred_breadth .- df.prey_breadth))
)

# observed edge-level mismatch vectors (computed once)
obs_edgevals = Dict{Symbol,Vector{Float64}}()
obs_edgevals[:overlap_pred] = collect(edges.overlap_pred)
obs_edgevals[Symbol("absdiff_Tmin")] = abs.(edges.pred_Tmin .- edges.prey_Tmin)
obs_edgevals[Symbol("absdiff_Tmax")] = abs.(edges.pred_Tmax .- edges.prey_Tmax)
obs_edgevals[Symbol("absdiff_Tmid")] = abs.(edges.pred_Tmid .- edges.prey_Tmid)
obs_edgevals[Symbol("absdiff_breadth")] = abs.(edges.pred_breadth .- edges.prey_breadth)

# ---------------------------
# 8) PLOT HELPERS
# ---------------------------

function pearson_r(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    (length(x) == length(y) && length(x) >= 3) || return NaN
    return cor(collect(x), collect(y))
end

function fit_line(x::Vector{Float64}, y::Vector{Float64})
    x̄ = mean(x); ȳ = mean(y)
    denom = sum((x .- x̄).^2)
    denom == 0 && return (ȳ, 0.0)
    b = sum((x .- x̄) .* (y .- ȳ)) / denom
    a = ȳ - b*x̄
    return a, b
end

function scatter_with_fit(x, y; title="", xlabel="", ylabel="")
    x = Float64.(x); y = Float64.(y)
    good = isfinite.(x) .& isfinite.(y)
    x = x[good]; y = y[good]
    n = length(x)

    r = (n >= 3) ? pearson_r(x, y) : NaN

    fig = Figure(size=(900, 650))
    ax = Axis(fig[1,1],
        title = string(title, " (Pearson=", round(r, digits=2), ", n=", n, ")"),
        xlabel = xlabel,
        ylabel = ylabel
    )

    scatter!(ax, x, y, markersize=10)

    if n >= 3 && std(x) > 0
        a,b = fit_line(x,y)
        xs = range(minimum(x), maximum(x), length=100)
        lines!(ax, xs, a .+ b .* xs, linewidth=2)
    end

    return fig
end

function hist_plot(vals; title="", xlabel="", ylabel="Count", bins=30)
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1,1], title=title, xlabel=xlabel, ylabel=ylabel)
    hist!(ax, vals, bins=bins)
    return fig
end

function null_hist_plot(null_vals, obs; title="", xlabel="", bins=30, footer="")
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1,1], title=title, xlabel=xlabel, ylabel="Count")
    hist!(ax, null_vals, bins=bins)
    vlines!(ax, [obs], linewidth=3, linestyle=:dash)
    if !isempty(footer)
        Label(fig[2,1], footer, fontsize=12)
    end
    return fig
end

# ---------------------------
# 9) PLOTS: SCATTERS (same as before, but auto)
# ---------------------------

pred_grp = groupby(edges, :pred)
pred_level = combine(pred_grp,
    :pred_Tmin     => first => :pred_Tmin,
    :pred_Tmax     => first => :pred_Tmax,
    :pred_Tmid     => first => :pred_Tmid,
    :pred_breadth  => first => :pred_breadth,
    :prey_Tmin     => mean  => :mean_prey_Tmin,
    :prey_Tmax     => mean  => :mean_prey_Tmax,
    :prey_Tmid     => mean  => :mean_prey_Tmid,
    :prey_breadth  => mean  => :mean_prey_breadth,
    nrow => :nprey
)

pred_pairs = [
    (:pred_breadth, :mean_prey_breadth, "Predator vs mean(prey): breadth", "Predator breadth", "Mean prey breadth", "pred_vs_meanprey_breadth.png"),
    (:pred_Tmid,    :mean_prey_Tmid,    "Predator vs mean(prey): Tmid",    "Predator Tmid",    "Mean prey Tmid",    "pred_vs_meanprey_Tmid.png"),
    (:pred_Tmax,    :mean_prey_Tmax,    "Predator vs mean(prey): Tmax",    "Predator Tmax",    "Mean prey Tmax",    "pred_vs_meanprey_Tmax.png"),
    (:pred_Tmin,    :mean_prey_Tmin,    "Predator vs mean(prey): Tmin",    "Predator Tmin",    "Mean prey Tmin",    "pred_vs_meanprey_Tmin.png"),
]

for (xcol,ycol,ttl,xlab,ylab,fname) in pred_pairs
    fig = scatter_with_fit(pred_level[!,xcol], pred_level[!,ycol];
        title=ttl, xlabel=xlab, ylabel=ylab
    )
    save(joinpath(outdir, fname), fig)
    display(fig)
end

edge_pairs = [
    (:pred_breadth, :prey_breadth, "Edge-level: predator breadth vs prey breadth", "Predator breadth", "Prey breadth", "edge_pred_breadth_vs_prey_breadth.png"),
    (:pred_Tmid,    :prey_Tmid,    "Edge-level: predator Tmid vs prey Tmid",       "Predator Tmid",    "Prey Tmid",    "edge_pred_Tmid_vs_prey_Tmid.png"),
    (:pred_Tmax,    :prey_Tmax,    "Edge-level: predator Tmax vs prey Tmax",       "Predator Tmax",    "Prey Tmax",    "edge_pred_Tmax_vs_prey_Tmax.png"),
    (:pred_Tmin,    :prey_Tmin,    "Edge-level: predator Tmin vs prey Tmin",       "Predator Tmin",    "Prey Tmin",    "edge_pred_Tmin_vs_prey_Tmin.png"),
]

for (xcol,ycol,ttl,xlab,ylab,fname) in edge_pairs
    fig = scatter_with_fit(edges[!,xcol], edges[!,ycol];
        title=ttl, xlabel=xlab, ylabel=ylab
    )
    save(joinpath(outdir, fname), fig)
    display(fig)
end

# ---------------------------
# 10) PLOTS + NULL TESTS FOR *ALL* STATS
# ---------------------------

results = DataFrame(
    stat = String[],
    obs = Float64[],
    null_mean = Float64[],
    null_sd = Float64[],
    z = Float64[],
    p_one_sided = Float64[]
)

for (k, meta) in stats
    statname = Symbol(k)

    # observed statistic
    obs = meta.stat_fn(edges)

    # null distribution
    @info "Running null model for" statname
    null_vals = null_stats_for(edges, gt_lookup, meta.stat_fn; nperm=1000, nsweeps=10, seed=1)

    μ0 = mean(null_vals)
    σ0 = std(null_vals)
    z  = (obs - μ0) / (σ0 > 0 ? σ0 : NaN)

    # one-sided p-value depends on direction:
    #   overlap: enrichment => obs > null is "interesting"
    #   mismatch: depletion => obs < null is "interesting"
    p = if statname == :mean_overlap_pred
        mean(null_vals .>= obs)  # enrichment
    else
        mean(null_vals .<= obs)  # smaller mismatch than null
    end

    push!(results, (string(statname), obs, μ0, σ0, z, p))

    # ---- plot observed edge-level distribution for that stat
    edgetag = meta.obs_edge
    vals = obs_edgevals[edgetag]
    fig_obs = hist_plot(vals;
        title = "Observed edge distribution: $(meta.title)",
        xlabel = (statname == :mean_overlap_pred ? "Directed overlap (predator-side)" : meta.xlabel),
        bins = 30
    )
    save(joinpath(outdir, "observed_edge_$(statname).png"), fig_obs)
    display(fig_obs)

    # ---- plot null distribution for the mean statistic
    footer =
        "Observed = $(round(obs, digits=4)) | " *
        "Null mean = $(round(μ0, digits=4)) ± $(round(σ0, digits=4)) | " *
        "z = $(round(z, digits=2)) | p(one-sided) = $(round(p, digits=4))"

    fig_null = null_hist_plot(null_vals, obs;
        title = "Null distribution: $(meta.title)",
        xlabel = "Mean under degree-preserving null",
        bins = 30,
        footer = footer
    )
    save(joinpath(outdir, "null_mean_$(statname).png"), fig_null)
    display(fig_null)
end

CSV.write(joinpath(outdir, "null_test_summary.csv"), results)

@info "Saved plots + summary to" outdir
println("\n=== DONE ===")
println("Wrote summary table: ", joinpath(outdir, "null_test_summary.csv"))
println("Stats:")
show(results, allrows=true, allcols=true)
println()

# ---------------------------
# 9B) NODE-LEVEL CORRELATIONS
#     Predator trait vs mean prey trait per predator
#     (accounts for predators having multiple preys)
# ---------------------------
function scatter_with_fit(x, y; title="", xlabel="", ylabel="")
    x = Float64.(x); y = Float64.(y)
    good = isfinite.(x) .& isfinite.(y)
    x = x[good]; y = y[good]
    n = length(x)

    r = (n >= 3) ? pearson_r(x, y) : NaN

    has_fit = (n >= 3 && std(x) > 0)
    a, b = has_fit ? fit_line(x, y) : (NaN, NaN)

    # put stats on a NEW LINE to avoid cropping
    stats_line = has_fit ?
        "r=$(round(r, digits=2)), n=$(n), y=$(round(a, digits=2)) + $(round(b, digits=2))x" :
        "r=$(round(r, digits=2)), n=$(n)"

    fig = Figure(size=(900, 650), figure_padding=(20, 20, 40, 20))  # extra top padding

    ax = Axis(fig[1,1],
        title  = string(title, "\n", stats_line),
        xlabel = xlabel,
        ylabel = ylabel
    )

    scatter!(ax, x, y, markersize=10)

    if has_fit
        xs = range(minimum(x), maximum(x), length=100)
        lines!(ax, xs, a .+ b .* xs, linewidth=2)
    end

    return fig
end

@info "Building predator-level (node-level) summaries"

# You already built pred_level above:
# pred_level = combine(groupby(edges, :pred), ...)

# Optional: exclude predators with very few prey (to avoid 1-prey noise)
min_prey_per_pred = 1 # set to 1 to include all without filtering
pred_level_filt = filter(row -> row.nprey >= min_prey_per_pred, pred_level)

@info "Predator-level table sizes" (
    all_preds = nrow(pred_level),
    filtered_preds = nrow(pred_level_filt),
    min_prey_per_pred = min_prey_per_pred
)

# Helper to run and save both "all" and "filtered" versions
function save_nodelevel_pair(df::DataFrame, xcol::Symbol, ycol::Symbol,
                             ttl::String, xlab::String, ylab::String,
                             outpath::String)
    fig = scatter_with_fit(df[!, xcol], df[!, ycol];
        title=ttl, xlabel=xlab, ylabel=ylab
    )
    save(outpath, fig)
    display(fig)
end

# Node-level pairs: predator trait vs mean prey trait (same trait)
node_pairs = [
    (:pred_Tmin,    :mean_prey_Tmin,    "Node-level: predator Tmin vs mean(prey Tmin)",    "Predator Tmin",    "Mean prey Tmin"),
    (:pred_Tmax,    :mean_prey_Tmax,    "Node-level: predator Tmax vs mean(prey Tmax)",    "Predator Tmax",    "Mean prey Tmax"),
    (:pred_Tmid,    :mean_prey_Tmid,    "Node-level: predator Tmid vs mean(prey Tmid)",    "Predator Tmid",    "Mean prey Tmid"),
    (:pred_breadth, :mean_prey_breadth, "Node-level: predator breadth vs mean(prey breadth)", "Predator breadth", "Mean prey breadth"),
]

# --- Plot for ALL predators
for (xcol, ycol, ttl, xlab, ylab) in node_pairs
    fname = "nodelevel_all_$(String(xcol))_vs_$(String(ycol)).png"
    save_nodelevel_pair(pred_level, xcol, ycol,
        ttl * " (all predators)",
        xlab, ylab,
        joinpath(outdir, fname)
    )
end

# --- Plot for FILTERED predators (nprey >= min_prey_per_pred)
if nrow(pred_level_filt) >= 3
    for (xcol, ycol, ttl, xlab, ylab) in node_pairs
        fname = "nodelevel_filt_nprey_ge_$(min_prey_per_pred)_$(String(xcol))_vs_$(String(ycol)).png"
        save_nodelevel_pair(pred_level_filt, xcol, ycol,
            ttl * " (nprey ≥ $(min_prey_per_pred))",
            xlab, ylab,
            joinpath(outdir, fname)
        )
    end
else
    @warn "Not enough predators after filtering to compute node-level correlations" nrow(pred_level_filt)
end

# --- Optional: export predator-level table for inspection
CSV.write(joinpath(outdir, "predator_level_traits.csv"), pred_level)
