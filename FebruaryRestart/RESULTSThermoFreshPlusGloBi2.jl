# ============================================================
# ThermoFresh + ThermoTol + GloBI fish metaweb
# Metric-agnostic analysis of predator–prey thermal alignment
#
# CHANGE ONLY THE METRIC NAME BELOW
# ============================================================

##############################
# 0) USER CHOICE (ONLY CHANGE THIS)
##############################
METRIC = "lt50"   # "ctmax", "ctmin", "lt50", "ltmax", "ltmin"

##############################
# 1) DEPENDENCIES
##############################
using CSV, DataFrames, Statistics, Random
using CairoMakie

##############################
# 2) PATHS
##############################
script_dir   = @__DIR__
project_root = script_dir

thermtol_csv = joinpath(project_root, "FebruaryRestart", "thermtol_comb_final.csv")
metaweb_csv  = joinpath(project_root, "FebruaryRestart", "thermofresh_globi_metaweb_fish_predators.csv")

outdir = joinpath(project_root, "outputs_empirical_thermofresh_globi_fish_$(METRIC)")
isdir(outdir) || mkpath(outdir)

##############################
# 3) LOAD DATA
##############################
tf_df = CSV.read(thermtol_csv, DataFrame; missingstring="", ntasks=1, pool=false)
mw    = CSV.read(metaweb_csv, DataFrame; missingstring="", ntasks=1, pool=false)

##############################
# 4) CANONICAL TAXON FUNCTION
##############################
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

    return uppercase(genus[1]) * lowercase(genus[2:end]) * " " * lowercase(species)
end

##############################
# 5) THERMAL METRIC LOOKUP
##############################
species_col = :species
metric_col  = :metric
value_col   = :tol

metric_df = tf_df[string.(tf_df[!, metric_col]) .== METRIC, :]
nrow(metric_df) == 0 && error("No rows found for metric = $METRIC")
canon_or_missing(x) =
    x === missing ? missing : canon_taxon(String(x))

as_float(x) =
    x === missing ? NaN :
    x isa Real ? Float64(x) :
    x isa AbstractString ? (y = tryparse(Float64, x); y === nothing ? NaN : y) :
    NaN

metric_df = transform(
    metric_df,
    species_col => (s -> canon_or_missing.(s)) => :canon,
    value_col   => (v -> as_float.(v))          => :value
)

metric_df = metric_df[
    (.!ismissing.(metric_df.canon)) .&
    (.!isempty.(metric_df.canon)) .&
    isfinite.(metric_df.value),
    :
]

metric_sum = combine(groupby(metric_df, :canon), :value => mean => :value)

metric_lookup = Dict(r.canon => r.value for r in eachrow(metric_sum))

@info "Metric lookup built" (metric=METRIC, species=length(metric_lookup))

##############################
# 6) MAP META-WEB
##############################
mw2 = DataFrame(
    pred = canon_taxon.(String.(mw.predator)),
    prey = canon_taxon.(String.(mw.prey))
)

filter!(r -> haskey(metric_lookup, r.pred) && haskey(metric_lookup, r.prey), mw2)

@info "Matched edges" nrow(mw2)

##############################
# 7) EDGE TABLE
##############################
E = nrow(mw2)

pred_val = Vector{Float64}(undef, E)
prey_val = Vector{Float64}(undef, E)
absdiff  = Vector{Float64}(undef, E)

for i in 1:E
    p = metric_lookup[mw2.pred[i]]
    q = metric_lookup[mw2.prey[i]]
    pred_val[i] = p
    prey_val[i] = q
    absdiff[i]  = abs(p - q)
end

edges = DataFrame(
    pred = mw2.pred,
    prey = mw2.prey,
    pred_val = pred_val,
    prey_val = prey_val,
    absdiff  = absdiff
)

##############################
# 8) NODE-LEVEL TABLE
##############################
pred_grp = groupby(edges, :pred)
pred_level = combine(pred_grp,
    :pred_val => first => :pred_val,
    :prey_val => mean  => :mean_prey_val,
    nrow => :nprey
)

CSV.write(joinpath(outdir, "predator_level_$(METRIC).csv"), pred_level)

##############################
# 9) NULL MODEL (PREDATOR SWAP)
##############################
function edge_swap(preds, preys, rng; nsweeps=10)
    preds = copy(preds)
    E = length(preds)
    for _ in 1:(nsweeps * E)
        i, j = rand(rng, 1:E, 2)
        i == j && continue
        preds[i], preds[j] = preds[j], preds[i]
    end
    return preds
end

function null_mean_absdiff(edges, lookup; nperm=1000, nsweeps=10, seed=1)
    rng = MersenneTwister(seed)
    base_preds = collect(edges.pred)
    base_preys = collect(edges.prey)
    E = length(base_preds)

    vals = Float64[]
    sizehint!(vals, nperm)

    for _ in 1:nperm
        p2 = edge_swap(base_preds, base_preys, rng; nsweeps=nsweeps)
        s = 0.0
        for i in 1:E
            s += abs(lookup[p2[i]] - lookup[base_preys[i]])
        end
        push!(vals, s / E)
    end
    return vals
end

obs = mean(edges.absdiff)
null_vals = null_mean_absdiff(edges, metric_lookup)

μ0 = mean(null_vals)
σ0 = std(null_vals)
z  = (obs - μ0) / (σ0 > 0 ? σ0 : NaN)
p  = mean(null_vals .<= obs)

##############################
# 10) CORRELATION HELPERS
##############################
pearson_r(x,y) = (length(x) ≥ 3 && std(x)>0 && std(y)>0) ? cor(x,y) : NaN

function fit_line(x,y)
    b = cov(x,y)/var(x)
    a = mean(y) - b*mean(x)
    return a,b
end

##############################
# 11) PLOTS
##############################
scatter_with_fit(x,y,title,xlab,ylab,file) = begin
    r = pearson_r(x,y)
    a,b = fit_line(x,y)
    fig = Figure(size=(900,650))
    ax = Axis(fig[1,1], title=title, xlabel=xlab, ylabel=ylab)
    scatter!(ax,x,y)
    xs = range(minimum(x), maximum(x), length=200)
    lines!(ax, xs, a .+ b .* xs)
    Label(fig[2,1], "r=$r, n=$(length(x))", fontsize=13)
    save(file, fig)
end

scatter_with_fit(
    edges.pred_val, edges.prey_val,
    "Edge-level: predator vs prey ($METRIC)",
    "Predator $METRIC", "Prey $METRIC",
    joinpath(outdir, "edge_pred_vs_prey_$(METRIC).png")
)

scatter_with_fit(
    pred_level.pred_val, pred_level.mean_prey_val,
    "Node-level: predator vs mean(prey) ($METRIC)",
    "Predator $METRIC", "Mean prey $METRIC",
    joinpath(outdir, "node_pred_vs_meanprey_$(METRIC).png")
)

##############################
# 12) NULL HISTOGRAM
##############################
fig = Figure(size=(900,650))
ax = Axis(fig[1,1], title="Null distribution: mean |Δ$METRIC|",
          xlabel="Mean |Δ$METRIC|", ylabel="Count")
hist!(ax, null_vals, bins=30)
vlines!(ax, [obs], linewidth=3)
Label(fig[2,1],
    "obs=$(obs) | null=$(μ0)±$(σ0) | z=$(z) | p=$(p)",
    fontsize=13
)
save(joinpath(outdir, "null_mean_absdiff_$(METRIC).png"), fig)
display(fig)
##############################
# 13) SUMMARY
##############################
summary = DataFrame(
    metric = METRIC,
    obs_mean_absdiff = obs,
    null_mean = μ0,
    null_sd = σ0,
    z = z,
    p_one_sided = p,
    edges = nrow(edges),
    predators = length(unique(edges.pred))
)

CSV.write(joinpath(outdir, "summary_$(METRIC).csv"), summary)

@info "DONE" summary
