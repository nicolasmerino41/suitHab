#!/usr/bin/env julia
# ============================================================
# ThermoFresh × ThermoTol × GloBI fish metaweb
# FULL METRIC PIPELINE
#
# • Predators: Actinopteri only
# • Prey: any ThermoFresh species
# • Metrics: ALL metrics present in thermtol_comb_final.csv
# • Outputs:
#     outputs_empirical_thermofresh_globi_fish/<metric>/
# ============================================================

using CSV, DataFrames, Statistics, Random
using CairoMakie

# ============================================================
# 0) PATHS
# ============================================================

script_dir   = @__DIR__
project_root = script_dir

thermtol_csv = joinpath(project_root, "thermtol_comb_final.csv")
metaweb_csv  = joinpath(project_root,
                        "thermofresh_globi_metaweb_fish_predators.csv")

OUTROOT = joinpath(project_root, "outputs_empirical_thermofresh_globi_fish")
isdir(OUTROOT) || mkpath(OUTROOT)

# ============================================================
# 1) CANONICAL TAXON
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
    lowercase(species) in ("sp","spp") && return ""

    return uppercase(genus[1]) * lowercase(genus[2:end]) * " " * lowercase(species)
end

# ============================================================
# 2) LOAD DATA
# ============================================================
tf_df = CSV.read(thermtol_csv, DataFrame; missingstring="", ntasks=1, pool=false)
mw    = CSV.read(metaweb_csv, DataFrame; missingstring="", ntasks=1, pool=false)

mw2 = DataFrame(
    pred = canon_taxon.(String.(mw.predator)),
    prey = canon_taxon.(String.(mw.prey))
)

# ============================================================
# 3) METRIC DISCOVERY
# ============================================================
metrics = sort(unique(String.(skipmissing(tf_df.metric))))
@info "Metrics detected" metrics

# ============================================================
# 4) HELPERS
# ============================================================

as_float(x) =
    x === missing ? NaN :
    x isa Real ? Float64(x) :
    x isa AbstractString ? (y=tryparse(Float64,x); y===nothing ? NaN : y) :
    NaN

function edge_swap(preds, rng; nsweeps=10)
    preds = copy(preds)
    E = length(preds)
    for _ in 1:(nsweeps * E)
        i,j = rand(rng,1:E,2)
        i==j && continue
        preds[i], preds[j] = preds[j], preds[i]
    end
    return preds
end

# ============================================================
# 5) MAIN METRIC LOOP
# ============================================================
pearson_r(x,y) = (length(x) ≥ 3 && std(x)>0 && std(y)>0) ? cor(x,y) : NaN

function fit_line(x,y)
    b = cov(x,y)/var(x)
    a = mean(y) - b*mean(x)
    return a,b
end

function scatter_with_fit(x,y,title,xlab,ylab,file)
    good = isfinite.(x) .& isfinite.(y)
    sum(good) < 2 && return

    x = x[good]; y = y[good]

    r = pearson_r(x,y)
    a,b = fit_line(x,y)

    fig = Figure(size=(900,570))
    ax = Axis(fig[1,1], title=title, xlabel=xlab, ylabel=ylab)
    scatter!(ax,x,y)

    xs = range(minimum(x), maximum(x), length=200)
    lines!(ax, xs, a .+ b .* xs)

    Label(fig[2,1:2], "r=$(round(r,digits=2)), n=$(length(x))", fontsize=13)
    save(file, fig)
    display(fig)
end

for METRIC in metrics

    @info "Running metric" METRIC
    outdir = joinpath(OUTROOT, METRIC)
    isdir(outdir) || mkpath(outdir)

    # ---------------------------
    # Metric lookup
    # ---------------------------
    df = tf_df[String.(tf_df.metric) .== METRIC, :]
    isempty(df) && continue

    canon_or_missing(x) = x === missing ? missing : canon_taxon(String(x))

    df = transform(
        df,
        :species => (s -> canon_or_missing.(s)) => :canon,
        :tol     => (v -> as_float.(v))          => :value
    )

    df = df[
        (.!ismissing.(df.canon)) .&
        (.!isempty.(coalesce.(df.canon, ""))) .&
        isfinite.(df.value),
        :
    ]

    metric_sum = combine(groupby(df,:canon), :value => mean => :value)
    lookup = Dict(r.canon => r.value for r in eachrow(metric_sum))

    # ---------------------------
    # Map metaweb
    # ---------------------------
    edges0 = mw2[haskey.(Ref(lookup), mw2.pred) .&
                 haskey.(Ref(lookup), mw2.prey), :]

    nrow(edges0) < 5 && @warn "Few edges for metric" METRIC nrow(edges0)

    # ---------------------------
    # Edge table
    # ---------------------------
    E = nrow(edges0)
    pred_val = [lookup[p] for p in edges0.pred]
    prey_val = [lookup[p] for p in edges0.prey]
    absdiff  = abs.(pred_val .- prey_val)

    edges = DataFrame(
        pred = edges0.pred,
        prey = edges0.prey,
        pred_val = pred_val,
        prey_val = prey_val,
        absdiff = absdiff
    )

    # ---------------------------
    # Node-level table
    # ---------------------------
    edges_clean = edges[
        isfinite.(edges.pred_val) .&
        isfinite.(edges.prey_val),
        :
    ]

    pred_grp = groupby(edges_clean, :pred)
    pred_level = combine(pred_grp,
        :pred_val => mean => :pred_val,
        :prey_val => mean => :mean_prey_val,
        nrow => :nprey
    )

    CSV.write(joinpath(outdir,"predator_level_$(METRIC).csv"), pred_level)

    # ---------------------------
    # Null model
    # ---------------------------
    rng = MersenneTwister(1)
    base_preds = collect(edges.pred)
    base_preys = collect(edges.prey)

    function null_mean_absdiff(edges, lookup; nperm=1000, nsweeps=10, seed=1)
        rng = MersenneTwister(seed)
        base_preds = collect(edges.pred)
        base_preys = collect(edges.prey)
        E = length(base_preds)

        vals = Float64[]
        sizehint!(vals, nperm)

        for _ in 1:nperm
            p2 = edge_swap(base_preds, rng; nsweeps=nsweeps)
            s = 0.0
            for i in 1:E
                s += abs(lookup[p2[i]] - lookup[base_preys[i]])
            end
            push!(vals, s / E)
        end
        return vals
    end

    null_vals = null_mean_absdiff(edges, lookup)

    obs = mean(edges.absdiff)
    μ0,σ0 = mean(null_vals), std(null_vals)
    z = (obs-μ0)/(σ0>0 ? σ0 : NaN)
    p = mean(null_vals .<= obs)

    # ---------------------------
    # PLOTS
    # ---------------------------
    clean = filter(isfinite, null_vals)

    μ0 = isempty(clean) ? NaN : mean(clean)
    σ0 = isempty(clean) ? NaN : std(clean)
    
    scatter_with_fit(
        edges.pred_val, edges.prey_val,
        "Edge-level: predator vs prey ($METRIC)",
        "Predator $METRIC",
        "Prey $METRIC",
        joinpath(outdir, "edge_pred_vs_prey_$(METRIC).png")
    )

    scatter_with_fit(
        pred_level.pred_val,
        pred_level.mean_prey_val,
        "Node-level: predator vs mean(prey) ($METRIC)",
        "Predator $METRIC",
        "Mean prey $METRIC",
        joinpath(outdir, "node_pred_vs_meanprey_$(METRIC).png")
    )

    fig = Figure(size=(900,570))
    ax = Axis(fig[1,1],
        title="Null: mean |Δ$METRIC|",
        xlabel="Mean |Δ$METRIC|",
        ylabel="Count"
    )
    clean = filter(isfinite, null_vals)
    isempty(clean) || hist!(ax, clean, bins=30)

    vlines!(ax,[obs])
    save(joinpath(outdir,"null_mean_absdiff.png"), fig)
    display(fig)
    # ---------------------------
    # SUMMARY
    # ---------------------------
    summary = DataFrame(
        metric = METRIC,
        obs = obs,
        null_mean = μ0,
        null_sd = σ0,
        z = z,
        p_one_sided = p,
        edges = nrow(edges),
        predators = length(unique(edges.pred))
    )

    CSV.write(joinpath(outdir,"summary_$(METRIC).csv"), summary)
end

@info "ALL METRICS COMPLETE"
