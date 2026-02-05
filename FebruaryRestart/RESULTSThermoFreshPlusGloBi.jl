# ============================================================
# ThermoFresh names + ThermoTol combined (ctmax) + GloBI fish metaweb
#   - Uses ONLY metric == "ctmax"
#   - Tests predator–prey Tmax alignment:
#       * correlations (edge-level + node-level)
#       * null test on mean |Tmax(pred) - Tmax(prey)|
# ============================================================

using CSV, DataFrames, Statistics, Random
using CairoMakie

# ---------------------------
# 0) PATHS
# ---------------------------
script_dir   = @__DIR__
project_root = script_dir

# Species name table (you said: only for species names; optional for sanity checks)
tf_names_csv = joinpath(project_root, "FebruaryRestart", "thermtol_taxonomy_final.csv")

# Real thermal data (long format, contains column :metric)
thermtol_csv = joinpath(project_root, "FebruaryRestart", "thermtol_comb_final.csv")

# Your new fish-predator metaweb
metaweb_csv  = joinpath(project_root, "FebruaryRestart", "thermofresh_globi_metaweb_fish_predators.csv")

outdir = joinpath(project_root, "outputs_empirical_thermofresh_globi_fish_ctmax")
isdir(outdir) || mkpath(outdir)

# ---------------------------
# 1) LOAD DATA
# ---------------------------

tf_names = CSV.read(tf_names_csv, DataFrame; missingstring="", ntasks=1, pool=false)
tf_df    = CSV.read(thermtol_csv, DataFrame; missingstring="", ntasks=1, pool=false)
mw       = CSV.read(metaweb_csv, DataFrame; missingstring="", ntasks=1, pool=false)

("predator" in names(mw) && "prey" in names(mw)) || error("Metaweb must have columns: predator, prey")

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

    (isempty(genus) || isempty(species)) && return ""
    lowercase(species) in ("sp","spp") && return ""

    genus2   = uppercase(genus[1]) * lowercase(genus[2:end])
    species2 = lowercase(species)
    return "$genus2 $species2"
end

# ---------------------------
# 3) THERMAL CTMAX LOOKUP (metric == "ctmax")
# ---------------------------

# Detect likely columns in thermtol_comb_final
# We need:
#   - species name col
#   - metric col (you said it's literally called :metric)
#   - value col (often :value, :val, :estimate, :temp, etc.)

function get_col(df::DataFrame, candidates::Vector{Symbol})
    for c in candidates
        c in names(df) && return c
    end
    return nothing
end

species_col = "species"
metric_col  = "metric"
value_col   = "tol"

species_col === nothing && error("Couldn't find a species column in thermtol_comb_final.csv")
metric_col  === nothing && error("Couldn't find a metric column in thermtol_comb_final.csv (expected :metric)")
value_col   === nothing && error("Couldn't find a value column in thermtol_comb_final.csv")

@info "Columns used in thermtol_comb_final" (species_col=species_col, metric_col=metric_col, value_col=value_col)

# Filter to ctmax
ctmax = tf_df[string.(tf_df[!, metric_col]) .== "ctmax", :]
nrow(ctmax) == 0 && error("After filtering metric == \"ctmax\", there are 0 rows. Check spelling/case of metric values.")

# Canonical species + numeric Tmax
function as_float(x)
    x === missing && return NaN
    x isa Real && return Float64(x)
    x isa AbstractString && (y=tryparse(Float64, x); y !== nothing ? y : NaN)
    return NaN
end
to_string_or_missing(x) =
    x === missing ? missing : String(x)
canon_or_missing(x) =
    x === missing ? missing : canon_taxon(x)

ctmax = transform(ctmax,
    species_col => (s -> canon_or_missing.(to_string_or_missing.(s))) => :canon,
    value_col   => (v -> as_float.(v))                                 => :Tmax
)

# Keep valid
ctmax = ctmax[
    (.!ismissing.(ctmax.canon)) .&
    (.!isempty.(ctmax.canon)) .&
    isfinite.(ctmax.Tmax),
    :
]
nrow(ctmax) == 0 && error("No valid ctmax rows after parsing species/value.")

# If multiple measurements per species: average
ctmax_sum = combine(groupby(ctmax, :canon), :Tmax => mean => :Tmax)

tmax_lookup = Dict(r.canon => r.Tmax for r in eachrow(ctmax_sum))

@info "ctmax lookup built" (species=length(tmax_lookup), rows=nrow(ctmax_sum))

# ---------------------------
# 4) MAP META-WEB TO CTMAX
# ---------------------------
mw2 = DataFrame(
    pred = canon_taxon.(String.(mw.predator)),
    prey = canon_taxon.(String.(mw.prey))
)

filter!(r -> haskey(tmax_lookup, r.pred) && haskey(tmax_lookup, r.prey), mw2)

@info "Matched interactions (have ctmax for both)" nrow(mw2)

nrow(mw2) < 5 && @warn "Very few matched edges; results will be noisy." nrow(mw2)

# ---------------------------
# 5) EDGE TABLE (Tmax only)
# ---------------------------
E = nrow(mw2)
pred_Tmax = Vector{Float64}(undef, E)
prey_Tmax = Vector{Float64}(undef, E)
absdiff   = Vector{Float64}(undef, E)

for i in 1:E
    p = tmax_lookup[mw2.pred[i]]
    q = tmax_lookup[mw2.prey[i]]
    pred_Tmax[i] = p
    prey_Tmax[i] = q
    absdiff[i]   = abs(p - q)
end

edges = DataFrame(
    pred = mw2.pred,
    prey = mw2.prey,
    pred_Tmax = pred_Tmax,
    prey_Tmax = prey_Tmax,
    absdiff_Tmax = absdiff
)

# ---------------------------
# 6) DEGREE-PRESERVING NULL (swap predator identities)
# ---------------------------

function edge_swap_predprey(preds::Vector{String}, preys::Vector{String}, rng; nsweeps=10)
    preds = copy(preds)
    preys = copy(preys)
    E = length(preds)

    for _ in 1:(nsweeps * E)
        i, j = rand(rng, 1:E, 2)
        i == j && continue
        preds[i] == preds[j] && continue
        preys[i] == preys[j] && continue
        preds[i], preds[j] = preds[j], preds[i]
    end
    return preds, preys
end

function null_stats_mean_absdiff_Tmax(edges::DataFrame, tmax_lookup;
                                      nperm=1000, nsweeps=10, seed=1)
    rng = MersenneTwister(seed)
    base_preds = collect(edges.pred)
    base_preys = collect(edges.prey)
    E = length(base_preds)

    vals = Float64[]
    sizehint!(vals, nperm)

    for _ in 1:nperm
        p2, q2 = edge_swap_predprey(base_preds, base_preys, rng; nsweeps=nsweeps)
        # compute mean abs diff under perm
        s = 0.0
        for i in 1:E
            s += abs(tmax_lookup[p2[i]] - tmax_lookup[q2[i]])
        end
        push!(vals, s / E)
    end
    return vals
end

# Observed statistic
obs_mean_absdiff = mean(edges.absdiff_Tmax)

@info "Observed mean |ΔTmax|" obs_mean_absdiff

# Null distribution
@info "Running null model for mean |ΔTmax|"
null_vals = null_stats_mean_absdiff_Tmax(edges, tmax_lookup; nperm=1000, nsweeps=10, seed=1)

μ0 = mean(null_vals)
σ0 = std(null_vals)
z  = (obs_mean_absdiff - μ0) / (σ0 > 0 ? σ0 : NaN)

# One-sided p-value for *smaller mismatch than null* (alignment)
p_one = mean(null_vals .<= obs_mean_absdiff)

@info "Null test summary (mean |ΔTmax|)" (obs=obs_mean_absdiff, null_mean=μ0, null_sd=σ0, z=z, p_one_sided=p_one)

# ---------------------------
# 7) PLOT HELPERS (title + visible label)
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
    a = ȳ - b * x̄
    return a, b
end

function scatter_with_fit(x, y; title="", xlabel="", ylabel="", outfile=nothing)
    x = Float64.(x); y = Float64.(y)
    good = isfinite.(x) .& isfinite.(y)
    x = x[good]; y = y[good]
    n = length(x)

    r = (n >= 3) ? pearson_r(x, y) : NaN
    has_fit = (n >= 3 && std(x) > 0)
    a, b = has_fit ? fit_line(x, y) : (NaN, NaN)

    fig = Figure(size=(900, 650))
    ax = Axis(fig[1, 1],
        title  = title,
        xlabel = xlabel,
        ylabel = ylabel
    )

    scatter!(ax, x, y, markersize=10)

    if has_fit
        xs = range(minimum(x), maximum(x), length=200)
        lines!(ax, xs, a .+ b .* xs, linewidth=2)
    end

    # IMPORTANT: put values as a visible label (titles sometimes render but you asked "nowhere")
    txt = has_fit ?
        "r=$(round(r,digits=2)), n=$n, fit: y=$(round(a,digits=2)) + $(round(b,digits=2))x" :
        "r=$(round(r,digits=2)), n=$n"

    Label(fig[2, 1], txt, fontsize=13)

    outfile !== nothing && save(outfile, fig)
    display(fig)
    return fig
end

function hist_plot(vals; title="", xlabel="", ylabel="Count", bins=30, outfile=nothing)
    fig = Figure(size=(900, 650))
    ax = Axis(fig[1,1], title=title, xlabel=xlabel, ylabel=ylabel)
    hist!(ax, vals, bins=bins)
    outfile !== nothing && save(outfile, fig)
    display(fig)
    return fig
end

function null_hist_plot(null_vals, obs; title="", xlabel="", bins=30, footer="", outfile=nothing)
    fig = Figure(size=(900, 650))
    ax = Axis(fig[1,1], title=title, xlabel=xlabel, ylabel="Count")
    hist!(ax, null_vals, bins=bins)
    vlines!(ax, [obs], linewidth=3, linestyle=:dash)
    Label(fig[2,1], footer, fontsize=13)
    outfile !== nothing && save(outfile, fig)
    display(fig)
    return fig
end

# ---------------------------
# 8) PLOTS: EDGE-LEVEL CORRELATION (Tmax)
# ---------------------------

scatter_with_fit(
    edges.pred_Tmax, edges.prey_Tmax;
    title="Edge-level: predator Tmax vs prey Tmax",
    xlabel="Predator Tmax",
    ylabel="Prey Tmax",
    outfile=joinpath(outdir, "edge_pred_Tmax_vs_prey_Tmax.png")
)

# ---------------------------
# 9) PLOTS: NODE-LEVEL CORRELATION (Tmax)
#     predator Tmax vs mean(prey Tmax) per predator
# ---------------------------
pred_grp = groupby(edges, :pred)
pred_level = combine(pred_grp,
    :pred_Tmax => first => :pred_Tmax,
    :prey_Tmax => mean  => :mean_prey_Tmax,
    nrow => :nprey
)

# optional filtering (avoid 1-prey predators dominating)
min_prey_per_pred = 3
pred_level_filt = filter(r -> r.nprey >= min_prey_per_pred, pred_level)

@info "Node-level table sizes" (all=nrow(pred_level), filtered=nrow(pred_level_filt), min_prey=min_prey_per_pred)

scatter_with_fit(
    pred_level.pred_Tmax, pred_level.mean_prey_Tmax;
    title="Node-level: predator Tmax vs mean(prey Tmax) (all predators)",
    xlabel="Predator Tmax",
    ylabel="Mean prey Tmax",
    outfile=joinpath(outdir, "nodelevel_all_pred_Tmax_vs_meanprey_Tmax.png")
)

if nrow(pred_level_filt) >= 3
    scatter_with_fit(
        pred_level_filt.pred_Tmax, pred_level_filt.mean_prey_Tmax;
        title="Node-level: predator Tmax vs mean(prey Tmax) (nprey ≥ $(min_prey_per_pred))",
        xlabel="Predator Tmax",
        ylabel="Mean prey Tmax",
        outfile=joinpath(outdir, "nodelevel_filt_pred_Tmax_vs_meanprey_Tmax.png")
    )
end

CSV.write(joinpath(outdir, "predator_level_ctmax.csv"), pred_level)

# ---------------------------
# 10) PLOTS: OBSERVED EDGE DISTRIBUTION + NULL DISTRIBUTION
# ---------------------------

hist_plot(
    edges.absdiff_Tmax;
    title="Observed edge distribution: |Tmax(pred) - Tmax(prey)|",
    xlabel="Abs difference (°C)",
    bins=30,
    outfile=joinpath(outdir, "observed_edge_absdiff_Tmax.png")
)

footer =
    "Observed mean = $(round(obs_mean_absdiff, digits=4)) | " *
    "Null mean = $(round(μ0, digits=4)) ± $(round(σ0, digits=4)) | " *
    "z = $(round(z, digits=2)) | p(one-sided, obs ≤ null) = $(round(p_one, digits=4))"

null_hist_plot(
    null_vals, obs_mean_absdiff;
    title="Null distribution: mean |Tmax(pred) - Tmax(prey)| (degree-preserving predator swap)",
    xlabel="Mean |ΔTmax| under null",
    bins=30,
    footer=footer,
    outfile=joinpath(outdir, "null_mean_absdiff_Tmax.png")
)

# ---------------------------
# 11) SUMMARY
# ---------------------------
summary = DataFrame(
    stat = ["mean_absdiff_Tmax"],
    obs = [obs_mean_absdiff],
    null_mean = [μ0],
    null_sd = [σ0],
    z = [z],
    p_one_sided = [p_one],
    edges_matched = [nrow(edges)],
    predators_matched = [length(unique(edges.pred))],
    preys_matched = [length(unique(edges.prey))]
)

CSV.write(joinpath(outdir, "null_test_summary_ctmax.csv"), summary)

@info "Saved outputs to" outdir
println("\n=== DONE ===")
show(summary, allrows=true, allcols=true)
println()
