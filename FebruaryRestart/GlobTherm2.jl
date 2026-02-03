#!/usr/bin/env julia
# Run after your loader script that defines:
#   gt :: DataFrame   (GlobalTherm)
#   web :: DataFrame  (predator, prey)
#
# Output: CSV summaries + plots in OUTDIR.

using DataFrames, Statistics, Printf
using CairoMakie

# =========================
# 0) Output directory
# =========================
OUTDIR = joinpath(pwd(), "empirical_globtherm_metaweb_outputs")
isdir(OUTDIR) || mkpath(OUTDIR)

# =========================
# 1) Utilities
# =========================
function normalize_name(s::AbstractString)
    # 1) Force materialization (kills SubString issues)
    str = String(s)

    # 2) Replace NBSP and all Unicode space separators
    str = replace(str, '\u00A0' => ' ')
    str = replace(str, r"\p{Z}+" => " ")

    # 3) Safe lowercase: operate char-by-char
    buf = IOBuffer()
    for c in str
        if c == ' '
            write(buf, c)
        else
            # lowercase only if valid
            try
                write(buf, lowercase(c))
            catch
                # drop invalid characters silently
                # (or replace with space if you prefer)
            end
        end
    end

    return strip(String(take!(buf)))
end

# Pick best value from a list of columns in order
function pick_first_nonmissing(row, cols::Vector{Symbol})
    for c in cols
        if c ∈ propertynames(row)
            v = row[c]
            if v !== missing && !(v isa AbstractString && isempty(strip(v)))
                return v
            end
        end
    end
    return missing
end

# Safe Pearson correlation
function cor_pearson(x::Vector{Float64}, y::Vector{Float64})
    (length(x) < 3) && return (missing, length(x))
    (std(x) == 0 || std(y) == 0) && return (missing, length(x))
    return (cor(x, y), length(x))
end

# Rank correlation (Spearman) without external packages
function rankdata(v::Vector{Float64})
    p = sortperm(v)
    r = similar(v)
    # average ranks for ties
    i = 1
    while i <= length(v)
        j = i
        while j < length(v) && v[p[j+1]] == v[p[i]]
            j += 1
        end
        avg = 0.5 * (i + j)
        for k in i:j
            r[p[k]] = avg
        end
        i = j + 1
    end
    return r
end

function cor_spearman(x::Vector{Float64}, y::Vector{Float64})
    (length(x) < 3) && return (missing, length(x))
    rx = rankdata(x)
    ry = rankdata(y)
    (std(rx) == 0 || std(ry) == 0) && return (missing, length(x))
    return (cor(rx, ry), length(x))
end

# Interval overlap: |A∩B| / |A| (directed, prey-support style)
function directed_overlap(a_lo::Float64, a_hi::Float64, b_lo::Float64, b_hi::Float64)
    lenA = max(0.0, a_hi - a_lo)
    lenA == 0.0 && return 0.0
    lo = max(a_lo, b_lo)
    hi = min(a_hi, b_hi)
    inter = max(0.0, hi - lo)
    return inter / lenA
end

# Safely extract a single string from messy taxonomy fields
function scalar_string(x)
    x === missing && return missing
    x isa AbstractString && return x
    x isa AbstractVector && !isempty(x) && return x[1]
    return missing
end

function sanitize_string_column!(df::DataFrame, col::Symbol)
    col ∈ propertynames(df) || return

    v = df[!, col]
    eltype(v) === String || return

    n = length(v)
    newv = Vector{Union{Missing,String}}(undef, n)

    for i in 1:n
        if !isassigned(v, i)
            newv[i] = missing
        else
            x = v[i]
            if x isa String
                newv[i] = x
            elseif x isa AbstractVector{<:AbstractString} && !isempty(x)
                newv[i] = x[1]          # 🔑 scalarize here
            else
                newv[i] = missing
            end
        end
    end

    # replace column storage directly (no coercion)
    df.columns[df.colindex[col]] = newv
end

for c in (:Genus, :Species, :Class, :Order, :Family)
    sanitize_string_column!(gt, c)
end

# =========================
# 2) Build GlobalTherm species traits table
# =========================
# GlobalTherm columns you showed include:
#   "Genus","Species","Tmax","Tmax_2","tmin","tmin_2", plus metadata.
using DataFrames, Statistics

# Your existing normalize_name + pick_first_nonmissing can stay as-is, but keep ONE copy only.
# I'll assume you already defined:
#   normalize_name(::AbstractString)
#   pick_first_nonmissing(row, cols)

@inline function safe_get(v, i)
    (1 <= i <= length(v) && isassigned(v, i)) ? v[i] : missing
end

function tofloat(x)
    x === missing && return missing
    x isa Real && return Float64(x)
    if x isa AbstractString
        t = tryparse(Float64, strip(String(x)))
        return t === nothing ? missing : t
    end
    return missing
end

function first_scalar(v)
    for x in v
        x === missing && continue
        if x isa AbstractVector{<:AbstractString} && !isempty(x)
            return x[1]
        elseif x isa AbstractString
            return x
        end
    end
    return missing
end

function build_gt_traits(gt::DataFrame)
    g = copy(gt)
    n = nrow(g)

    genus   = g[!, :Genus]
    species = g[!, :Species]

    binomial      = fill!(Vector{Union{Missing,String}}(undef, n), missing)
    binomial_norm = fill!(Vector{Union{Missing,String}}(undef, n), missing)

    for i in 1:n
        gen = safe_get(genus, i)
        sp  = safe_get(species, i)

        if gen === missing || sp === missing
            continue
        end

        name = string(gen, " ", sp)
        binomial[i]      = name
        binomial_norm[i] = normalize_name(name)
    end

    g.binomial = binomial
    g.binomial_norm = binomial_norm

    Tmax_cols = [:Tmax, :Tmax_2]
    Tmin_cols = [:tmin, :tmin_2]

    Tmax_best = fill!(Vector{Union{Missing,Float64}}(undef, n), missing)
    Tmin_best = fill!(Vector{Union{Missing,Float64}}(undef, n), missing)

    for i in 1:n
        row = g[i, :]
        Tmax_best[i] = tofloat(pick_first_nonmissing(row, Tmax_cols))
        Tmin_best[i] = tofloat(pick_first_nonmissing(row, Tmin_cols))
    end

    g.Tmax_best = Tmax_best
    g.Tmin_best = Tmin_best

    gb = groupby(g, :binomial_norm; skipmissing=true)

    specs = Any[
        :binomial   => first_scalar => :binomial,
        :Tmax_best  => (v -> isempty(skipmissing(v)) ? missing : mean(skipmissing(v))) => :Tmax,
        :Tmin_best  => (v -> isempty(skipmissing(v)) ? missing : mean(skipmissing(v))) => :Tmin,
        :Class      => first_scalar => :Class,
        :Order      => first_scalar => :Order,
        :Family     => first_scalar => :Family
    ]

    if :N ∈ propertynames(g)
        push!(specs, :N => (v -> isempty(skipmissing(v)) ? missing : maximum(skipmissing(v))) => :Nmax)
    end

    gt_sp = combine(gb, specs...)

    m = nrow(gt_sp)
    gt_sp.Tmid    = fill!(Vector{Union{Missing,Float64}}(undef, m), missing)
    gt_sp.breadth = fill!(Vector{Union{Missing,Float64}}(undef, m), missing)

    for i in 1:m
        tmax = gt_sp.Tmax[i]
        tmin = gt_sp.Tmin[i]
        if tmax === missing || tmin === missing
            continue
        end
        gt_sp.Tmid[i] = 0.5 * (tmax + tmin)
        gt_sp.breadth[i] = tmax - tmin
    end

    return gt_sp
end

gt_traits = build_gt_traits(gt)

@info "GlobalTherm trait coverage" (
    species_rows = nrow(gt_traits),
    has_Tmax = count(!ismissing, gt_traits.Tmax),
    has_Tmin = count(!ismissing, gt_traits.Tmin),
    has_Tmid = count(!ismissing, gt_traits.Tmid),
    has_breadth = count(!ismissing, gt_traits.breadth)
)

# =========================
# 3) Match EU metaweb names to GlobalTherm
# =========================
# Best-effort exact normalized match on binomial.
# If EU names contain subspecies/author strings, you’ll need a smarter parser later.

gt_lookup = Dict{String,NamedTuple}()
for r in eachrow(gt_traits)
    gt_lookup[r.binomial_norm] = (;
        binomial = r.binomial,
        Tmax = r.Tmax,
        Tmin = r.Tmin,
        Tmid = r.Tmid,
        breadth = r.breadth
    )
end

# Attach traits to each interaction (edge-level table)
function attach_traits_to_web(web::DataFrame, gt_lookup::Dict{String,NamedTuple})
    w = copy(web)
    w.pred_norm = normalize_name.(w.predator)
    w.prey_norm = normalize_name.(w.prey)

    # initialize columns
    for col in [:pred_Tmax,:pred_Tmin,:pred_Tmid,:pred_breadth,
                :prey_Tmax,:prey_Tmin,:prey_Tmid,:prey_breadth]
        w[!, col] = Vector{Union{Missing,Float64}}(missing, nrow(w))
    end

    for i in 1:nrow(w)
        p = w.pred_norm[i]
        q = w.prey_norm[i]
        if haskey(gt_lookup, p)
            t = gt_lookup[p]
            w.pred_Tmax[i] = t.Tmax
            w.pred_Tmin[i] = t.Tmin
            w.pred_Tmid[i] = t.Tmid
            w.pred_breadth[i] = t.breadth
        end
        if haskey(gt_lookup, q)
            t = gt_lookup[q]
            w.prey_Tmax[i] = t.Tmax
            w.prey_Tmin[i] = t.Tmin
            w.prey_Tmid[i] = t.Tmid
            w.prey_breadth[i] = t.breadth
        end
    end
    return w
end

webT = attach_traits_to_web(web, gt_lookup)

edge_cov = (
    edges = nrow(webT),
    edges_with_Tmid_both = count(i -> (webT.pred_Tmid[i] !== missing && webT.prey_Tmid[i] !== missing), 1:nrow(webT)),
    edges_with_Tmax_both = count(i -> (webT.pred_Tmax[i] !== missing && webT.prey_Tmax[i] !== missing), 1:nrow(webT)),
    edges_with_Tmin_both = count(i -> (webT.pred_Tmin[i] !== missing && webT.prey_Tmin[i] !== missing), 1:nrow(webT)),
)
@info "Edge-level trait coverage" edge_cov

# =========================
# 4) Correlation analyses
# =========================

# 4.1 Edge-level correlation (pred trait vs prey trait over interactions)
function edge_level_corr(webT::DataFrame, trait::Symbol)
    pred_col = Symbol("pred_", trait)
    prey_col = Symbol("prey_", trait)
    x = Float64[]
    y = Float64[]
    for r in eachrow(webT)
        px = r[pred_col]
        py = r[prey_col]
        if px !== missing && py !== missing
            push!(x, Float64(px))
            push!(y, Float64(py))
        end
    end
    rp, n = cor_pearson(x, y)
    rs, _ = cor_spearman(x, y)
    return (pearson=rp, spearman=rs, n=n, x=x, y=y)
end

# 4.2 Predator-vs-mean(prey) correlation (mechanistic proxy)
function predator_vs_preymean_corr(webT::DataFrame, trait::Symbol; min_prey_with_data::Int=3)
    pred_col = Symbol("pred_", trait)
    prey_col = Symbol("prey_", trait)

    predators = unique(webT.predator)
    X = Float64[]  # predator trait
    Y = Float64[]  # mean prey trait

    # Also store per-predator diagnostics
    rows = NamedTuple[]

    for p in predators
        sub = webT[webT.predator .== p, :]
        pvals = sub[!, pred_col]
        # predator trait (mean across edges, should be same, but safe)
        pred_trait = (all(ismissing, pvals) ? missing : mean(skipmissing(pvals)))

        prey_vals = sub[!, prey_col]
        prey_vals_ok = collect(skipmissing(prey_vals))

        if pred_trait === missing || length(prey_vals_ok) < min_prey_with_data
            continue
        end

        prey_mean = mean(prey_vals_ok)
        push!(X, Float64(pred_trait))
        push!(Y, Float64(prey_mean))

        push!(rows, (; predator=p, pred_trait=Float64(pred_trait),
                      prey_mean=Float64(prey_mean), n_prey=length(prey_vals_ok)))
    end

    rp, n = cor_pearson(X, Y)
    rs, _ = cor_spearman(X, Y)
    return (pearson=rp, spearman=rs, n=n, X=X, Y=Y, table=DataFrame(rows))
end

# 4.3 Overlap-based “support alignment”: predator interval overlap with mean prey interval
# This is NOT the same as correlation; it's a mechanistic alignment score in [0,1].
function predator_prey_overlap(webT::DataFrame; min_prey_with_data::Int=3)
    predators = unique(webT.predator)
    rows = NamedTuple[]

    for p in predators
        sub = webT[webT.predator .== p, :]

        pTmin = (all(ismissing, sub.pred_Tmin) ? missing : mean(skipmissing(sub.pred_Tmin)))
        pTmax = (all(ismissing, sub.pred_Tmax) ? missing : mean(skipmissing(sub.pred_Tmax)))

        preyTmin_ok = collect(skipmissing(sub.prey_Tmin))
        preyTmax_ok = collect(skipmissing(sub.prey_Tmax))

        if pTmin === missing || pTmax === missing
            continue
        end
        if length(preyTmin_ok) < min_prey_with_data || length(preyTmax_ok) < min_prey_with_data
            continue
        end

        preyTmin_mean = mean(preyTmin_ok)
        preyTmax_mean = mean(preyTmax_ok)

        align = directed_overlap(Float64(pTmin), Float64(pTmax), Float64(preyTmin_mean), Float64(preyTmax_mean))
        push!(rows, (; predator=p, align=align,
                      pred_Tmin=Float64(pTmin), pred_Tmax=Float64(pTmax),
                      prey_Tmin_mean=preyTmin_mean, prey_Tmax_mean=preyTmax_mean,
                      n_prey=min(length(preyTmin_ok), length(preyTmax_ok))))
    end
    return DataFrame(rows)
end

# Run correlations for key thermal traits
traits = [:Tmin, :Tmax, :Tmid, :breadth]

edge_results = DataFrame(trait=String[], pearson=Any[], spearman=Any[], n=Int[])
predmean_results = DataFrame(trait=String[], pearson=Any[], spearman=Any[], n=Int[])

edge_xy = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()
predmean_xy = Dict{Symbol,Tuple{Vector{Float64},Vector{Float64}}}()
predmean_tables = Dict{Symbol,DataFrame}()

for tr in traits
    e = edge_level_corr(webT, tr)
    push!(edge_results, (string(tr), e.pearson, e.spearman, e.n))
    edge_xy[tr] = (e.x, e.y)

    pm = predator_vs_preymean_corr(webT, tr; min_prey_with_data=3)
    push!(predmean_results, (string(tr), pm.pearson, pm.spearman, pm.n))
    predmean_xy[tr] = (pm.X, pm.Y)
    predmean_tables[tr] = pm.table
end

overlap_tbl = predator_prey_overlap(webT; min_prey_with_data=3)

# Save tables
CSV.write(joinpath(OUTDIR, "edge_level_correlations.csv"), edge_results)
CSV.write(joinpath(OUTDIR, "predator_vs_preymean_correlations.csv"), predmean_results)
CSV.write(joinpath(OUTDIR, "predator_vs_preymean_table_Tmid.csv"), predmean_tables[:Tmid])
CSV.write(joinpath(OUTDIR, "predator_prey_alignment_overlap.csv"), overlap_tbl)

@info "Saved summaries to" OUTDIR

# =========================
# 5) Plots
# =========================

function scatter_with_fit(x::Vector{Float64}, y::Vector{Float64}; title="", xlabel="", ylabel="", outpath=nothing)
    f = Figure(size=(800, 650))
    ax = Axis(f[1,1], title=title, xlabel=xlabel, ylabel=ylabel)

    scatter!(ax, x, y)

    # simple least squares fit line
    if length(x) >= 3 && std(x) > 0
        X = hcat(ones(length(x)), x)
        β = X \ y
        xs = range(minimum(x), maximum(x), length=100)
        ys = β[1] .+ β[2] .* xs
        lines!(ax, xs, ys)
    end

    if outpath !== nothing
        save(outpath, f)
    end
    return f
end

# Edge-level scatter
for tr in traits
    x, y = edge_xy[tr]
    rp, n = cor_pearson(x, y)
    ttl = @sprintf("Edge-level: predator %s vs prey %s (Pearson=%s, n=%d)",
                   string(tr), string(tr),
                   rp === missing ? "NA" : @sprintf("%.2f", rp), n)
    f = scatter_with_fit(x, y;
        title=ttl,
        xlabel="Predator $(tr)",
        ylabel="Prey $(tr)",
        outpath=joinpath(OUTDIR, "scatter_edge_level_$(tr).png")
    )
    display(f)
end

# Predator vs mean(prey) scatter (mechanistic proxy)
for tr in traits
    x, y = predmean_xy[tr]
    rp, n = cor_pearson(x, y)
    ttl = @sprintf("Predator vs mean(prey): %s (Pearson=%s, n=%d)",
                   string(tr),
                   rp === missing ? "NA" : @sprintf("%.2f", rp), n)
    f = scatter_with_fit(x, y;
        title=ttl,
        xlabel="Predator $(tr)",
        ylabel="Mean prey $(tr)",
        outpath=joinpath(OUTDIR, "scatter_pred_vs_preymean_$(tr).png")
    )
    display(f)
end

# Overlap alignment distribution
if nrow(overlap_tbl) > 0
    f = Figure(size=(800, 650))
    ax = Axis(f[1,1], title="Predator–prey thermal alignment (directed overlap)", xlabel="Overlap (0–1)", ylabel="Count")
    hist!(ax, Float64.(overlap_tbl.align), bins=30)
    save(joinpath(OUTDIR, "hist_pred_prey_alignment_overlap.png"), f)
    display(f)
end

println("\n=== DONE ===")
println("Outputs in: $(OUTDIR)")
println("Key files:")
println(" - edge_level_correlations.csv")
println(" - predator_vs_preymean_correlations.csv")
println(" - predator_prey_alignment_overlap.csv")
