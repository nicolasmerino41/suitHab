#!/usr/bin/env julia
# GENUS-LEVEL ANALOGUE OF SPECIES-LEVEL SCRIPT
#
# Differences vs species script:
# - Traits aggregated by Genus
# - Metaweb matched on predator/prey genus
# - All downstream analyses identical

using DataFrames, Statistics, Printf
using CairoMakie

# =========================
# 0) Output directory
# =========================
OUTDIR = joinpath(pwd(), "empirical_globtherm_metaweb_outputs_genus")
isdir(OUTDIR) || mkpath(OUTDIR)

# =========================
# 1) Utilities
# =========================

function normalize_name(s::AbstractString)
    str = String(s)
    str = replace(str, '\u00A0' => ' ')
    str = replace(str, r"\p{Z}+" => " ")
    buf = IOBuffer()
    for c in str
        if c == ' '
            write(buf, c)
        else
            try
                write(buf, lowercase(c))
            catch
            end
        end
    end
    return strip(String(take!(buf)))
end

function pick_first_nonmissing(row, cols::Vector{Symbol})
    for c in cols
        if c ∈ propertynames(row)
            v = row[c]
            if v !== missing
                if v isa AbstractString
                    isempty(strip(v)) && continue
                    return v
                else
                    return v
                end
            end
        end
    end
    return missing
end

function cor_pearson(x::Vector{Float64}, y::Vector{Float64})
    (length(x) < 3 || std(x) == 0 || std(y) == 0) && return (missing, length(x))
    return (cor(x, y), length(x))
end

function rankdata(v::Vector{Float64})
    p = sortperm(v); r = similar(v)
    i = 1
    while i ≤ length(v)
        j = i
        while j < length(v) && v[p[j+1]] == v[p[i]]; j += 1; end
        for k in i:j; r[p[k]] = 0.5*(i+j); end
        i = j + 1
    end
    r
end

cor_spearman(x,y) = cor_pearson(rankdata(x), rankdata(y))

directed_overlap(a_lo,a_hi,b_lo,b_hi) =
    max(0.0, min(a_hi,b_hi)-max(a_lo,b_lo)) / max(a_hi-a_lo, eps())

# =========================
# 2) Build GlobalTherm GENUS traits
# =========================

function build_gt_traits_genus(gt::DataFrame)
    g = copy(gt)

    g.genus_norm = Union{Missing,String}[
        (x === missing ? missing : normalize_name(x)) for x in g.Genus
    ]
    function tofloat(x)
        x === missing && return missing

        if x isa Real
            return Float64(x)
        elseif x isa AbstractString
            s = strip(String(x))
            isempty(s) && return missing
            v = tryparse(Float64, s)
            return v === nothing ? missing : v
        else
            return missing
        end
    end

    Tmax_cols = [:Tmax, :Tmax_2]
    Tmin_cols = [:tmin, :tmin_2]

    g.Tmax_best = [
        tofloat(pick_first_nonmissing(g[i,:], Tmax_cols))
        for i in 1:nrow(g)
    ]

    g.Tmin_best = [
        tofloat(pick_first_nonmissing(g[i,:], Tmin_cols))
        for i in 1:nrow(g)
    ]

    gb = groupby(g, :genus_norm; skipmissing=true)

    gt_gen = combine(gb,
        :Genus     => first         => :Genus,
        :Tmax_best => mean ∘ skipmissing => :Tmax,
        :Tmin_best => mean ∘ skipmissing => :Tmin,
        :Class     => first         => :Class,
        :Order     => first         => :Order,
        :Family    => first         => :Family
    )

    gt_gen.Tmid    = 0.5 .* (gt_gen.Tmax .+ gt_gen.Tmin)
    gt_gen.breadth = gt_gen.Tmax .- gt_gen.Tmin

    return gt_gen
end

gt_traits = build_gt_traits_genus(gt)

@info "GlobalTherm genus trait coverage" (
    genera = nrow(gt_traits),
    has_Tmid = count(!ismissing, gt_traits.Tmid)
)

# =========================
# 3) Match METAWEB by GENUS
# =========================

gt_lookup = Dict{String,NamedTuple}()

for r in eachrow(gt_traits)
    gt_lookup[normalize_name(r.Genus)] = (
        Tmax=r.Tmax, Tmin=r.Tmin, Tmid=r.Tmid, breadth=r.breadth
    )
end

function genus_of_name(s::AbstractString)
    first(split(s, " "))
end

function attach_traits_to_web_genus(web::DataFrame, gt_lookup)
    w = copy(web)
    w.pred_gen = normalize_name.(genus_of_name.(w.predator))
    w.prey_gen = normalize_name.(genus_of_name.(w.prey))

    for col in [:pred_Tmax,:pred_Tmin,:pred_Tmid,:pred_breadth,
                :prey_Tmax,:prey_Tmin,:prey_Tmid,:prey_breadth]
        w[!, col] = Vector{Union{Missing,Float64}}(missing, nrow(w))
    end

    for i in 1:nrow(w)
        haskey(gt_lookup, w.pred_gen[i]) && begin
            t = gt_lookup[w.pred_gen[i]]
            w.pred_Tmax[i]=t.Tmax; w.pred_Tmin[i]=t.Tmin
            w.pred_Tmid[i]=t.Tmid; w.pred_breadth[i]=t.breadth
        end
        haskey(gt_lookup, w.prey_gen[i]) && begin
            t = gt_lookup[w.prey_gen[i]]
            w.prey_Tmax[i]=t.Tmax; w.prey_Tmin[i]=t.Tmin
            w.prey_Tmid[i]=t.Tmid; w.prey_breadth[i]=t.breadth
        end
    end
    w
end

webT = attach_traits_to_web_genus(web, gt_lookup)

# =========================
# 4) ANALYSES (UNCHANGED)
# =========================

traits = [:Tmin, :Tmax, :Tmid, :breadth]

function edge_level_corr(webT, trait)
    px = Symbol("pred_", trait); py = Symbol("prey_", trait)
    x=Float64[]; y=Float64[]
    for r in eachrow(webT)
        r[px]!==missing && r[py]!==missing && (push!(x,r[px]); push!(y,r[py]))
    end
    rp,n = cor_pearson(x,y); rs,_ = cor_spearman(x,y)
    (pearson=rp, spearman=rs, n=n, x=x, y=y)
end

edge_results = DataFrame(trait=String[], pearson=Any[], spearman=Any[], n=Int[])
for tr in traits
    e = edge_level_corr(webT,tr)
    push!(edge_results,(string(tr),e.pearson,e.spearman,e.n))
end

CSV.write(joinpath(OUTDIR,"edge_level_correlations_genus.csv"), edge_results)

# =========================
# 5) Plots (UNCHANGED)
# =========================

function scatter_with_fit(x,y; title="", xlabel="", ylabel="", outpath=nothing)
    f=Figure(size=(800,650)); ax=Axis(f[1,1],title=title,xlabel=xlabel,ylabel=ylabel)
    scatter!(ax,x,y)
    if length(x)>=3 && std(x)>0
        β=hcat(ones(length(x)),x)\y
        xs=range(minimum(x),maximum(x),length=100)
        lines!(ax,xs,β[1].+β[2].*xs)
    end
    outpath!==nothing && save(outpath,f)
    f
end

for tr in traits
    x,y = edge_level_corr(webT,tr).x, edge_level_corr(webT,tr).y
    f = scatter_with_fit(x,y;
        title="Genus-level: predator vs prey $tr",
        xlabel="Predator $tr",
        ylabel="Prey $tr",
        outpath=joinpath(OUTDIR,"scatter_edge_genus_$(tr).png")
    )
    display(f)
end

println("\n=== GENUS-LEVEL DONE ===")
println("Outputs in: $OUTDIR")
