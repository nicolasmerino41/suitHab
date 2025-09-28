# ==========================================================
# Figure pack (fixed): AM vs BAM — baseline, three grids
# ==========================================================
const Mke = CairoMakie

GRIDS = ["gradient","ridge","fractal"]
datadir = joinpath(@__DIR__, "data")
figdir  = joinpath(@__DIR__, "data", "figs/final")
mkpath(figdir)

# ---------------------------
# IO
# ---------------------------
function load_sweeps(grid::AbstractString)
    res1 = CSV.read(joinpath(datadir, "sweep_C_align_$(grid).csv"), DataFrame)
    res2 = CSV.read(joinpath(datadir, "sweep_R_sigma_$(grid).csv"), DataFrame)
    return res1, res2
end

res1_by_grid = Dict{String,DataFrame}()
res2_by_grid = Dict{String,DataFrame}()
for g in GRIDS
    res1_by_grid[g], res2_by_grid[g] = load_sweeps(g)
end

# ---------------------------
# helpers
# ---------------------------
# Dense grid for heatmaps
function gridZ(df::DataFrame; xname::Symbol, yname::Symbol, zname::Symbol)
    xs = sort(unique(df[!, xname]))
    ys = sort(unique(df[!, yname]))
    lut = Dict( (df[i, xname], df[i, yname]) => df[i, zname] for i in 1:nrow(df) )
    Z = [get(lut, (x,y), NaN) for x in xs, y in ys]
    @assert all(!isnan, Z) "Missing combinations for $(zname)."
    return xs, ys, Z
end

se(v) = std(v) / sqrt(length(v))

# binned P_suff → mean ΔA ± SE
function binned_xy(df::DataFrame; x::Symbol, y::Symbol, nbins::Int=25, xrange::Tuple{Float64,Float64}=nothing)
    xv = collect(skipmissing(df[!, x]))
    yv = collect(skipmissing(df[!, y]))
    if xrange === nothing
        xmin, xmax = minimum(xv), maximum(xv)
    else
        xmin, xmax = xrange
    end
    edges = range(xmin, xmax; length=nbins+1)
    centers = [(edges[i] + edges[i+1]) / 2 for i in 1:nbins]
    bins = [Int[] for _ in 1:nbins]
    for (i,xx) in enumerate(xv)
        b = clamp(searchsortedlast(edges, xx), 1, nbins)
        push!(bins[b], i)
    end
    means, ses = Float64[], Float64[]
    for b in 1:nbins
        idx = bins[b]
        if isempty(idx)
            push!(means, NaN); push!(ses, NaN)
        else
            vals = yv[idx]
            push!(means, mean(vals))
            push!(ses, se(vals))
        end
    end
    (x=centers, y=means, ylo=means .- ses, yhi=means .+ ses, edges=collect(edges))
end

function global_range(res_by_grid::Dict{String,DataFrame}, col::Symbol)
    vals = vcat([collect(skipmissing(df[!, col])) for df in values(res_by_grid)]...)
    (minimum(vals), maximum(vals))
end

# ----------------------------------------------------------
# Figure 01 — ΔArea heatmaps: C × alignment, per grid
# one colorbar **per** heatmap, horizontal below each
# ----------------------------------------------------------
begin
    fig = Figure(; size=(1350, 580))
    for (j,g) in enumerate(GRIDS)
        df = res1_by_grid[g]
        xs, ys, Z = gridZ(df; xname=:align, yname=:C, zname=:ΔAmean)

        # per-panel color range
        zmin, zmax = minimum(Z), maximum(Z)

        ax = Axis(fig[1,j], xlabel="Alignment", ylabel="Connectance C",
                  title="ΔArea (AM − BAM) — align × C ($(g))")
        # Orientation note:
        # Makie heatmap(x,y,Z) expects Z size (length(x), length(y)).
        # Our Z is built as [x,y] → already correct → NO permute here.
        heatmap!(ax, xs, ys, Z; colorrange=(zmin, zmax))
        Colorbar(fig[2,j], limits=(zmin, zmax), vertical=false, label="ΔArea")
    end
    colgap!(fig.layout, 20); rowgap!(fig.layout, 8)
    save(joinpath(figdir, "01_CxAlign_heatmaps_PERPANEL.png"), fig)
    display(fig)
end

# ----------------------------------------------------------
# Figure 02 — ΔArea heatmaps: R95 × σ, per grid
# one colorbar **per** heatmap, horizontal below each
# (fix: use the **unpermuted** Z so high R95 behaves correctly)
# ----------------------------------------------------------
begin
    fig = Figure(; size=(1350, 580))
    for (j,g) in enumerate(GRIDS)
        df = res2_by_grid[g]
        dfc = combine(groupby(df, [:R95, :σ]), :ΔAmean=>mean=>:ΔA)
        xs, ys, Z = gridZ(dfc; xname=:R95, yname=:σ, zname=:ΔA)
        zmin, zmax = minimum(Z), maximum(Z)

        ax = Axis(fig[1,j], xlabel="R95 (diet redundancy)", ylabel="niche breadth σ",
                  title="ΔArea (AM − BAM) — R95 × σ ($(g))")
        # IMPORTANT: use Z (no permute) so gradient along R95 is correct.
        heatmap!(ax, xs, ys, Z; colorrange=(zmin, zmax))
        Colorbar(fig[2,j], limits=(zmin, zmax), vertical=false, label="ΔArea")
    end
    colgap!(fig.layout, 20); rowgap!(fig.layout, 8)
    save(joinpath(figdir, "02_R95xSigma_heatmaps_PERPANEL.png"), fig)
    display(fig)
end

# ----------------------------------------------------------
# Figure 03 — Mediation (binned) for all grids
# ----------------------------------------------------------
begin
    fig = Figure(; size=(1350, 420))

    # use overlapping range across grids to ensure comparable bin edges
    Pmins = Float64[]; Pmaxs = Float64[]
    for g in GRIDS
        df = res1_by_grid[g]
        push!(Pmins, minimum(skipmissing(df.Pmean)))
        push!(Pmaxs, maximum(skipmissing(df.Pmean)))
    end
    Prng = (maximum([0.5, minimum(Pmins)]), minimum([0.95, maximum(Pmaxs)]))

    for (j,g) in enumerate(GRIDS)
        df = res1_by_grid[g]
        b = binned_xy(df; x=:Pmean, y=:ΔAmean, nbins=25, xrange=Prng)
        ax = Axis(fig[1,j], xlabel="P_suff", ylabel="ΔArea",
                  title="Mediation (binned) — $(g)")
        band!(ax, b.x, b.ylo, b.yhi; transparency=true)
        lines!(ax, b.x, b.y)
    end
    save(joinpath(figdir, "03_Mediation_binned_ALL.png"), fig)
    display(fig)
end

# ----------------------------------------------------------
# Figure 04 — Mechanism inserts (π_A, P_suff, π_A(1−P_suff)) vs alignment
#   (robust column-presence check; skip if missing)
# ----------------------------------------------------------
hascol(df, name::AbstractString) = name in names(df)
πA_available = all(hascol(res1_by_grid[g], "πAmean") for g in GRIDS)

if πA_available
    fig = Figure(; size=(1350, 900))
    for (i,g) in enumerate(GRIDS)
        df = sort(res1_by_grid[g], :align)

        ax1 = Axis(fig[i,1], xlabel="alignment", ylabel="π_A (consumers)",
                   title="Abiotic admission vs alignment — $(g)")
        lines!(ax1, df.align, df."πAmean")  # column has Unicode name

        ax2 = Axis(fig[i,2], xlabel="alignment", ylabel="P_suff",
                   title="Prey sufficiency vs alignment — $(g)")
        lines!(ax2, df.align, df.Pmean)

        ax3 = Axis(fig[i,3], xlabel="alignment", ylabel="ΔArea ≈ π_A × (1 − P_suff)",
                   title="Approx ΔArea from mediation — $(g)")
        proxy = df."πAmean" .* (1 .- df.Pmean)
        lines!(ax3, df.align, proxy)
    end
    save(joinpath(figdir, "04_Mechanism_alignment_ALL.png"), fig)
    display(fig)
else
    @info "Skipping Fig 04: column 'πAmean' not present in all res1 CSVs."
end

# ----------------------------------------------------------
# Figure 05 — Lines: ΔArea vs R95 conditioned on σ deciles
#   (fix legend: pass plotted series + labels)
# ----------------------------------------------------------
function sigma_decile_lines(df::DataFrame; ndeciles::Int=10)
    σvals = collect(skipmissing(df.σ))
    edges = quantile(σvals, range(0,1; length=ndeciles+1))
    cutσ(s) = clamp(searchsortedlast(edges, s), 1, ndeciles)
    df = transform(df, :σ => ByRow(cutσ) => :σdec)
    tmp = combine(groupby(df, [:σdec, :R95]), :ΔAmean=>mean=>:ΔA)

    xs = sort(unique(tmp.R95))
    linesY = Vector{Vector{Float64}}(undef, ndeciles)
    labels = Vector{String}(undef, ndeciles)
    for d in 1:ndeciles
        sub = tmp[tmp.σdec .== d, :]
        lut = Dict(sub.R95 .=> sub.ΔA)
        linesY[d] = [get(lut, x, NaN) for x in xs]
        labels[d] = "σ ∈ [" * @sprintf("%.3f", edges[d]) * "," * @sprintf("%.3f", edges[d+1]) * "]"
    end
    return xs, linesY, labels
end

begin
    fig = Figure(; size=(1350, 420))
    for (j,g) in enumerate(GRIDS)
        df = res2_by_grid[g]
        xs, linesY, labels = sigma_decile_lines(df; ndeciles=10)
        ax = Axis(fig[1,j], xlabel="R95 (diet redundancy)", ylabel="ΔArea (AM − BAM)",
                  title="ΔArea vs R95 | σ deciles — $(g)")
        series = Any[]
        for d in eachindex(linesY)
            push!(series, lines!(ax, xs, linesY[d]))
        end
        axislegend(ax, series, labels; position=:rb, nbanks=2, framevisible=false, labelsize=8)
    end
    save(joinpath(figdir, "05_R95_lines_by_sigmaDeciles_ALL.png"), fig)
    display(fig)
end

# ==========================================================
# SUPPLEMENTARY / SI FIGURES
# ==========================================================

# ---------- SI-1: Variability maps + fraction-above + ECDFs ----------
# We don’t have replicate-level ΔA; use proxy variability = half-band width.
var_proxy(df) = (df.ΔAhi .- df.ΔAlo) ./ 2

begin
    # (a) variability heatmaps for C×align and R95×σ, per grid
    figA = Figure(; size=(1350, 900))
    for (row,(grid, df)) in enumerate(res1_by_grid)
        xs, ys, Zμ = gridZ(df; xname=:C, yname=:align, zname=:ΔAmean)
        _,  _, Zv = gridZ(transform(df, :ΔAmean=>(_->nothing)=>:tmp, [:ΔAlo,:ΔAhi]=>((lo,hi)->(hi.-lo)./2)=>:var) ;
                          xname=:C, yname=:align, zname=:var)
        zmin, zmax = minimum(Zv), maximum(Zv)
        ax = Axis(figA[row,1], xlabel="C", ylabel="align", title="Var proxy (C×align) — $(grid)")
        heatmap!(ax, xs, ys, Zv; colorrange=(zmin, zmax))
        Colorbar(figA[row,2], limits=(zmin, zmax), vertical=false, label="proxy SD")
    end
    for (row,(grid, df)) in enumerate(res2_by_grid)
        dfc = transform(combine(groupby(df, [:R95,:σ]), :ΔAmean=>mean=>:ΔA, :ΔAlo=>mean=>:ΔAlo, :ΔAhi=>mean=>:ΔAhi),
                        [:ΔAlo,:ΔAhi]=>((lo,hi)->(hi.-lo)./2)=>:var)
        xs, ys, Zv = gridZ(dfc; xname=:R95, yname=:σ, zname=:var)
        zmin, zmax = minimum(Zv), maximum(Zv)
        ax = Axis(figA[row,3], xlabel="R95", ylabel="σ", title="Var proxy (R95×σ) — $(grid)")
        heatmap!(ax, xs, ys, Zv; colorrange=(zmin, zmax))
        Colorbar(figA[row,4], limits=(zmin, zmax), vertical=false, label="proxy SD")
    end
    save(joinpath(figdir, "SI01_Variability_heatmaps.png"), figA)
    display(figA)

    # (b) fraction-above-threshold maps (simple: mean>τ)
    τ = 0.30
    figB = Figure(; size=(1350, 580))
    for (j,g) in enumerate(GRIDS)
        df = res1_by_grid[g]
        df2 = transform(df, :ΔAmean => ByRow(>(τ)) => :frac)
        xs, ys, Z = gridZ(df2; xname=:C, yname=:align, zname=:frac)
        ax = Axis(figB[1,j], xlabel="C", ylabel="align", title="Frac(ΔArea>$(τ)) — C×align ($(g))")
        heatmap!(ax, xs, ys, Z)
        Colorbar(figB[2,j], vertical=false, label="fraction>τ")
    end
    save(joinpath(figdir, "SI01_FractionAbove_CxAlign.png"), figB)
    display(figB)

    # (c) ECDFs of ΔAmean across three regimes of (C,align): low, mid, high C
    #    (purely descriptive, based on grid of res1)
    function ecdf_data(v::AbstractVector)
        sv = sort(v); n = length(sv)
        (x = sv, y = collect(1:n) ./ n)
    end
    figC = Figure(; size=(1350, 420))
    for (j,g) in enumerate(GRIDS)
        df = res1_by_grid[g]
        low  = df[df.C .<= 0.06, :]
        mid  = df[(df.C .> 0.06) .& (df.C .<= 0.12), :]
        high = df[df.C .> 0.12, :]

        ax = Axis(figC[1,j], xlabel="ΔArea", ylabel="ECDF",
                  title="ECDF(ΔArea mean) by C regime — $(g)")
        series, labels = Any[], String[]
        for (sub, lab) in ((low,"low C"), (mid,"mid C"), (high,"high C"))
            e = ecdf_data(sub.ΔAmean)
            push!(series, lines!(ax, e.x, e.y))
            push!(labels, lab)
        end
        axislegend(ax, series, labels; position=:lt, framevisible=false, labelsize=9)
    end
    save(joinpath(figdir, "SI01_ECDF_regimes.png"), figC)
    display(figC)
end

# ---------- SI-2: Scaling diagnostics (optional: skip if CSVs not found) ----------
function maybe_read(path)
    isfile(path) ? CSV.read(path, DataFrame) : nothing
end

begin
    df_grid = maybe_read(joinpath(datadir, "scaling_grid_gradient.csv"))     # expects cols: cells, metric, value, grid
    df_S    = maybe_read(joinpath(datadir, "scaling_species_gradient.csv"))  # expects cols: S, metric, value, grid
    if !(df_grid === nothing)
        fig = Figure(; size=(1100, 380))
        for (j,g) in enumerate(GRIDS)
            sub = df_grid[df_grid.grid .== g, :]
            ax = Axis(fig[1,j], xscale=log10, xlabel="number of cells", ylabel="metric",
                      title="Finite-size convergence — $(g)")
            # plot ΔArea, ΔGini, P_suff if present
            for m in unique(sub.metric)
                s = sub[sub.metric .== m, :]
                lines!(ax, s.cells, s.value)
            end
        end
        save(joinpath(figdir, "SI02_Scaling_grid.png"), fig)
        display(fig)
    else
        @info "SI-2 grid scaling: file not found → skipped."
    end

    if !(df_S === nothing)
        fig = Figure(; size=(500, 380))
        ax = Axis(fig[1,1], xlabel="S", ylabel="metric value", title="Species scaling (all grids pooled)")
        for m in unique(df_S.metric)
            s = df_S[df_S.metric .== m, :]
            lines!(ax, s.S, s.value)
        end
        save(joinpath(figdir, "SI02_Scaling_S.png"), fig)
        display(fig)
    else
        @info "SI-2 species scaling: file not found → skipped."
    end
end

# ---------- SI-3: Mechanism inserts (again, SI wrapper for all grids) ----------
if πA_available
    fig = Figure(; size=(1350, 900))
    for (i,g) in enumerate(GRIDS)
        df = sort(res1_by_grid[g], :align)
        ax1 = Axis(fig[i,1], xlabel="alignment", ylabel="π_A", title="π_A vs alignment — $(g)")
        lines!(ax1, df.align, df."πAmean")
        ax2 = Axis(fig[i,2], xlabel="alignment", ylabel="P_suff", title="P_suff vs alignment — $(g)")
        lines!(ax2, df.align, df.Pmean)
        ax3 = Axis(fig[i,3], xlabel="alignment", ylabel="π_A(1−P_suff)", title="Proxy vs alignment — $(g)")
        lines!(ax3, df.align, df."πAmean" .* (1 .- df.Pmean))
    end
    save(joinpath(figdir, "SI03_Mechanisms_alignment_ALL.png"), fig)
    display(fig)
end

@info "All figures written to $(figdir)"
