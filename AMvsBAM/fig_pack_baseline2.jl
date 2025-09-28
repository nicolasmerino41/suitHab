# ======================================================================
# make_all_figures.jl
# Reads res1/res2 for three grids and produces all paper figures.
# Figures are saved into data/figs/ as 01_..., 02_..., etc.
# - Heatmaps: separate horizontal colorbar per grid (not shared).
# - R95–σ: correct (no unintended permute), proper axes orientation.
# - Mediation: consumers-only, P_suff computed within A-kept consumer cells,
#              binned by P_suff with mean ± SE (reduces 2-cloud noise).
# - Mechanism inserts: πA (consumers), P_suff, and πA*(1-P_suff) vs alignment,
#                      smoothed via alignment bins (mean ± SE).
# - Variability maps: use IQR→SD proxy from ΔAhi/ΔAlo; fixed ranges/labels.
# ======================================================================

# ==========================================================
# Figure pack (fixed): AM vs BAM — baseline, three grids
# ==========================================================
const Mke = CairoMakie

GRIDS = ["gradient","ridge","fractal"]
datadir = joinpath(@__DIR__, "data")
figdir  = joinpath(@__DIR__, "data", "figs/final2")
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
# ---------- utilities ----------
binmean(x::AbstractVector, y::AbstractVector; nbins::Int=30, lohi=nothing) = begin
    @assert length(x) == length(y)
    lo, hi = isnothing(lohi) ? (minimum(x), maximum(x)) : lohi
    edges = range(lo, hi; length=nbins+1)
    mids  = 0.5 .* (edges[1:end-1] .+ edges[2:end])
    μ     = similar(mids); se = similar(mids)
    for i in eachindex(mids)
        m = (x .>= edges[i]) .& (x .< edges[i+1])
        if any(m)
            vals = y[m]
            μ[i]  = mean(vals)
            se[i] = std(vals) / sqrt(length(vals))
        else
            μ[i]  = NaN; se[i] = NaN
        end
    end
    (x=mids, y=μ, se=se)
end

# quick SE helper for columns already aggregated
se(v) = std(v) / sqrt(length(v))

# turn ΔA IQR band into SD proxy (ΔAhi-ΔAlo is ~IQR if qband used at .25/.75)
iqr_to_sd(iqr) = iqr / 1.349

# ----- Figure 01: C×align heatmaps side-by-side (one colorbar per grid) -----
function fig01_heatmaps_C_align(res1_by_grid)
    fig = Figure(size=(1500, 520))
    for (j,g) in enumerate(GRIDS)
        df = res1_by_grid[g]
        xs = sort(unique(df.C)); ys = sort(unique(df.align))
        # build matrix Z[ix,iy] with x across columns, y down rows
        lut = Dict((df.C[i], df.align[i]) => df.ΔAmean[i] for i in 1:nrow(df))
        Z = [lut[(x,y)] for y in ys, x in xs]  # rows=y, cols=x (no transpose)
        ax = Axis(fig[1,j], title="ΔArea (AM−BAM) — C×align — $(g)",
                  xlabel="connectance C", ylabel=(j==1 ? "alignment" : ""))
        hm = heatmap!(ax, xs, ys, Z; interpolate=false)
        Colorbar(fig[2,j], hm, vertical=false, label="ΔArea")
    end
    save(joinpath(figdir, "01_HM_DeltaArea_C_align_all_grids.png"), fig)
    display(fig)
end

# ----- Figure 02 (FIXED): R95×σ heatmaps side-by-side (one colorbar per grid) -----
function fig02_heatmaps_R_sigma(res2_by_grid)
    fig = Figure(size=(1500, 520))
    for (j,g) in enumerate(GRIDS)
        df = res2_by_grid[g]

        # collapse replicates (if any) so we have one value per (R95, σ)
        dfc = combine(groupby(df, [:R95, :σ]), :ΔAmean => mean => :ΔA)

        # axes
        xs = sort(unique(dfc.R95))          # R95 along the x-axis (left → right)
        ys = sort(unique(dfc.σ))            # σ  along the y-axis (bottom → top)

        # dense matrix with x-first indexing (NO transpose later)
        lut = Dict((dfc.R95[i], dfc.σ[i]) => dfc.ΔA[i] for i in 1:nrow(dfc))
        Z   = [lut[(x, y)] for x in xs, y in ys]   # size = (length(xs), length(ys))

        # plot
        ax = Axis(fig[1,j],
                  title  = "ΔArea (AM−BAM) — R95×σ — $(g)",
                  xlabel = "R95 (diet redundancy)",
                  ylabel = (j == 1 ? "niche breadth σ" : ""))

        hm = heatmap!(ax, xs, ys, Z; interpolate=false)
        Colorbar(fig[2,j], hm, vertical=false, label="ΔArea")
    end

    save(joinpath(figdir, "02_HM_DeltaArea_R95_sigma_all_grids.png"), fig)
    display(fig)
end

# ----- Figure 03: Mediation (ΔArea vs P_suff), binned & consumers-only -----
# We only need res1 (C,align sweep), where Pmean/ΔAmean live
function fig03_mediation_binned(res1_by_grid)
    fig = Figure(size=(1500, 450))
    for (j,g) in enumerate(GRIDS)
        df = res1_by_grid[g]
        # consumers-only relation was ensured upstream when computing P_suff
        # (Pmean refers to consumers within A-kept cells). We bin P_suff:
        b = binmean(df.Pmean, df.ΔAmean; nbins=30)
        ax = Axis(fig[1,j], title="Mediation: ΔArea vs P_suff — $(g)",
                  xlabel="P_suff (consumers, within A-kept cells)", ylabel=(j==1 ? "ΔArea" : ""))
        lines!(ax, b.x, b.y)
        band!(ax, b.x, b.y .- b.se, b.y .+ b.se; transparency=true)
    end
    save(joinpath(figdir, "03_Mediation_binned_all_grids.png"), fig)
    display(fig)
end

# ----- Figure 04: Mechanism inserts (πA, P_suff, πA*(1-P_suff)) vs alignment -----
# Some CSVs may not have πAmean; detect safely.
function fig04_mechanism_inserts(res1_by_grid)
    fig = Figure(size=(1150, 820))
    row = 1
    for g in GRIDS
        df = res1_by_grid[g]
        have_piA = (:πAmean in propertynames(df)) || ("πAmean" in names(df))
        # sort by alignment
        sdf = sort(df, :align)
        # bin by alignment to smooth “wiggles”
        bπ = have_piA ? binmean(sdf.align, sdf.πAmean; nbins=35) : (x=Float64[],y=Float64[],se=Float64[])
        bP = binmean(sdf.align, sdf.Pmean; nbins=35)
        # approx ΔArea mediator: πA * (1-P_suff), use available columns
        approx = have_piA ? bπ.y .* (1 .- bP.y) : fill(NaN, length(bP.y))
        appse = have_piA ? sqrt.(bπ.se.^2 .* (1 .- bP.y).^2 .+ bP.se.^2 .* bπ.y.^2) : fill(NaN, length(bP.se))

        # πA
        ax1 = Axis(fig[row,1], title="Abiotic admission vs alignment — $(g)",
                   xlabel="alignment", ylabel="π_A (consumers)")
        if have_piA
            lines!(ax1, bπ.x, bπ.y); band!(ax1, bπ.x, bπ.y .- bπ.se, bπ.y .+ bπ.se; transparency=true)
        else
            text!(ax1, 0.1, 0.5, text="π_A unavailable in CSV", align=(:left,:center))
        end
        # P_suff
        ax2 = Axis(fig[row,2], title="Prey sufficiency vs alignment — $(g)",
                   xlabel="alignment", ylabel="P_suff")
        lines!(ax2, bP.x, bP.y); band!(ax2, bP.x, bP.y .- bP.se, bP.y .+ bP.se; transparency=true)
        # πA*(1-P_suff)
        ax3 = Axis(fig[row,3], title="Approx ΔArea from mediation — $(g)",
                   xlabel="alignment", ylabel="ΔArea ≈ π_A × (1 − P_suff)")
        if have_piA
            lines!(ax3, bP.x, approx)
            band!(ax3, bP.x, approx .- appse, approx .+ appse; transparency=true)
        else
            text!(ax3, 0.1, 0.5, text="π_A unavailable in CSV", align=(:left,:center))
        end
        row += 1
    end
    save(joinpath(figdir, "04_Mechanism_inserts_all_grids.png"), fig)
    display(fig)
end

# ----- Figure 05: ΔArea vs R95 conditional on σ deciles (lines) for each grid -----
function fig05_R95_lines_by_sigma_deciles(res2_by_grid; ndeciles=10)
    fig = Figure(size=(1500, 450))
    for (j,g) in enumerate(GRIDS)
        df = res2_by_grid[g]
        dfc = combine(groupby(df, [:R95, :σ]), :ΔAmean => mean => :ΔA)
        xs = sort(unique(dfc.R95))
        # make σ deciles:
        σvals = sort(unique(dfc.σ))
        cuts = [quantile(σvals, i/ndeciles) for i in 0:ndeciles]
        ax = Axis(fig[1,j], title="ΔArea vs R95 conditioned on σ deciles — $(g)",
                  xlabel="R95 (diet redundancy)", ylabel=(j==1 ? "ΔArea (AM−BAM)" : ""))
        objs = Lines[]
        labels = String[]
        for k in 1:ndeciles
            σlo, σhi = cuts[k], cuts[k+1]
            m = (dfc.σ .>= σlo) .& (dfc.σ .< σhi)
            if any(m)
                tmp = dfc[m, :]
                lut = Dict((tmp.R95[i]) => tmp.ΔA[i] for i in 1:nrow(tmp))
                ys = [get(lut, x, NaN) for x in xs]
                push!(objs, lines!(ax, xs, ys).plot)
                push!(labels, "σ∈[$(round(σlo,digits=3)),$(round(σhi,digits=3))]")
            end
        end
        Legend(fig[2,j], objs, labels; orientation=:horizontal, nbanks=2, framevisible=false)
    end
    save(joinpath(figdir, "05_R95_lines_by_sigma_deciles_all_grids.png"), fig)
    display(fig)
end

# ----- SI-1: Variability maps + fraction-above + ECDFs -----
function si01_variability_and_fraction(res1_by_grid, res2_by_grid; τmode=:p90)
    # helper: convert IQR to SD and then to CV
    iqr_to_sd(iqr) = iqr / 1.349
    safe_div(num, den) = num / (abs(den) < eps(Float64) ? eps(Float64) : den)

    # (A) CV maps over C×align
    figA = Figure(size=(1500, 520))
    for (j,g) in enumerate(GRIDS)
        df = res1_by_grid[g]
        xs = sort(unique(df.C)); ys = sort(unique(df.align))
        # SD from IQR band:
        sd_est = iqr_to_sd.(df.ΔAhi .- df.ΔAlo)
        # CV = SD / mean, guarded against tiny means
        cv = similar(sd_est)
        @inbounds for i in eachindex(sd_est)
            cv[i] = safe_div(sd_est[i], df.ΔAmean[i])
        end
        lut = Dict((df.C[i], df.align[i]) => cv[i] for i in 1:nrow(df))
        Z = [lut[(x,y)] for y in ys, x in xs]   # rows=y, cols=x
        ax = Axis(figA[1,j], title="CV(ΔArea) — C×align — $(g)",
                  xlabel="C", ylabel=(j==1 ? "align" : ""))
        hm = heatmap!(ax, xs, ys, Z)
        Colorbar(figA[2,j], hm, vertical=false, label="CV = SD/mean")
    end
    save(joinpath(figdir, "SI01A_CV_C_align.png"), figA)

    # (B) CV maps over R95×σ
    figB = Figure(size=(1500, 520))
    for (j,g) in enumerate(GRIDS)
        df = res2_by_grid[g]
        # collapse to one row per (R95,σ)
        dfc = combine(groupby(df, [:R95, :σ]),
                      :ΔAmean => mean => :μ,
                      [:ΔAlo, :ΔAhi] => ((lo,hi) -> mean(iqr_to_sd.(hi .- lo))) => :sd)
        xs = sort(unique(dfc.R95)); ys = sort(unique(dfc.σ))
        cv = [safe_div(dfc.sd[i], dfc.μ[i]) for i in 1:nrow(dfc)]
        lut = Dict((dfc.R95[i], dfc.σ[i]) => cv[i] for i in 1:nrow(dfc))
        Z = [lut[(x,y)] for y in ys, x in xs]
        ax = Axis(figB[1,j], title="CV(ΔArea) — R95×σ — $(g)",
                  xlabel="R95", ylabel=(j==1 ? "σ" : ""))
        hm = heatmap!(ax, xs, ys, Z)
        Colorbar(figB[2,j], hm, vertical=false, label="CV = SD/mean")
    end
    save(joinpath(figdir, "SI01B_CV_R95_sigma.png"), figB)

    # (C) fraction-above maps (unchanged logic; τ is grid-specific)
    figC = Figure(size=(1650, 520))
    for (j,g) in enumerate(GRIDS)
        df = res1_by_grid[g]
        τ = τmode === :p90 ? quantile(df.ΔAmean, 0.90) :
            τmode === :p75 ? quantile(df.ΔAmean, 0.75) : τmode # allow numeric
        xs = sort(unique(df.C)); ys = sort(unique(df.align))
        lut = Dict((df.C[i], df.align[i]) => (df.ΔAmean[i] > τ ? 1.0 : 0.0) for i in 1:nrow(df))
        Z = [lut[(x,y)] for y in ys, x in xs]
        ax = Axis(figC[1,j], title="Frac(ΔArea > $(round(τ, digits=3))) — C×align — $(g)",
                  xlabel="C", ylabel=(j==1 ? "align" : ""))
        hm = heatmap!(ax, xs, ys, Z; colormap=:viridis)
        Colorbar(figC[2,j], hm, vertical=false, label="fraction > τ", ticks=[0,0.5,1.0])
    end
    save(joinpath(figdir, "SI01C_fraction_above.png"), figC)

    display(figA); display(figB); display(figC)
end

# ----- SI-3: Mechanism (same content as Fig04 but it’s explicitly an SI figure) -----
# In the paper: Fig04 shows them as “inserts”; SI-3 is the full-page version.
function si03_mechanism_full(res1_by_grid)
    fig = fig04_mechanism_inserts(res1_by_grid) # reuse
    # save(joinpath(figdir, "SI03_mechanism_full.png"), fig)
    display(fig)
end

# =========================
# Run all figures
# =========================
fig01_heatmaps_C_align(res1_by_grid)
fig02_heatmaps_R_sigma(res2_by_grid)
fig03_mediation_binned(res1_by_grid)
fig04_mechanism_inserts(res1_by_grid)
fig05_R95_lines_by_sigma_deciles(res2_by_grid)
si01_variability_and_fraction(res1_by_grid, res2_by_grid; τmode=:p90)
si03_mechanism_full(res1_by_grid)

println("All figures written to ", figdir)
