# ============================================================
# 13) MECHANISM VIGNETTES (4 param cells)
#     Fixed: netfamily=:heavytail, regime="Narrow + HighVar"
#     Cells: (C low/high) × (r low/high)
#     Outputs per cell:
#       - mismatch distribution with mean + q90
#       - prey-support vs mismatch scatter
#       - maps (A vs AB) for median & q90 species
# ============================================================

# --- choose the fixed regime index (Narrow + HighVar is regimes[2] in your list)
const VIGNETTE_REGIME_INDEX = 2
const VIGNETTE_NET = :heavytail

# Helper: nearest index in a vector
nearest_index(vals::Vector{Float64}, target::Float64) = argmin(abs.(vals .- target))

# Pick “low/high” targets (tweak if you want)
C_low_target  = CONNECTANCE_RANGE[1]
C_high_target = CONNECTANCE_RANGE[2]
r_low_target  = CORR_RANGE[1]
r_high_target = CORR_RANGE[2]

cind_lo = nearest_index(Cvals, C_low_target)
cind_hi = nearest_index(Cvals, C_high_target)
rind_lo = nearest_index(Rvals, r_low_target)
rind_hi = nearest_index(Rvals, r_high_target)

vignette_cells = [
    (cind_lo, rind_lo, "lowC_lowr"),
    (cind_hi, rind_hi, "highC_highr"),
    (cind_hi, rind_lo, "highC_lowr"),
    (cind_lo, rind_hi, "lowC_highr")
]

# ---------------------------
# Per-species diagnostics
# ---------------------------

# Per-species mismatch vector (consumers-only, excluding |A_i|=0 as before)
function per_species_mismatch(A::Vector{BitVector}, AB::Vector{BitVector}, basal_mask::BitVector)
    m = fill(NaN, S)
    for i in 1:S
        basal_mask[i] && continue
        Ai = A[i]
        count(Ai) == 0 && continue
        inter = count(Ai .& AB[i])
        uni   = count(Ai .| AB[i])
        J = uni == 0 ? 1.0 : (inter / uni)
        m[i] = 1 - J
    end
    return m
end

# Per-species prey support:
# s_i = |A_i ∩ (⋃_{j in prey(i)} A_j)| / |A_i|
function per_species_prey_support(A_raw::Vector{BitVector}, prey::Vector{Vector{Int}}, basal_mask::BitVector)
    s = fill(NaN, S)
    for i in 1:S
        basal_mask[i] && continue
        Ai = A_raw[i]
        a = count(Ai)
        (a == 0 || isempty(prey[i])) && continue
        u = BitVector(falses(NCELLS))
        @inbounds for j in prey[i]
            u .|= A_raw[j]
        end
        s[i] = count(Ai .& u) / a
    end
    return s
end

# Small utility: choose exemplar indices (median + q90 among valid consumers)
function pick_exemplars(mismatch_vec::Vector{Float64}, basal_mask::BitVector)
    idx = [i for i in 1:S if !basal_mask[i] && isfinite(mismatch_vec[i])]
    isempty(idx) && return (nothing, nothing)
    vals = mismatch_vec[idx]
    ord = sortperm(vals)
    med_i = idx[ord[clamp(Int(round(0.50*length(ord))), 1, length(ord))]]
    q90_i = idx[ord[clamp(Int(round(0.90*length(ord))), 1, length(ord))]]
    return (med_i, q90_i)
end

# ---------------------------
# A single replicate that RETURNS the objects we need
# (this is like simulate_one!, but it exposes prey, A_raw, A, AB)
# ---------------------------
function simulate_one_with_details!(
    rng::AbstractRNG,
    ws::CCWorkspace,
    envkind::Symbol,
    netfamily::Symbol,
    regime::BreadthRegime,
    C::Float64,
    target_r::Float64
)
    nb, basal_mask, consumers = consumers_and_basal()

    # Environment
    E = make_environment(rng, envkind)

    # Network
    prey = build_metaweb_heavytail(rng, C, basal_mask)

    # Niche parameters
    σ = draw_sigmas(rng, regime)
    μ = Vector{Float64}(undef, S)
    for i in 1:nb
        μ[i] = rand(rng) * (E_MAX - E_MIN) + E_MIN
    end
    for i in (nb+1):S
        μ[i] = rand(rng) * (E_MAX - E_MIN) + E_MIN
    end

    achieved_r = assign_mus_with_target_corr!(rng, μ, prey, basal_mask, target_r)

    # A (raw climatic)
    A_raw = Vector{BitVector}(undef, S)
    for i in 1:S
        A_raw[i] = suitability_mask_1d(E, μ[i], σ[i], SUIT_THRESH)
    end

    # AB (raw)
    AB_raw = fixed_point_AB(A_raw, prey, basal_mask)

    # Connectivity filter post-hoc
    A  = Vector{BitVector}(undef, S)
    AB = Vector{BitVector}(undef, S)
    for i in 1:S
        A[i]  = apply_connectivity_filter(ws, A_raw[i], Emin_patch)
        AB[i] = apply_connectivity_filter(ws, AB_raw[i], Emin_patch)
    end

    return (A_raw=A_raw, A=A, AB=AB, prey=prey, basal_mask=basal_mask,
            achieved_r=achieved_r)
end

# ---------------------------
# Plot helpers (Makie)
# ---------------------------
function hist_with_mean_q90(vals::Vector{Float64}; title="", xlabel="", outfile=nothing)
    v = vals[isfinite.(vals)]

    fig = Figure(size = (900, 650))
    
    # --- MAIN AXIS (creates row 1) ---
    ax = Axis(
        fig[1, 1];
        title  = title,
        xlabel = xlabel,
        ylabel = "Count",
        limits = (-0.01, 1.01, nothing, nothing)
    )

    hist!(ax, v, bins = 30)

    # --- ANNOTATION (THIS CREATES ROW 2) ---
    if !isempty(v)
        μ = mean(v)
        q = quantile(v, 0.90)

        vlines!(ax, [μ], linewidth = 3)
        vlines!(ax, [q], linewidth = 3, linestyle = :dash)

        Label(
            fig[2, 1],
            "mean=$(round(μ,digits=3)) | q90=$(round(q,digits=3)) | n=$(length(v))",
            fontsize = 12,
            halign = :center
        )
    end

    # --- NOW it is SAFE to resize ---
    colsize!(fig.layout, 1, Relative(1))
    rowsize!(fig.layout, 1, Relative(1))
    rowsize!(fig.layout, 2, Auto())

    outfile !== nothing && save(outfile, fig)
    return fig
end

function scatter_support_vs_mismatch(
    support::Vector{Float64},
    mismatch::Vector{Float64};
    title = "",
    outfile = nothing
)
    good = isfinite.(support) .& isfinite.(mismatch)
    x = support[good]
    y = mismatch[good]

    fig = Figure(size = (900, 650))

    ax = Axis(fig[1, 1];
        title  = title,
        xlabel = "Prey support in A (fraction)",
        ylabel = "Mismatch 1 - J(A,AB)",
        limits = (0, 1, 0, 1)   # 🔑 THIS fixes the layout
    )

    scatter!(
        ax,
        x,
        y;
        markersize = 8,
        color = (:dodgerblue, 0.6)
    )

    ylims!(ax, -0.05, 1.05)
    xlims!(ax, -0.05, 1.05)

    # Label(fig[2, 1], "n=$(length(x))", fontsize = 12)

    if outfile !== nothing
        save(outfile, fig)
    end

    return fig
end


# map plot for one species
function map_A_AB(Ai::BitVector, ABi::BitVector; title="", outfile=nothing)
    # reshape to NY×NX for plotting
    A2  = reshape(Float32.(Ai), NX, NY)'
    AB2 = reshape(Float32.(ABi), NX, NY)'

    fig = Figure(size=(1100, 450))
    ax1 = Axis(fig[1,1], title=title*" — A",  xticksvisible=false, yticksvisible=false)
    ax2 = Axis(fig[1,2], title=title*" — AB", xticksvisible=false, yticksvisible=false)

    heatmap!(ax1, A2)
    heatmap!(ax2, AB2)

    outfile !== nothing && save(outfile, fig)
    return fig
end

# ---------------------------
# Run vignettes
# ---------------------------

vdir = joinpath(OUTDIR, "vignettes_heavytail_narrowHighVar")
isdir(vdir) || mkpath(vdir)

for env in ENVKINDS
    envtag = String(env)

    for (cind, rind, tag) in vignette_cells
        C = Cvals[cind]
        r = Rvals[rind]

        # Use ONE representative replicate for the spatial maps,
        # but average distributions over several reps to stabilize hist/scatter.
        # (You can set reps_dist = NREP if you want.)
        reps_dist = max(6, min(12, NREP))

        tid = Threads.threadid()
        ws = WSS[tid]
        rng = MersenneTwister(BASE_SEED + 99_998 + 17*tid)

        # collect per-species vectors across reps (stacked)
        all_m = Float64[]
        all_s = Float64[]

        # store one "example replicate" for maps
        example = nothing

        for rep in 1:reps_dist
            seed = BASE_SEED + 1_000_000*(env==:autocorr ? 2 : 1) + 10_000*rep + 100*cind + rind
            Random.seed!(rng, seed)

            out = simulate_one_with_details!(rng, ws, env, VIGNETTE_NET, regimes[VIGNETTE_REGIME_INDEX], C, r)

            m = per_species_mismatch(out.A, out.AB, out.basal_mask)
            s = per_species_prey_support(out.A_raw, out.prey, out.basal_mask)

            append!(all_m, m[isfinite.(m)])
            # keep matched pairs only for scatter
            good = isfinite.(m) .& isfinite.(s)
            append!(all_s, s[good])

            example === nothing && (example = (out=out, m=m))
        end

        # --- distribution fig
        f1 = hist_with_mean_q90(all_m;
            title="Mismatch distribution ($envtag) — $tag  (C=$(round(C,digits=3)), r=$(round(r,digits=2)))",
            xlabel="1 - J(A,AB)",
            outfile=joinpath(vdir, "dist_$(envtag)_$(tag).png")
        )
        display(f1)

        # --- mediator scatter (using the example replicate’s per-species vectors is also fine,
        # but we’ll use the latest example’s vectors for matching support+mismatch)
        out_ex = example.out
        m_ex = example.m
        s_ex = per_species_prey_support(out_ex.A_raw, out_ex.prey, out_ex.basal_mask)

        f2 = scatter_support_vs_mismatch(s_ex, m_ex;
            title="Prey support vs mismatch ($envtag) — $tag  (C=$(round(C,digits=3)), r=$(round(r,digits=2)))",
            outfile=joinpath(vdir, "support_vs_mismatch_$(envtag)_$(tag).png")
        )
        display(f2)

        # --- spatial maps for median and q90 consumers
        med_i, q90_i = pick_exemplars(m_ex, out_ex.basal_mask)
        if med_i !== nothing
            f3 = map_A_AB(out_ex.A[med_i], out_ex.AB[med_i];
                title="Species $(med_i) (median mismatch) — $envtag $tag",
                outfile=joinpath(vdir, "maps_median_$(envtag)_$(tag).png")
            )
            display(f3)
        end
        if q90_i !== nothing
            f4 = map_A_AB(out_ex.A[q90_i], out_ex.AB[q90_i];
                title="Species $(q90_i) (q90 mismatch) — $envtag $tag",
                outfile=joinpath(vdir, "maps_q90_$(envtag)_$(tag).png")
            )
            display(f4)
        end

        @info "Saved vignette set" env=envtag tag=tag C=C r=r outdir=vdir
    end
end
