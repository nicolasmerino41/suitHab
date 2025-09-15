# --- adjacency and component mask (4-neighbour) -----------------------------
function compmask4(keep::BitVector, nx::Int, ny::Int, T::Int)
    C = length(keep)
    mark = falses(C)
    good = falses(C)
    # neighbour helper
    neighbors4(ix) = begin
        i = ((ix - 1) % nx) + 1
        j = ((ix - 1) ÷ nx) + 1
        out = Int[]
        if i > 1;  push!(out, ix-1);   end
        if i < nx; push!(out, ix+1);   end
        if j > 1;  push!(out, ix-nx);  end
        if j < ny; push!(out, ix+nx);  end
        out
    end
    for s in 1:C
        (keep[s] && !mark[s]) || continue
        # BFS this component
        q = [s]; ptr = 1
        comp = Int[]
        mark[s] = true
        while ptr <= length(q)
            v = q[ptr]; ptr += 1
            push!(comp, v)
            for nb in neighbors4(v)
                if keep[nb] && !mark[nb]
                    mark[nb] = true
                    push!(q, nb)
                end
            end
        end
        if length(comp) >= T
            @inbounds for v in comp; good[v] = true; end
        end
    end
    return good
end

# --- per-species summaries ---------------------------------------------------
# A_s  : fraction of ALL cells climatically suitable (on the kept grid for post-HL)
# B_s  : among A-cells, share with ≥1 prey present (movement OFF for assembling prey)
# M_s  : among A-cells, share that are inside components of size ≥ T (if M:on; else 1)
function summarize_ABM_per_species(pool, grid, A::Matrix{Float64}, keep::BitVector;
                                   B_level::Symbol, M_level::Symbol,
                                   τA::Float64, τocc::Float64, T_frac_on::Float64)
    C = grid.C
    # climate pass on kept cells only
    Z = falses(pool.S, C)
    for s in 1:pool.S, i in 1:C
        if keep[i]
            Z[s,i] = A[s,i] ≥ τA
        end
    end
    # assemble once with movement OFF to get prey layers
    bam  = bam_from_axes(; B_level=B_level, M_level=:off, τA=τA, τocc=τocc).bam
    mp   = MovementParams(; mode=:none, T=8)
    P, _, _, _ = assemble_BAM(pool, grid, A, keep; bam=bam, mp=mp)

    # movement component mask (species-independent)
    Mcomp = if M_level === :on
        T = max(1, round(Int, T_frac_on * C))
        compmask4(keep, grid.nx, grid.ny, T)
    else
        trues(C)
    end

    A_s = zeros(Float64, pool.S)
    B_s = zeros(Float64, pool.S)
    M_s = zeros(Float64, pool.S)

    for s in 1:pool.S
        if pool.basal[s]; continue; end  # consumers only (you can include basal if you want)
        denomA = count(@view Z[s, :])   # number of suitable cells on kept grid
        A_s[s] = denomA / C
        if denomA == 0
            B_s[s] = 0.0; M_s[s] = 0.0
            continue
        end
        # B: share of suitable cells with ≥1 prey present
        okB = 0
        for i in 1:C
            if Z[s,i]
                has = false
                @inbounds for q in pool.E[s]
                    if P[q,i]; has = true; break; end
                end
                okB += has ? 1 : 0
            end
        end
        B_s[s] = okB / denomA
        # M: share of suitable cells that are in big components
        okM = 0
        for i in 1:C
            if Z[s,i] && Mcomp[i]; okM += 1; end
        end
        M_s[s] = okM / denomA
    end
    return (; A_s, B_s, M_s)
end

# --- simple HL keep masks ----------------------------------------------------
function front_keep_mask(grid::Grid; keep_frac::Float64, axis::Symbol=:x, noise::Float64=0.03)
    t = axis === :x ? grid.xy[1,:] : grid.xy[2,:]
    thr = quantile(t, keep_frac)
    keep = t .>= (thr .- noise .* (rand(length(t)) .- 0.5))
    return BitVector(keep)
end

# --- Shapley decomposition of ΔF over (A,B,M) for given β, μ ----------------
G(A,B,M,β,μ) = A * (B^β) * (M^μ)

function shapley_deltaF(A0,B0,M0, A1,B1,M1; β::Float64=1.0, μ::Float64=1.0)
    perms = ((:A,:B,:M), (:A,:M,:B), (:B,:A,:M), (:B,:M,:A), (:M,:A,:B), (:M,:B,:A))
    dA = 0.0; dB = 0.0; dM = 0.0
    for p in perms
        Acur,Bcur,Mcur = A0,B0,M0
        v0 = G(Acur,Bcur,Mcur,β,μ)
        # first
        if p[1] === :A; Acur = A1
        elseif p[1] === :B; Bcur = B1
        else; Mcur = M1; end
        v1 = G(Acur,Bcur,Mcur,β,μ)
        if p[1] === :A; dA += (v1 - v0)
        elseif p[1] === :B; dB += (v1 - v0)
        else; dM += (v1 - v0); end
        # second
        if p[2] === :A; Acur = A1
        elseif p[2] === :B; Bcur = B1
        else; Mcur = M1; end
        v2 = G(Acur,Bcur,Mcur,β,μ)
        if p[2] === :A; dA += (v2 - v1)
        elseif p[2] === :B; dB += (v2 - v1)
        else; dM += (v2 - v1); end
        # third
        if p[3] === :A; Acur = A1
        elseif p[3] === :B; Bcur = B1
        else; Mcur = M1; end
        v3 = G(Acur,Bcur,Mcur,β,μ)
        if p[3] === :A; dA += (v3 - v2)
        elseif p[3] === :B; dB += (v3 - v2)
        else; dM += (v3 - v2); end
    end
    k = 1/6
    return (k*dA, k*dB, k*dM)  # exact split; sums to ΔF
end

"""
plot_beta_sweep(; A_level, B_level, M_level, hl_kind, loss, μ, β_grid, nx, ny, S, basal_frac, seeds_pool, seeds_mask, τA, τocc, T_frac_on, front_axis, front_noise)

Produces a line plot of ΔF components vs β (A, B, M), with a horizontal dashed line for ΔBSH (mean consumers).
"""
function plot_beta_sweep(; A_level::Symbol=:divergent, B_level::Symbol=:soft, M_level::Symbol=:on,
        hl_kind::Symbol=:front, loss::Float64=0.6, μ::Float64=1.0, β_grid = 0.0:0.1:1.5,
        nx::Int=60, ny::Int=60, S::Int=150, basal_frac::Float64=0.45,
        seeds_pool=1:3, seeds_mask=1:3, τA::Float64=0.5, τocc::Float64=0.35,
        T_frac_on::Float64=0.02, front_axis::Symbol=:x, front_noise::Float64=0.03, sim_seed::Int=1234)

    grid = make_grid(nx, ny; seed=42)
    keep_frac = 1.0 - loss

    # accumulate means across seeds
    dAβ = zeros(length(β_grid)); dBβ = similar(dAβ); dMβ = similar(dAβ)
    ΔBSH_vals = Float64[]

    for ps in seeds_pool, ms in seeds_mask
        rng = MersenneTwister(hash((sim_seed,:pool,ps)))
        pool = build_pool_from_axes(rng; S, basal_frac, A_level, B_level)
        A = abiotic_matrix(pool, grid)

        # baseline (no HL): keep everything
        keep0 = trues(grid.C)
        base = summarize_ABM_per_species(pool, grid, A, keep0; B_level, M_level, τA, τocc, T_frac_on)

        # HL keep mask
        keep1 = begin
            if hl_kind === :random
                random_mask(grid.C, keep_frac; seed=ms)
            elseif hl_kind === :clustered
                clustered_mask(grid, keep_frac; nseeds=6, seed=ms)
            elseif hl_kind === :front
                front_keep_mask(grid; keep_frac=keep_frac, axis=front_axis, noise=front_noise)
            else
                error("Unknown hl_kind=$hl_kind")
            end
        end
        post = summarize_ABM_per_species(pool, grid, A, keep1; B_level, M_level, τA, τocc, T_frac_on)

        # ΔBSH (ground truth) with your full BAM assembly (movement on/off as M_level)
        # we re-use your sweep machinery for a single f:
        rr = sweep_dBSH_axes(; grid,S,basal_frac,A_level,B_level,M_level,hl_kind,
                              loss_fracs=[loss], seeds_pool=[ps], seeds_mask=[ms], sim_seed,
                              front_axis, front_noise, τA, τocc, T_frac_on)
        push!(ΔBSH_vals, rr.y[1])

        # aggregate means over consumers (basal ignored)
        cons = .!pool.basal
        A0 = mean(base.A_s[cons]);  B0 = mean(base.B_s[cons]);  M0 = mean(base.M_s[cons])
        A1 = mean(post.A_s[cons]);  B1 = mean(post.B_s[cons]);  M1 = mean(post.M_s[cons])

        for (k,β) in enumerate(β_grid)
            dA, dB, dM = shapley_deltaF(A0,B0,M0, A1,B1,M1; β=β, μ=μ)
            dAβ[k] += dA; dBβ[k] += dB; dMβ[k] += dM
        end
    end

    nrep = length(seeds_pool)*length(seeds_mask)
    dAβ ./= nrep; dBβ ./= nrep; dMβ ./= nrep
    ΔBSH_mean = mean(ΔBSH_vals)

    # plot
    fig = Figure(; size=(800, 420))
    ax = Axis(fig[1,1],
              title=@sprintf("β-sweep — A=%s, B=%s, M=%s | %s, loss=%.2f", string(A_level), string(B_level), string(M_level), string(hl_kind), loss),
              xlabel="β (biotic requirement in F)", ylabel="ΔF (mean cons.)")

    lines!(ax, collect(β_grid), dAβ; label="Abiotic", linewidth=2)
    lines!(ax, collect(β_grid), dBβ; label="Biotic", linewidth=2)
    lines!(ax, collect(β_grid), dMβ; label="Movement", linewidth=2)
    hlines!(ax, [ΔBSH_mean]; color=:gray, linestyle=:dash, label="ΔBSH (ground truth)")

    axislegend(ax; position=:rb, framevisible=false)
    display(fig)
    return fig, (; β=collect(β_grid), dA=dAβ, dB=dBβ, dM=dMβ, ΔBSH=ΔBSH_mean)
end

# Example: show how the biotic share grows under a front when B is strong
_ = plot_beta_sweep(; A_level=:divergent, B_level=:strong, M_level=:on,
                    hl_kind=:front, loss=0.6, β_grid=0.0:0.1:1.5)
