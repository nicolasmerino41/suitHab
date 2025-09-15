# ----- helpers already in your code base assumed available -----
# make_grid, build_pool_from_axes, bam_from_axes, abiotic_matrix,
# assemble_BAM, bsh1_per_species, shapley_per_species,
# random_mask(rng, C, keep_frac), clustered_mask(rng, grid, keep_frac; nseeds=...),
# frontlike_mask(rng, grid, keep_frac; axis=:x, noise=...), consumer_mask(pool)

# colors/labels
const _glabel = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")
const _gcolor = Dict(:random=>RGBf(0.30,0.45,0.85), :clustered=>RGBf(0.95,0.60,0.15), :front=>RGBf(0.20,0.65,0.35))

_make_keepmask = function(rng::AbstractRNG, kind::Symbol, grid::Grid, keep_frac::Float64;
                          nseeds_cluster::Int, front_axis::Symbol, front_noise::Float64)
    if kind === :random
        return random_mask(rng, grid.C, keep_frac)
    elseif kind === :clustered
        return clustered_mask(rng, grid, keep_frac; nseeds=nseeds_cluster)
    elseif kind === :front
        return frontlike_mask(rng, grid, keep_frac; axis=front_axis, noise=front_noise)
    else
        error("Unknown hl_kind=$kind")
    end
end

# Reuse your AM vs ABM sweep to get both curves for a given geometry
function _sweep_AM_ABM_onegeom(; grid::Grid, S::Int, basal_frac::Float64,
        A_level::Symbol, B_level::Symbol, M_level::Symbol, hl_kind::Symbol,
        loss_fracs, seeds_pool, seeds_mask, sim_seed::Int,
        nseeds_cluster::Int, front_axis::Symbol, front_noise::Float64,
        τA::Float64, τocc::Float64, T_frac_on::Float64)

    xs   = Float64[]
    AMs  = Float64[]   # climate-only BSH (B ignored)
    ABMs = Float64[]   # full BSH (B on)

    for f in loss_fracs
        keepfrac = 1 - f
        vals_AM  = Float64[]
        vals_ABM = Float64[]

        for ps in seeds_pool, ms in seeds_mask
            # pool + abiotic
            rng_pool = MersenneTwister(hash((sim_seed,:pool,ps)))
            pool = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
            A    = abiotic_matrix(pool, grid)

            # parameters
            pars_on  = bam_from_axes(; B_level=B_level, M_level=M_level, τA=τA, τocc=τocc)
            pars_off = bam_from_axes(; B_level=:none,        M_level=M_level, τA=τA, τocc=τocc)

            # mask
            rng_mask = MersenneTwister(hash((sim_seed,:mask,ms,hl_kind,f)))
            keep = _make_keepmask(rng_mask, hl_kind, grid, keepfrac;
                                  nseeds_cluster=nseeds_cluster, front_axis=front_axis, front_noise=front_noise)

            # full
            P_on, _, _, _ = assemble_BAM(pool, grid, A, keep; bam=pars_on.bam, mp=pars_on.mp)
            bsh_on = bsh1_per_species(P_on, (A .>= τA), pool)[consumer_mask(pool)]
            push!(vals_ABM, mean(bsh_on))

            # climate-only
            P_off, _, _, _ = assemble_BAM(pool, grid, A, keep; bam=pars_off.bam, mp=pars_off.mp)
            bsh_off = bsh1_per_species(P_off, (A .>= τA), pool)[consumer_mask(pool)]
            push!(vals_AM, mean(bsh_off))
        end

        push!(xs, f); push!(AMs, mean(vals_AM)); push!(ABMs, mean(vals_ABM))
    end

    return (; x=xs, AM=AMs, ABM=ABMs)
end

# Shapley shares at one f* (averaged over seeds) for a geometry
function _shares_at_fstar(; grid::Grid, S::Int, basal_frac::Float64,
        A_level::Symbol, B_level::Symbol, M_level::Symbol, hl_kind::Symbol,
        loss_pick::Float64, seeds_pool, seeds_mask, sim_seed::Int,
        nseeds_cluster::Int, front_axis::Symbol, front_noise::Float64,
        τA::Float64, τocc::Float64)

    keepfrac = 1 - loss_pick
    ΔA_vec = Float64[]; ΔB_vec = Float64[]; ΔM_vec = Float64[]

    for ps in seeds_pool, ms in seeds_mask
        rng_pool = MersenneTwister(hash((sim_seed,:pool,ps)))
        pool     = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
        A        = abiotic_matrix(pool, grid)
        pars     = bam_from_axes(; B_level, M_level, τA, τocc)
        bam, mp  = pars.bam, pars.mp

        rng_mask = MersenneTwister(hash((sim_seed,:mask,ms,hl_kind,loss_pick)))
        keep = _make_keepmask(rng_mask, hl_kind, grid, keepfrac;
                              nseeds_cluster=nseeds_cluster, front_axis=front_axis, front_noise=front_noise)

        shap, _, _ = shapley_per_species(pool, grid, A, trues(grid.C), keep; bam=bam, mp=mp)
        cons = consumer_mask(pool)
        push!(ΔA_vec, mean(abs, getfield.(shap[cons], :dA)))
        push!(ΔB_vec, mean(abs, getfield.(shap[cons], :dB)))
        push!(ΔM_vec, mean(abs, getfield.(shap[cons], :dM)))
    end

    ΔA, ΔB, ΔM = mean(ΔA_vec), mean(ΔB_vec), mean(ΔM_vec)
    denom = max(ΔA + ΔB + ΔM, eps())
    return (; S_A=ΔA/denom, S_B=ΔB/denom, S_M=ΔM/denom)
end

# ternary helpers
_ternary_xy(sA,sB,sM) = (sB + 0.5*sM, (√3/2)*sM)
function _ternary_frame!(ax; title="")
    xs = [0.0,1.0,0.5,0.0]; ys=[0.0,0.0,√3/2,0.0]
    lines!(ax, xs, ys; color=:black)
    text!(ax, 0.02, -0.04, text="Abiotic", align=(:left,:bottom))
    text!(ax, 0.98, -0.04, text="Biotic",  align=(:right,:bottom))
    text!(ax, 0.50,  √3/2+0.03, text="Movement", align=(:center,:bottom))
    ax.title = title
    hidespines!(ax, :r, :t)
end

# --- degree-preserving rewiring (prey must be lighter than predator) ---
function rewire_metaweb_inplace!(rng::AbstractRNG, pool::SpeciesPool)
    order = sortperm(pool.masses)            # light→heavy
    rank  = zeros(Int, pool.S); rank[order] .= 1:pool.S
    for s in 1:pool.S
        pool.basal[s] && continue
        k = length(pool.E[s]); k==0 && continue
        cand = order[1:rank[s]-1]
        if isempty(cand)
            pool.E[s] = Int[]; continue
        end
        if length(cand) <= k
            pool.E[s] = copy(cand)
        else
            pool.E[s] = sample(rng, cand, k; replace=false)
        end
    end
    return pool
end

# --- helper: sweep AM (B off) & ABM (B on) for one geometry; placebo version rewires pool ---
function _sweep_AM_ABM_onegeom_placebo(; grid::Grid, S::Int, basal_frac::Float64,
        A_level::Symbol, B_level::Symbol, M_level::Symbol, hl_kind::Symbol,
        loss_fracs, seeds_pool, seeds_mask, sim_seed::Int,
        nseeds_cluster::Int, front_axis::Symbol, front_noise::Float64,
        τA::Float64, τocc::Float64)

    xs   = Float64[]
    AMs  = Float64[]   # climate-only (B off)
    ABMs = Float64[]   # full (B on) but with rewired metaweb

    for f in loss_fracs
        keepfrac = 1 - f
        vals_AM  = Float64[]
        vals_ABM = Float64[]

        for ps in seeds_pool, ms in seeds_mask
            rng_pool = MersenneTwister(hash((sim_seed,:pool,ps)))
            pool     = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
            # rewire AFTER building pool/metaweb
            rewire_metaweb_inplace!(MersenneTwister(hash((sim_seed,:rewire,ps))), pool)

            A    = abiotic_matrix(pool, grid)
            pars_on  = bam_from_axes(; B_level=B_level, M_level=M_level, τA=τA, τocc=τocc)
            pars_off = bam_from_axes(; B_level=:none,        M_level=M_level, τA=τA, τocc=τocc)

            rng_mask = MersenneTwister(hash((sim_seed,:mask,ms,hl_kind,f)))
            keep = if hl_kind === :random
                random_mask(rng_mask, grid.C, keepfrac)
            elseif hl_kind === :clustered
                clustered_mask(rng_mask, grid, keepfrac; nseeds=nseeds_cluster)
            else
                frontlike_mask(rng_mask, grid, keepfrac; axis=front_axis, noise=front_noise)
            end

            # Full BSH with rewired metaweb
            P_on, _, _, _  = assemble_BAM(pool, grid, A, keep; bam=pars_on.bam, mp=pars_on.mp)
            bsh_on         = bsh1_per_species(P_on, (A .>= τA), pool)[consumer_mask(pool)]
            push!(vals_ABM, mean(bsh_on))

            # Climate-only (unchanged by rewiring)
            P_off, _, _, _ = assemble_BAM(pool, grid, A, keep; bam=pars_off.bam, mp=pars_off.mp)
            bsh_off        = bsh1_per_species(P_off, (A .>= τA), pool)[consumer_mask(pool)]
            push!(vals_AM, mean(bsh_off))
        end

        push!(xs, f); push!(AMs, mean(vals_AM)); push!(ABMs, mean(vals_ABM))
    end
    return (; x=xs, AM=AMs, ABM=ABMs)
end

"""
did_geometry_dashboard(; grid, S, basal_frac, A_level, B_level, M_level,
                       loss_fracs, loss_pick, seeds_pool, seeds_mask, sim_seed,
                       nseeds_cluster, front_axis, front_noise, τA, τocc)

→ makes a 2×2 figure:
  (1) ΔBSH vs loss, B OFF   (A×M)
  (2) ΔBSH vs loss, B ON    (A×B×M)
  (3) DiD(f) = Δ_bio - Δ_non  by geometry
  (4) Elasticity triangle at f* annotated with S_B for each geometry
"""
# --- UPDATED did_geometry_dashboard: adds placebo DiD overlay in panel (3) ---
function did_geometry_dashboard(; grid::Grid, S::Int=150, basal_frac::Float64=0.45,
        A_level::Symbol=:intermediate, B_level::Symbol=:soft, M_level::Symbol=:on,
        loss_fracs=0.2:0.1:0.8, loss_pick::Float64=0.6,
        seeds_pool=1:6, seeds_mask=1:6, sim_seed::Int=1234,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        τA::Float64=0.5, τocc::Float64=0.35)

    geoms = (:random,:clustered,:front)
    curves = Dict{Symbol,NamedTuple}()
    curves_pl = Dict{Symbol,NamedTuple}()

    for hk in geoms
        c  = _sweep_AM_ABM_onegeom(; grid,S,basal_frac,A_level,B_level,M_level,
                                   hl_kind=hk, loss_fracs,seeds_pool,seeds_mask,sim_seed,
                                   nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on=0.0)
        cp = _sweep_AM_ABM_onegeom_placebo(; grid,S,basal_frac,A_level,B_level,M_level,
                                   hl_kind=hk, loss_fracs,seeds_pool,seeds_mask,sim_seed,
                                   nseeds_cluster,front_axis,front_noise,τA,τocc)

        Δnon  = c.AM  .- first(c.AM)
        Δbio  = c.ABM .- first(c.ABM)
        did   = Δbio .- Δnon

        Δnon_p = cp.AM .- first(cp.AM)     # same AM baseline
        Δbio_p = cp.ABM .- first(cp.ABM)
        did_p  = Δbio_p .- Δnon_p

        curves[hk]    = (; x=c.x, Δnon, Δbio, did, AM=c.AM, ABM=c.ABM)
        curves_pl[hk] = (; x=cp.x, did=did_p)
    end

    # Elasticity shares at f* (B ON)
    shares = Dict{Symbol,NamedTuple}()
    for hk in geoms
        shares[hk] = _shares_at_fstar(; grid,S,basal_frac,A_level,
                                      B_level=B_level, M_level, hl_kind=hk,
                                      loss_pick, seeds_pool,seeds_mask,sim_seed,
                                      nseeds_cluster,front_axis,front_noise,τA,τocc)
    end

    # ---------- Figure ----------
    fig = Figure(; size=(1400, 900))
    Label(fig[0,:],
          @sprintf("HL effect with vs without biotic — A=%s, B(on)=%s, M=%s   (loss* = %.2f, placebo dashed)",
                   String(A_level), String(B_level), String(M_level), loss_pick);
          fontsize=18, padding=(0,0,8,0))

    # (1) ΔBSH vs loss, B OFF
    ax1 = Axis(fig[1,1], xlabel="Area lost (fraction)", ylabel="ΔBSH (mean cons.)",
               title="Climate-only (B OFF: A×M)")
    for hk in geoms
        r = curves[hk]
        lines!(ax1, r.x, r.Δnon; color=_gcolor[hk], label=_glabel[hk])
    end
    axislegend(ax1; position=:lb, framevisible=false, labelsize=10)

    # (2) ΔBSH vs loss, B ON
    ax2 = Axis(fig[1,2], xlabel="Area lost (fraction)", ylabel="ΔBSH (mean cons.)",
               title=@sprintf("Full BSH (B ON: A×B×M; B=%s)", String(B_level)))
    for hk in geoms
        r = curves[hk]
        lines!(ax2, r.x, r.Δbio; color=_gcolor[hk], label=_glabel[hk])
    end
    axislegend(ax2; position=:lb, framevisible=false, labelsize=10)

    # (3) DiD with placebo overlay
    ax3 = Axis(fig[2,1], xlabel="Area lost (fraction)", ylabel="DiD(f) = Δbio − Δnon",
               title="Incremental impact of HL due to biotic (solid) vs placebo (dashed)")
    for hk in geoms
        r  = curves[hk]
        rp = curves_pl[hk]
        lines!(ax3, r.x,  r.did;  color=_gcolor[hk], label=_glabel[hk])
        lines!(ax3, rp.x, rp.did; color=_gcolor[hk], linestyle=:dash)
    end
    hlines!(ax3, [0.0]; linestyle=:dot, color=:gray)
    axislegend(ax3; position=:lb, framevisible=false, labelsize=10)

    # (4) Elasticity triangle at f*
    ax4 = Axis(fig[2,2])
    xs = [0.0,1.0,0.5,0.0]; ys=[0.0,0.0,√3/2,0.0]
    lines!(ax4, xs, ys; color=:black)
    text!(ax4, 0.02, -0.04, text="Abiotic", align=(:left,:bottom))
    text!(ax4, 0.98, -0.04, text="Biotic",  align=(:right,:bottom))
    text!(ax4, 0.50,  √3/2+0.03, text="Movement", align=(:center,:bottom))
    hidespines!(ax4, :r, :t)
    ax4.title = @sprintf("Elasticity shares at f* = %.2f (B ON)", loss_pick)
    for hk in geoms
        s = shares[hk]
        x,y = (s.S_B + 0.5*s.S_M, (√3/2)*s.S_M)
        scatter!(ax4, [x],[y]; color=_gcolor[hk], markersize=14, label=_glabel[hk])
        text!(ax4, x, y; text=@sprintf(" B=%.2f ", s.S_B), align=(:left,:bottom),
              fontsize=10, color=_gcolor[hk])
    end
    axislegend(ax4; position=:rb, framevisible=false)

    display(fig)
    return (; fig, curves, curves_placebo=curves_pl, shares)
end

grid = make_grid(60,60; seed=42)
_ = did_geometry_dashboard(; grid,
      S=150, basal_frac=0.45,
      A_level=:divergent, B_level=:strong, M_level=:on,
      loss_fracs=0.2:0.1:0.8, loss_pick=0.60,
      seeds_pool=1:6, seeds_mask=1:6, sim_seed=1234,
      nseeds_cluster=6, front_axis=:x, front_noise=0.04,
      τA=0.5, τocc=0.35)