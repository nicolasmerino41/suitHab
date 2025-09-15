
# Map geometry symbol to label/color
const _glabel = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")
const _gcolor = Dict(:random=>RGBf(0.30,0.45,0.85), :clustered=>RGBf(0.95,0.60,0.15), :front=>RGBf(0.20,0.65,0.35))

# Make a keep-mask with RNG-first signatures
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

"""
geometry_elasticity_shares_all_geoms(; grid, S, basal_frac, A_level, B_level, M_level,
                                     loss_pick, seeds_pool, seeds_mask, sim_seed,
                                     nseeds_cluster, front_axis, front_noise, τA, τocc)

Returns Dict{Symbol,NamedTuple} keyed by :random/:clustered/:front with fields:
  ΔA, ΔB, ΔM    (mean over consumers & seeds)
  S_A, S_B, S_M (elasticity shares = |Δ•| / sum |Δ•|)
  ΔF            (mean ΔF)
"""
function geometry_elasticity_shares_all_geoms(; grid::Grid, S::Int, basal_frac::Float64,
        A_level::Symbol, B_level::Symbol, M_level::Symbol,
        loss_pick::Float64, seeds_pool, seeds_mask, sim_seed::Int,
        nseeds_cluster::Int, front_axis::Symbol, front_noise::Float64,
        τA::Float64, τocc::Float64)

    keepfrac = 1 - loss_pick
    geoms = (:random, :clustered, :front)
    out = Dict{Symbol,NamedTuple}()

    for hk in geoms
        ΔA_vec = Float64[]; ΔB_vec = Float64[]; ΔM_vec = Float64[]; ΔF_vec = Float64[]

        for ps in seeds_pool, ms in seeds_mask
            rng_pool = MersenneTwister(hash((sim_seed,:pool,ps)))
            pool     = build_pool_from_axes(rng_pool; S, basal_frac, A_level, B_level)
            A        = abiotic_matrix(pool, grid)
            pars     = bam_from_axes(; B_level, M_level, τA, τocc)
            bam, mp  = pars.bam, pars.mp

            rng_mask = MersenneTwister(hash((sim_seed,:mask,ms,hk,loss_pick)))
            keep     = _make_keepmask(rng_mask, hk, grid, keepfrac;
                                      nseeds_cluster=nseeds_cluster, front_axis=front_axis, front_noise=front_noise)

            shap, _, _ = shapley_per_species(pool, grid, A, trues(grid.C), keep; bam=bam, mp=mp)
            cons = consumer_mask(pool)
            push!(ΔA_vec, mean(abs, getfield.(shap[cons], :dA)))
            push!(ΔB_vec, mean(abs, getfield.(shap[cons], :dB)))
            push!(ΔM_vec, mean(abs, getfield.(shap[cons], :dM)))
            push!(ΔF_vec, mean(getfield.(shap[cons], :dF)))
        end

        ΔA = mean(ΔA_vec); ΔB = mean(ΔB_vec); ΔM = mean(ΔM_vec); ΔF = mean(ΔF_vec)
        denom = max(ΔA + ΔB + ΔM, eps())
        out[hk] = (; ΔA, ΔB, ΔM, ΔF, S_A=ΔA/denom, S_B=ΔB/denom, S_M=ΔM/denom)
    end

    return out
end

# Ternary helpers
_ternary_xy(sA, sB, sM) = (sB + 0.5*sM, (√3/2)*sM)  # A on left base, B on right base, M to apex
function _ternary_frame!(ax; title="")
    xs = [0.0, 1.0, 0.5, 0.0]; ys = [0.0, 0.0, √3/2, 0.0]
    lines!(ax, xs, ys; color=:black)
    text!(ax, 0.02, -0.04, text="Abiotic", align=(:left,:bottom))
    text!(ax, 0.98, -0.04, text="Biotic", align=(:right,:bottom))
    text!(ax, 0.50,  √3/2 + 0.03, text="Movement", align=(:center,:bottom))
    ax.title = title
    hidespines!(ax, :r, :t)
end

"""
plot_geometry_effects_dual(; grid, S, basal_frac, A_level, B_level_on, M_level,
                            loss_fracs, loss_pick, seeds_pool, seeds_mask, sim_seed,
                            nseeds_cluster, front_axis, front_noise, τA, τocc)

Top row: ΔBSH vs loss for three geometries
  left  = B:none (climate-only, A×M)
  right = B=B_level_on (full BSH, A×B×M)

Bottom row: two elasticity triangles at loss_pick
  left  = elasticity shares with B:none
  right = elasticity shares with B=B_level_on
Points are colored by geometry; label next to each point shows S_B (biotic share).

Returns the figure and the two share Dicts.
"""
function plot_geometry_effects_dual(; grid::Grid, S::Int=150, basal_frac::Float64=0.45,
        A_level::Symbol=:intermediate, B_level_on::Symbol=:soft, M_level::Symbol=:on,
        loss_fracs=0.2:0.1:0.8, loss_pick::Float64=0.6,
        seeds_pool=1:6, seeds_mask=1:6, sim_seed::Int=1234,
        nseeds_cluster::Int=6, front_axis::Symbol=:x, front_noise::Float64=0.04,
        τA::Float64=0.5, τocc::Float64=0.35, T_frac_on::Float64=0.02)

    # --- ΔBSH vs loss (reuse your correct sweep)
    noB  = sweep_AM_vs_ABM(; grid,S,basal_frac, A_level, B_level=:none,       M_level,
                           hl_kind=:random,   loss_fracs,seeds_pool,seeds_mask,sim_seed,
                           nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)
    noB_C = sweep_AM_vs_ABM(; grid,S,basal_frac, A_level, B_level=:none,       M_level,
                            hl_kind=:clustered,loss_fracs,seeds_pool,seeds_mask,sim_seed,
                            nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)
    noB_F = sweep_AM_vs_ABM(; grid,S,basal_frac, A_level, B_level=:none,       M_level,
                            hl_kind=:front,    loss_fracs,seeds_pool,seeds_mask,sim_seed,
                            nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)

    yesB  = sweep_AM_vs_ABM(; grid,S,basal_frac, A_level, B_level=B_level_on, M_level,
                            hl_kind=:random,   loss_fracs,seeds_pool,seeds_mask,sim_seed,
                            nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)
    yesB_C = sweep_AM_vs_ABM(; grid,S,basal_frac, A_level, B_level=B_level_on, M_level,
                             hl_kind=:clustered,loss_fracs,seeds_pool,seeds_mask,sim_seed,
                             nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)
    yesB_F = sweep_AM_vs_ABM(; grid,S,basal_frac, A_level, B_level=B_level_on, M_level,
                             hl_kind=:front,    loss_fracs,seeds_pool,seeds_mask,sim_seed,
                             nseeds_cluster,front_axis,front_noise,τA,τocc,T_frac_on)

    # --- Elasticity shares at loss_pick (B off vs on)
    shares_off = geometry_elasticity_shares_all_geoms(; grid,S,basal_frac,A_level,
                       B_level=:none, M_level, loss_pick,
                       seeds_pool,seeds_mask,sim_seed,
                       nseeds_cluster,front_axis,front_noise,τA,τocc)

    shares_on  = geometry_elasticity_shares_all_geoms(; grid,S,basal_frac,A_level,
                       B_level=B_level_on, M_level, loss_pick,
                       seeds_pool,seeds_mask,sim_seed,
                       nseeds_cluster,front_axis,front_noise,τA,τocc)

    # --- Figure
    fig = Figure(; size=(1400, 900))
    Label(fig[0,:],
          @sprintf("Geometry effects without vs with biotic — A=%s, B(on)=%s, M=%s    (loss_pick=%.2f)",
                   String(A_level), String(B_level_on), String(M_level), loss_pick);
          fontsize=18, padding=(0,0,8,0))

    # ΔBSH curves (climate-only)
    ax1 = Axis(fig[1,1], title="ΔBSH vs loss — Biotic OFF (A×M)",
               xlabel="Area lost (fraction)", ylabel="Mean BSH (consumers)")
    lines!(ax1, noB.x, noB.mean_AM;  color=_gcolor[:random],   label="Random")
    lines!(ax1, noB_C.x, noB_C.mean_AM; color=_gcolor[:clustered], label="Clustered")
    lines!(ax1, noB_F.x, noB_F.mean_AM; color=_gcolor[:front],     label="Front-like")
    axislegend(ax1; position=:lb, framevisible=false, labelsize=10)

    # ΔBSH curves (full)
    ax2 = Axis(fig[1,2], title=@sprintf("ΔBSH vs loss — Biotic ON (A×B×M; B=%s)", String(B_level_on)),
               xlabel="Area lost (fraction)", ylabel="Mean BSH (consumers)")
    lines!(ax2, yesB.x, yesB.mean_ABM;    color=_gcolor[:random],   label="Random")
    lines!(ax2, yesB_C.x, yesB_C.mean_ABM;color=_gcolor[:clustered],label="Clustered")
    lines!(ax2, yesB_F.x, yesB_F.mean_ABM;color=_gcolor[:front],    label="Front-like")
    axislegend(ax2; position=:lb, framevisible=false, labelsize=10)

    # Elasticity triangles
    ax3 = Axis(fig[2,1]); _ternary_frame!(ax3; title="Elasticity shares at loss_pick — Biotic OFF")
    for hk in (:random,:clustered,:front)
        s = shares_off[hk]
        x,y = _ternary_xy(s.S_A, s.S_B, s.S_M)
        scatter!(ax3, [x],[y]; color=_gcolor[hk], markersize=12, label=_glabel[hk])
        text!(ax3, x, y; text=@sprintf(" B=%.2f ", s.S_B), align=(:left,:bottom), fontsize=10, color=_gcolor[hk])
    end
    axislegend(ax3; position=:rb, framevisible=false)

    ax4 = Axis(fig[2,2]); _ternary_frame!(ax4; title="Elasticity shares at loss_pick — Biotic ON")
    for hk in (:random,:clustered,:front)
        s = shares_on[hk]
        x,y = _ternary_xy(s.S_A, s.S_B, s.S_M)
        scatter!(ax4, [x],[y]; color=_gcolor[hk], markersize=12, label=_glabel[hk])
        text!(ax4, x, y; text=@sprintf(" B=%.2f ", s.S_B), align=(:left,:bottom), fontsize=10, color=_gcolor[hk])
    end
    axislegend(ax4; position=:rb, framevisible=false)

    display(fig)
    return (; fig, shares_off, shares_on)
end

grid = make_grid(60,60; seed=42)

res = plot_geometry_effects_dual(; grid,
        S=150, basal_frac=0.45,
        A_level=:divergent, B_level_on=:strong, M_level=:on,
        loss_fracs=0.2:0.1:0.8, loss_pick=0.6,
        seeds_pool=1:6, seeds_mask=1:6, sim_seed=1234,
        nseeds_cluster=6, front_axis=:x, front_noise=0.04,
        τA=0.5, τocc=0.35)
