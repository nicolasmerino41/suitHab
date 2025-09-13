"""
plot_audit_phase_rule_basic(audit; title)
- audit :: Vector{<:NamedTuple} returned by audit_phase_rule(...)
Draws:
  A) Confusion (observed worst vs predicted worst)
  B) Geometry spread bar (max-min ΔBSH) per combo, sorted
  C) Per-combo ΔBSH by geometry (lollipop)
"""
function plot_audit_phase_rule_basic(audit; title="Audit of phase rule (baseline predictors)")
    geoms = (:random, :clustered, :front)
    glabel = Dict(:random=>"Random", :clustered=>"Clustered", :front=>"Front-like")
    col = Dict(:random=>RGBf(0.25,0.45,0.85), :clustered=>RGBf(0.95,0.60,0.15), :front=>RGBf(0.20,0.65,0.35))

    # table helpers
    obs_worst = [a.obs_order[1] for a in audit]
    pred_worst = [a.pred_worst for a in audit]
    spread = [maximum([a.d_random,a.d_clustered,a.d_front]) - minimum([a.d_random,a.d_clustered,a.d_front]) for a in audit]
    names = [Symbol("A_"*String(a.A)*"_B_"*String(a.B)*"_M_"*String(a.M)) for a in audit]

    # confusion counts
    idx = Dict(:random=>1, :clustered=>2, :front=>3)
    M = zeros(Int, 3, 3)  # rows = observed, cols = predicted
    for i in eachindex(audit)
        M[idx[obs_worst[i]], idx[pred_worst[i]]] += 1
    end

    # figure
    fig = Figure(; size=(1400, 600))
    Label(fig[0,:], title; fontsize=18, padding=(0,0,8,0))

    # A) Confusion
    axA = Axis(fig[1,1], title="Observed worst vs Predicted worst", xlabel="Predicted", ylabel="Observed",
               xticks=(1:3, [glabel[g] for g in geoms]), yticks=(1:3, [glabel[g] for g in geoms]),
               xticklabelrotation=π/6)
    heatmap!(axA, M; colormap=:Blues)
    for r in 1:3, c in 1:3
        text!(axA, c, r, text=string(M[r,c]), align=(:center,:center), color=:black)
    end

    # B) Spread bars
    ord = sortperm(spread; rev=true)
    axB = Axis(fig[1,2], title="Geometry spread = max(ΔBSH) - min(ΔBSH)", ylabel="Spread (ΔBSH)", xticklabelrotation=π/6)
    barplot!(axB, 1:length(audit), spread[ord])
    axB.xticks = (1:length(audit), string.(names[ord]))
    # rotate!(axB.xticklabelrotation, π/2)

    # C) Per-combo ΔBSH by geometry (lollipop)
    axC = Axis(
        fig[2,1:2], 
        title="ΔBSH by geometry (each combo)", xlabel="Combo", ylabel="ΔBSH",
        xticklabelrotation=π/6
    )
    xs = 1:length(audit)
    # draw stems from min to max
    mins = [minimum((audit[i].d_random, audit[i].d_clustered, audit[i].d_front)) for i in 1:length(audit)]
    maxs = [maximum((audit[i].d_random, audit[i].d_clustered, audit[i].d_front)) for i in 1:length(audit)]
    lines!(axC, xs, mins; color=RGBAf(0,0,0,0.15))
    lines!(axC, xs, maxs; color=RGBAf(0,0,0,0.15))
    # place points
    r = [a.d_random for a in audit]; c = [a.d_clustered for a in audit]; f = [a.d_front for a in audit]
    scatter!(axC, xs, r; color=col[:random],    markersize=7)
    scatter!(axC, xs, c; color=col[:clustered], markersize=7)
    scatter!(axC, xs, f; color=col[:front],     markersize=7)
    axislegend(axC, [scatter!(Point2f[]; color=col[:random]),
                     scatter!(Point2f[]; color=col[:clustered]),
                     scatter!(Point2f[]; color=col[:front])],
               ["Random","Clustered","Front-like"]; position=:rb, framevisible=false)

    # x-ticks last (ordered by spread)
    axC.xticks = (xs, string.(names))

    display(fig)
    return fig
end

plot_audit_phase_rule_basic(A)