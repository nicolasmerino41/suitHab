Fs = 0.2:0.02:0.8
geom = :front
best = (f=NaN, KS=-Inf)
for f in Fs
    rAM, rBAM = BSH.per_species_relative_loss(rng, pool, grid, pars; fstar=f, geometry=geom, seed_A=1)
    m = Metrics.dist_metrics(rAM, rBAM; tail_cut=0.8)
    if m.KS > best.KS; best = (f=f, KS=m.KS); end
end
@info "Best f* by KS (geometry=$(geom)) => f*=$(round(best.f,digits=2)), KS=$(round(best.KS,digits=3))"
