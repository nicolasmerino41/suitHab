# climate-normalized prey hotspot to downweight climate alignment
function prey_hotspot_residual_mask(grid::Grid, pool::SpeciesPool;
                                    τ::Float64, keep_frac::Float64,
                                    power::Float64=2.5, seed::Int=0)
    Z = climate_pass(pool, grid; τ=τ)
    cons = .!pool.basal
    ncons = max(1, count(cons))
    cons_clim_score = vec(sum(@view Z[cons, :]; dims=1)) ./ ncons  # mean consumer climate presence

    _, prey_score = _consumer_and_prey_scores(pool, grid; τ=τ)
    eps = 1e-6
    residual = prey_score ./ (cons_clim_score .+ eps)               # prey “per unit climate”
    return _top_mask_by_score(residual, keep_frac; power=power, seed=seed)
end

function consumer_hotspot_residual_mask(grid::Grid, pool::SpeciesPool;
                                    τ::Float64, keep_frac::Float64,
                                    power::Float64=2.5, seed::Int=0)
    Z = climate_pass(pool, grid; τ=τ)
    cons = .!pool.basal
    ncons = max(1, count(cons))
    cons_clim_score = vec(sum(@view Z[cons, :]; dims=1)) ./ ncons  # mean consumer climate presence

    _, prey_score = _consumer_and_prey_scores(pool, grid; τ=τ)
    eps = 1e-6
    residual = cons_clim_score ./ (prey_score .+ eps)               # consumer “per unit prey”
    return _top_mask_by_score(residual, keep_frac; power=power, seed=seed)
end
    