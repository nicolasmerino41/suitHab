using Logging
using Serialization

const S_VALUES = collect(150:50:400)
const S_BASE   = 250
const C_BASE   = 0.055

function run_sensitivity_richness_khat(; tag::Int=41)
    # baseline mean degree (directed): k̂0 = C0 * S0
    Khat0 = C_BASE * S_BASE
    @info "k̂-constant richness sweep" S_BASE C_BASE Khat0

    # Keep the same r-grid, only overwrite C each time
    r_values = [cr[2] for cr in CELLS_CR]
    CELLS_CR_ORIG = copy(CELLS_CR)

    outfiles = String[]
    for (si, Sp) in enumerate(S_VALUES)
        outfile = joinpath(OUTDIR, "sens_S_khat_S$(Sp).jls")
        push!(outfiles, outfile)

        if isfile(outfile)
            @info "Skipping (already computed)" Sp si=si n=length(S_VALUES)
            continue
        end

        Cnew = clamp(Khat0 / Sp, 1e-6, 0.5)
        @info "Running S point" Sp Cnew si=si n=length(S_VALUES)

        # mutate CELLS_CR in-place (length unchanged => compatible with run_sensitivity_threaded)
        @inbounds for i in eachindex(CELLS_CR)
            CELLS_CR[i] = (Cnew, r_values[i])
        end

        # run just this one S as a 1-point sensitivity
        res_one = run_sensitivity_threaded([Sp], x -> begin
                grid = build_grid(GRID_BASE, GRID_BASE)
                Emin = Emin_patch_base
                suit = SUIT_THRESH_BASE
                return (grid, Emin, suit, Int(x))
            end; tag=tag, label_x="S")

        serialize(outfile, res_one)
        @info "Saved" Sp outfile
    end

    # restore original CELLS_CR contents
    CELLS_CR .= CELLS_CR_ORIG

    # ---- merge per-S results into a full results dict over all S_VALUES
    first_idx = findfirst(isfile, outfiles)
    first_idx === nothing && error("No results files found; nothing to merge.")
    first_res = deserialize(outfiles[first_idx])

    results_merged = Dict{Tuple{Symbol,Symbol,Int,Symbol,Int}, Vector{Float64}}()
    for (k, v) in first_res
        results_merged[k] = fill(NaN, length(S_VALUES))
    end

    for (si, Sp) in enumerate(S_VALUES)
        outfile = outfiles[si]
        isfile(outfile) || continue
        res_one = deserialize(outfile)
        for (k, v) in res_one
            results_merged[k][si] = v[1]
        end
    end

    return results_merged
end

res_S = run_sensitivity_richness_khat()

# Plot later with a Windows-safe expname (ASCII, no slashes/colons)
plot_sensitivity_lines(S_VALUES, res_S;
    expname="Richness sweep (khat fixed; C=khat_over_S)",
    xlabel="Regional richness S"
)
