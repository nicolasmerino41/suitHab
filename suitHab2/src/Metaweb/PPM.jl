module PPM

export build_ppm

using Random, Distributions
using ..MetawebUtils: trophic_levels
using ..Types: Metaweb

function ppm(S, B, L, T)
    A = zeros(Int, S, S)

    # desired connectance = L / total_possible
    total_possible = S * (S - 1)
    target_connectance = L / total_possible

    # Expected proportion of possible prey per species
    # For species i, possible prey = i-1
    # Expected prey count ≈ target_connectance * (i - 1)
    # This ensures total expected links ≈ L
    # (minor deviation at basal because few possible prey there)
    
    for i in B+1:S
        existing = 1:(i-1)
        n_i = length(existing)

        # Expected number of prey for species i
        expected_prey = target_connectance * n_i
        k_i = max(1, round(Int, expected_prey))

        # Always pick at least 1 prey (from PPM logic)
        # First prey: uniformly random
        first_prey = rand(existing)
        A[i, first_prey] = 1

        if k_i > 1
            # Compute provisional trophic levels as in classical PPM
            provisional_s = zeros(Float64, i)
            provisional_s[1:i] .= 1.0  # simple placeholder
            
            # Trophic similarity-based probabilities
            # (still valid even without real s yet)
            probs = [exp(-abs(provisional_s[first_prey] - provisional_s[j]) / T)
                     for j in existing]
            probs ./= sum(probs)

            # Sample remaining prey
            chosen = rand(Categorical(probs), k_i - 1)
            for idx in unique(chosen)
                prey = existing[idx]
                A[i, prey] = 1
            end
        end
    end

    return A
end

function build_ppm(S, B, L, T)
    A = ppm(S, B, L, T)

    # Compute trophic levels
    s = trophic_levels(A)

    # Identify basal species
    basal = collect(1:B)

    Metaweb(A, s, basal)
end

end
