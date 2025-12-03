module PPM

export build_ppm

using Random, Distributions
using ..MetawebUtils: trophic_levels
using ..Types: Metaweb

function ppm(S, B, L, T)
    A = zeros(Int, S, S)

    β = (S^2 - B^2) / (2L - 1)
    beta_dist = Beta(β, β)

    prey_count = zeros(Int, S)

    current = B

    while current < S
        current += 1
        i = current

        existing = 1:(i-1)
        n_i = length(existing)

        # 1) First prey
        j = rand(existing)
        A[i,j] = 1
        prey_count[i] += 1

        # provisional TLs:
        s_hat = 1 .+ prey_count[1:i]

        # 2) Expected prey count
        k_exp = rand(beta_dist) * n_i
        k_total = max(1, round(Int, k_exp))
        k_extra = k_total - 1

        if k_extra > 0
            # probabilities based on provisional TLs
            probs = [exp(-abs(s_hat[j] - s_hat[ℓ]) / T) for ℓ in existing]
            probs ./= sum(probs)

            chosen = rand(Distributions.Categorical(probs), k_extra)
            for idx in unique(chosen)
                prey = existing[idx]
                A[i, prey] = 1
                prey_count[i] += 1
            end
        end
    end

    # final, correct TL computation
    s = trophic_levels(A)

    return A, s
end

function build_ppm(S, B, L, T)
    A, s = ppm(S, B, L, T)

    # # Compute trophic levels
    # s = trophic_levels(A)
    
    # Identify basal species
    basal = collect(1:B)

    Metaweb(A, s, basal)
end

end
