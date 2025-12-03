module RealizedNetwork

export realized_metaweb

function realized_metaweb(metaweb, F)
    S = length(F)
    A = metaweb.A
    present = [any(F[i]) for i in 1:S]

    R = zeros(Int, S, S)
    for i in 1:S, j in 1:S
        if A[i,j] == 1 && present[i] && present[j]
            R[i,j] = 1
        end
    end
    R
end

end
