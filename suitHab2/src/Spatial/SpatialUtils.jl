module SpatialUtils

export norm01, gaussian_blur, label_components

using LinearAlgebra, ImageFiltering

norm01(A) = (A .- minimum(A)) ./ max(eps(), maximum(A) - minimum(A))

function gaussian_blur(A, σ)
    kernel = Kernel.gaussian(σ)
    imfilter(A, kernel)
end

function label_components(mask::BitMatrix)
    nx, ny = size(mask)
    labels = fill(0, nx, ny)
    comp = 0

    for i in 1:nx, j in 1:ny
        if mask[i,j] && labels[i,j] == 0
            comp += 1
            stack = [(i,j)]
            labels[i,j] = comp

            while !isempty(stack)
                (x,y) = pop!(stack)
                for (dx,dy) in ((1,0),(-1,0),(0,1),(0,-1))
                    xx = x + dx
                    yy = y + dy
                    if 1 ≤ xx ≤ nx && 1 ≤ yy ≤ ny
                        if mask[xx,yy] && labels[xx,yy] == 0
                            labels[xx,yy] = comp
                            push!(stack, (xx,yy))
                        end
                    end
                end
            end
        end
    end

    labels, comp
end

end
