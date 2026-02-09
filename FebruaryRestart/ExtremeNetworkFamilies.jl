#!/usr/bin/env julia
using Random
using CairoMakie

# ============================================================
# EXTREME toy parameters (for visual clarity)
# ============================================================

const S = 300                 # divisible by modules
const BASAL_FRAC = 0.20
const C = 0.18              # denser for clearer patterns
const RNG_SEED = 22

# Exaggerated structure knobs
const N_MODULES = 3
const MODULAR_IN_BIAS = 25.0        # HUGE modular bias
const HEAVYTAIL_GAMMA = 1.6         # very heavy tail
const CASCADE_LAMBDA = 8.0          # strong hierarchy

rng = MersenneTwister(RNG_SEED)

# ============================================================
# Helpers
# ============================================================

function consumers_and_basal()
    nb = round(Int, BASAL_FRAC * S)
    basal = falses(S)
    basal[1:nb] .= true
    consumers = findall(!, basal)
    return basal, consumers
end

function prey_to_adj(prey)
    A = zeros(Float64, S, S)
    for i in 1:S, j in prey[i]
        A[i,j] = 1.0
    end
    return A
end

function outdegrees(prey, basal)
    d = Int[]
    for i in 1:S
        basal[i] && continue
        push!(d, length(prey[i]))
    end
    return d
end

# ============================================================
# Network builders
# ============================================================

function build_random(rng, basal)
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal)
    Ltarget = round(Int, C * S^2)

    for i in consumers
        j = rand(rng, 1:S); j==i && (j=mod1(j+1,S))
        push!(prey[i], j)
    end

    L = length(consumers)
    while L < Ltarget
        i = consumers[rand(rng, 1:end)]
        j = rand(rng, 1:S)
        (j!=i && j∉prey[i]) || continue
        push!(prey[i], j); L+=1
    end
    return prey
end

function build_modular(rng, basal)
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal)

    MODULE = [1 + (i-1) ÷ (S÷N_MODULES) for i in 1:S]
    Ltarget = round(Int, C * S^2)

    function sample_prey(i)
        mi = MODULE[i]
        inmod = [j for j in 1:S if MODULE[j]==mi && j!=i]
        outmod = [j for j in 1:S if MODULE[j]!=mi]
        rand(rng) < MODULAR_IN_BIAS/(MODULAR_IN_BIAS+1) ? rand(rng,inmod) : rand(rng,outmod)
    end

    for i in consumers
        push!(prey[i], sample_prey(i))
    end

    L = length(consumers)
    while L < Ltarget
        i = consumers[rand(rng,1:end)]
        j = sample_prey(i)
        j∉prey[i] || continue
        push!(prey[i], j); L+=1
    end
    return prey
end

function build_heavytail(rng, basal)
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal)
    nC = length(consumers)
    Ltarget = round(Int, C * S^2)

    w = [rand(rng)^(-1/(HEAVYTAIL_GAMMA-1)) for _ in 1:nC]
    w ./= sum(w)

    deg = ones(Int,nC)
    for _ in 1:(Ltarget-nC)
        k = findfirst(cumsum(w) .>= rand(rng))
        deg[k]+=1
    end

    for (kk,i) in enumerate(consumers)
        cand = shuffle(rng,[j for j in 1:S if j!=i])
        prey[i] = cand[1:min(deg[kk],length(cand))]
    end
    return prey
end

function build_cascade(rng, basal)
    prey = [Int[] for _ in 1:S]
    consumers = findall(!, basal)
    ranks = rand(rng,S)
    Ltarget = round(Int,C*S^2)

    function sample_lower(i)
        lower=[j for j in 1:S if ranks[j]<ranks[i] && j!=i]
        isempty(lower) && return rand(rng,[j for j in 1:S if j!=i])
        w=exp.(-CASCADE_LAMBDA.*(ranks[i].-ranks[lower]))
        sw=sum(w); sw<=0 && return rand(rng,lower)
        w./=sw
        r=rand(rng); acc=0.0
        for k in 1:length(lower)
            acc+=w[k]; r<=acc && return lower[k]
        end
        return lower[end]
    end

    for i in consumers
        push!(prey[i], sample_lower(i))
    end

    L=length(consumers)
    while L<Ltarget
        i=consumers[rand(rng,1:end)]
        j=sample_lower(i)
        j∉prey[i] || continue
        push!(prey[i],j); L+=1
    end
    return prey,ranks
end

# ============================================================
# Build networks
# ============================================================

basal,_=consumers_and_basal()

nets=Dict(
    "Random"=>build_random(rng,basal),
    "Modular"=>build_modular(rng,basal),
    "Heavy-tail"=>build_heavytail(rng,basal),
)

cascade,ranks=build_cascade(rng,basal)
nets["Cascade"]=cascade

# ============================================================
# Plot
# ============================================================
begin
    f=Figure(size=(1200,900))
    Label(f[0,:],"Extreme network families (adjacency + out-degree)",fontsize=20)

    names=collect(keys(nets))

    for (i,name) in enumerate(names[3:4])
        prey=nets[name]
        A=prey_to_adj(prey)
        deg=outdegrees(prey,basal)

        # reorder cascade by rank for clean triangle
        if name=="Cascade"
            ord=sortperm(ranks,rev=true)
            A=A[ord,ord]
        end

        axM=Axis(f[1,i],title=name)
        heatmap!(axM,A)
        axM.xlabel="Prey"; axM.ylabel="Consumer"

        axD=Axis(f[2,i])
        hist!(axD,deg,bins=8)
        axD.xlabel="Out-degree"
    end

    save("network_families_extreme.png",f)
    display(f)
end
println("Saved: network_families_extreme.png")
