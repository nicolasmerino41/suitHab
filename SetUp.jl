@time begin
    PC = ENV["USERNAME"]
    using Pkg
    using Random, CairoMakie, StatsBase, Statistics
    Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio2"))
    cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio2"))
    meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")
    import Base.Threads: @threads
end