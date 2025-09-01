@time begin
    PC = ENV["USERNAME"]
    using Random, CairoMakie, StatsBase, Statistics
    using Pkg
    Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio2"))
    cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\simBio2"))
    meta_path = joinpath("C:\\Users", PC, "OneDrive\\PhD\\Metaweb Modelling")
    import Base.Threads: @threads
end