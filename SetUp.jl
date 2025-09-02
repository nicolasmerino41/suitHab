@time begin
    PC = ENV["USERNAME"]
    using Pkg
    using Random, CairoMakie, StatsBase, Statistics
    Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\SuitHab"))
    cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\SuitHab"))
    import Base.Threads: @threads
end