@time begin
    PC = ENV["USERNAME"]
    using Pkg
    using Random, CairoMakie, StatsBase, Statistics
    Pkg.activate(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\suitHab"))
    cd(joinpath("C:\\Users", PC, "OneDrive\\PhD\\GitHub\\suitHab"))
    import Base.Threads: @threads
end