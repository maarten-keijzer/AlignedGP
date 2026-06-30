using AlignedGP
using Test

@testset "Interval Alignment" begin
    include("test_interval_alignment.jl")
end
@testset "Functions" begin
    include("test_functions.jl")
end
@testset "Trees" begin
    include("test_trees.jl")
end
