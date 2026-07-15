using Test

@testset "Reverse Functions" begin 
    include("test_rev_functions.jl")
end

@testset "Intervals" begin 
    include("test_Intervals.jl")
end

using AlignedGP

@testset "Interval Alignment" begin
    include("test_interval_alignment.jl")
end
@testset "Functions" begin
    include("test_functions.jl")
end
@testset "Trees" begin
    include("test_trees.jl")
end
@testset "Alignment Surrogate" begin
    include("test_alignment_surrogate.jl")
end
@testset "Initialization" begin
    include("test_init.jl")
end
@testset "Crossover" begin
    include("test_crossover.jl")
end
@testset "Linear Scaling" begin
    include("test_linear_scaling.jl")
end
