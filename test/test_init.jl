using AlignedGP
using Test

@testset "init" begin
    symboltable = SymbolTable(2, [sqrt], [+])
    tree = init(symboltable, 10)
    evals = evaluate(tree, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    @test evals isa Vector{Float64}
end


@testset "test length distribution" begin
    symboltable = SymbolTable(5, [sqrt], [+])
    for sz = 1:15
        tree = init(symboltable, sz)
        @test length(tree) == sz
    end
end