using Test
using AlignedGP
using AlignedGP.ReverseIntervals: intervaltype, IntervalType

# A Constant whose value lives in an equivalence region [lo, hi].
region_const(lo, hi, value) = Constant(AddedValue([intervaltype(lo, hi)], value))

@testset "simplify — leaf symbol" begin
    @test simplify(Var(1)) == "x1"
end

@testset "simplify — constant" begin
    @test simplify(Constant(3.0)) == "3"        # integer-valued float prints cleanly
    @test simplify(Constant(0.5)) == "0.500000000000000"
end

@testset "simplify — additive offset" begin
    @test simplify(Var(1, 2.0)) == "x1 + 2"
    @test simplify(Var(1))      == "x1"          # zero offset adds nothing
end

@testset "simplify — binary operators" begin
    @test simplify(BinaryNode(+, Var(1), Var(2))) == "x1 + x2"
    @test simplify(BinaryNode(-, Var(1), Var(2))) == "x1 - x2"
    @test simplify(BinaryNode(*, Var(1), Var(2))) == "x1*x2"
    @test simplify(BinaryNode(/, Var(1), Var(2))) == "x1/x2"
end

@testset "simplify — actually simplifies & nests" begin
    @test simplify(BinaryNode(+, Var(1), Var(1))) == "2*x1"
    @test simplify(BinaryNode(-, Var(1), Var(1))) == "0"
    @test simplify(UnaryNode(log, UnaryNode(exp, Var(1)))) == "x1"   # real assumption
    # offset on an internal node
    @test simplify(BinaryNode(*, Var(1), Var(2), AddedValue(1.0))) == "x1*x2 + 1"
end

@testset "simplify — region picks simplest rational" begin
    @test simplify(region_const(0.49, 0.51, 0.5)) == "1/2"
    @test simplify(region_const(1.9, 2.4, 2.0))   == "2"     # integer falls out
    # no region -> honest float, no rationalization
    @test simplify(Constant(0.5)) == "0.500000000000000"
end

@testset "simplify — region reveals constant multiples" begin
    @test simplify(region_const(1.5707, 1.5709, π/2)) == "pi/2"
    @test simplify(region_const(3.1415, 3.1417, π))   == "pi"
    # reciprocal form: 1/pi is not a rational *multiple* of pi
    @test simplify(region_const(0.3183, 0.3184, 1/π)) == "1/pi"
end

@testset "simplify — Tree delegates to root" begin
    tree = Tree(BinaryNode(+, Var(1), Var(1)), BitVector())
    @test simplify(tree) == "2*x1"
    @test simplify(tree, form=:str) == "2*x1"
end

@testset "simplify — unsupported function errors" begin
    @test_throws Exception simplify(UnaryNode(abs, Var(1)))
    @test_throws Exception simplify(BinaryNode(^, Var(1), Var(2)))
end

@testset "simplify — output forms" begin
    node = UnaryNode(sqrt, Var(1))
    @test simplify(node) == "sqrt(x1)"                      # default :str
    @test simplify(node, form=:latex) == "\\sqrt{x_{1}}"
    pretty = simplify(node, form=:pretty)
    @test pretty isa String && occursin("x", pretty)        # 2D unicode form
    @test_throws Exception simplify(node, form=:bogus)
end

@testset "simplify — sign assumptions from inputs" begin
    node = UnaryNode(sqrt, BinaryNode(*, Var(1), Var(1)))   # sqrt(x1^2)
    @test simplify(node) == "Abs(x1)"                        # bare: real only
    @test simplify(node, [[1.0, 2.0, 3.0]]) == "x1"          # all > 0  -> positive
    @test simplify(node, [[0.0, 1.0, 2.0]]) == "x1"          # all >= 0 -> nonnegative
    @test simplify(node, [[-1.0, 1.0, 2.0]]) == "Abs(x1)"    # has < 0  -> real
end

@testset "simplify — unary functions" begin
    @test simplify(UnaryNode(-, Var(1)))    == "-x1"
    @test simplify(UnaryNode(inv, Var(1)))  == "1/x1"
    @test simplify(UnaryNode(sqrt, Var(1))) == "sqrt(x1)"
    @test simplify(UnaryNode(log, Var(1)))  == "log(x1)"
    @test simplify(UnaryNode(exp, Var(1)))  == "exp(x1)"
    @test simplify(UnaryNode(sin, Var(1)))  == "sin(x1)"
    @test simplify(UnaryNode(cos, Var(1)))  == "cos(x1)"
end
