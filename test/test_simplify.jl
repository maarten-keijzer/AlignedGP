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
    # Regression: an integer-valued offset above 2^53 used to throw InexactError
    # (Int overflow) in _num. It must render as a compact Float, not a 106-digit
    # exact integer whose trailing digits are meaningless binary artifacts.
    s = simplify(Var(9, -2.9016275607237895e105))
    @test occursin("e+105", s)                 # compact scientific form
    @test !occursin("29016275607237895", s)    # not the expanded 106-digit int
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

@testset "_simplest_rational — continued-fraction search" begin
    sr(lo, hi; cap=10_000) = AlignedGP._simplest_rational(lo, hi; den_cap=cap)

    # Basic hits: simplest (min-denominator) rational in the interval.
    @test sr(0.3, 0.4)         == (1, 3)
    @test sr(0.49, 0.51)       == (1, 2)
    @test sr(1.9, 2.4)         == (2, 1)   # integer wins (denominator 1)
    @test sr(0.33333, 0.33334) == (1, 3)

    # Sign and zero handling.
    @test sr(-0.4, -0.3) == (-1, 3)
    @test sr(-0.1, 0.1)  == (0, 1)
    @test sr(3.2, Inf)   == (4, 1)         # half-infinite interval

    # Denominator cap: too-tight interval has no rational under the cap.
    @test sr(0.333331, 0.333332; cap=100) === nothing

    # Degenerate point intervals (lo == hi) used to recurse on their own
    # continued fraction forever and overflow the stack. They must terminate:
    # returning the exact rational when one fits under the cap, else `nothing`.
    @test sr(0.5, 0.5) == (1, 2)
    @test sr(0.1, 0.1) == (1, 10)
    for x in (Float64(π), exp(1.0), sqrt(2.0) - 1, 6.9652149972729935)
        @test sr(x, x) === nothing         # irrational-valued: no small rational
    end

    # Regression: a wide interval around a huge integer used to walk the
    # Stern-Brocot tree one mediant at a time (~1e105 steps) and hang, because
    # the denominator stayed 1 so the den_cap guard never fired. It must now
    # return promptly; the numerator is too large to keep exact -> `nothing`.
    local r
    t = @elapsed r = sr(2.6689566025148944e105, 3.348824162893223e105)
    @test r === nothing
    @test t < 1.0
end

@testset "simplify — region reveals constant multiples" begin
    @test simplify(region_const(1.5707, 1.5709, π/2)) == "pi/2"
    @test simplify(region_const(3.1415, 3.1417, π))   == "pi"
    # reciprocal form: 1/pi is not a rational *multiple* of pi
    @test simplify(region_const(0.3183, 0.3184, 1/π)) == "1/pi"
end

@testset "simplify — degenerate region does not overflow" begin
    # A point region (lo == hi) at an irrational value drove _cf_simplest into
    # unbounded recursion (StackOverflowError). It must resolve and fall back to
    # the raw value when no simple rational fits the region.
    @test simplify(region_const(Float64(π), Float64(π), Float64(π))) isa String
    @test simplify(region_const(0.5, 0.5, 0.5)) == "1/2"
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
