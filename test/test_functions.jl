using Test
using AlignedGP
using Intervals

const CI = Interval{Float64,Closed,Closed}

Base.isapprox(a::CI, b::CI; kwargs...) =
    isapprox(first(a), first(b); kwargs...) && isapprox(last(a), last(b); kwargs...)

@testset "evaluate" begin
    @test evaluate(+, 3.0, 2.0) == 5.0
    @test evaluate(-, 3.0, 2.0) == 1.0
    @test evaluate(*, 3.0, 2.0) == 6.0
    @test evaluate(/, 6.0, 2.0) == 3.0

    @test evaluate(+, [1.0, 2.0], [3.0, 4.0]) == [4.0, 6.0]
    @test evaluate(-, [3.0, 4.0], [1.0, 2.0]) == [2.0, 2.0]
    @test evaluate(*, [2.0, 3.0], [4.0, 5.0]) == [8.0, 15.0]
    @test evaluate(/, [6.0, 9.0], [2.0, 3.0]) == [3.0, 3.0]
end

@testset "leftinverse +" begin
    @test leftinverse(+, CI(3.0, 5.0), 2.0)  == CI(1.0, 3.0)
    @test leftinverse(+, CI(3.0, 5.0), -1.0) == CI(4.0, 6.0)
    @test leftinverse(+, CI(0.0, 0.0), 0.0)  == CI(0.0, 0.0)

    ts = [CI(3.0, 5.0), CI(4.0, 6.0)]
    @test leftinverse(+, ts, [1.0, 2.0]) == [CI(2.0, 4.0), CI(2.0, 4.0)]
    @test leftinverse(+, ts, 1.0)        == [CI(2.0, 4.0), CI(3.0, 5.0)]
end

@testset "rightinverse +" begin
    @test rightinverse(+, CI(3.0, 5.0), 1.0)  == CI(2.0, 4.0)
    @test rightinverse(+, CI(3.0, 5.0), -1.0) == CI(4.0, 6.0)

    ts = [CI(3.0, 5.0), CI(4.0, 6.0)]
    @test rightinverse(+, ts, [1.0, 2.0]) == [CI(2.0, 4.0), CI(2.0, 4.0)]
end

@testset "leftinverse -" begin
    # x - y = t  →  x = t + y
    @test leftinverse(-, CI(3.0, 5.0), 2.0)  == CI(5.0, 7.0)
    @test leftinverse(-, CI(3.0, 5.0), -2.0) == CI(1.0, 3.0)

    ts = [CI(1.0, 3.0), CI(2.0, 4.0)]
    @test leftinverse(-, ts, [1.0, 2.0]) == [CI(2.0, 4.0), CI(4.0, 6.0)]
end

@testset "rightinverse -" begin
    # x - y = t  →  y = x - t  (bounds flip)
    @test rightinverse(-, CI(3.0, 5.0), 7.0) == CI(2.0, 4.0)
    @test rightinverse(-, CI(1.0, 4.0), 0.0) == CI(-4.0, -1.0)

    ts = [CI(1.0, 3.0), CI(2.0, 4.0)]
    @test rightinverse(-, ts, [5.0, 6.0]) == [CI(2.0, 4.0), CI(2.0, 4.0)]
end

@testset "leftinverse *" begin
    # x * y = t  →  x = t / y  (bounds flip when y < 0)
    @test leftinverse(*, CI(2.0, 6.0), 2.0)   == CI(1.0, 3.0)
    @test leftinverse(*, CI(2.0, 6.0), -2.0)  == CI(-3.0, -1.0)
    @test leftinverse(*, CI(-6.0, -2.0), -2.0) == CI(1.0, 3.0)

    ts = [CI(2.0, 4.0), CI(6.0, 9.0)]
    @test leftinverse(*, ts, [2.0, 3.0]) == [CI(1.0, 2.0), CI(2.0, 3.0)]
end

@testset "rightinverse *" begin
    @test rightinverse(*, CI(2.0, 6.0), 2.0)  == CI(1.0, 3.0)
    @test rightinverse(*, CI(2.0, 6.0), -2.0) == CI(-3.0, -1.0)

    ts = [CI(2.0, 4.0), CI(6.0, 9.0)]
    @test rightinverse(*, ts, [2.0, 3.0]) == [CI(1.0, 2.0), CI(2.0, 3.0)]
end

@testset "leftinverse /" begin
    # x / y = t  →  x = t * y  (bounds flip when y < 0)
    @test leftinverse(/, CI(1.0, 3.0), 2.0)   == CI(2.0, 6.0)
    @test leftinverse(/, CI(1.0, 3.0), -2.0)  == CI(-6.0, -2.0)
    @test leftinverse(/, CI(-3.0, -1.0), -2.0) == CI(2.0, 6.0)

    ts = [CI(1.0, 2.0), CI(2.0, 3.0)]
    @test leftinverse(/, ts, [2.0, 3.0]) == [CI(2.0, 4.0), CI(6.0, 9.0)]
end

@testset "rightinverse /" begin
    # x / y = t  →  y = x / t  (x/[l,u] → [x/u, x/l] when x,l,u > 0)
    @test rightinverse(/, CI(1.0, 3.0), 6.0)  == CI(2.0, 6.0)
    @test rightinverse(/, CI(1.0, 2.0), -4.0) == CI(-4.0, -2.0)

    ts = [CI(1.0, 2.0), CI(2.0, 3.0)]
    @test rightinverse(/, ts, [4.0, 6.0]) == [CI(2.0, 4.0), CI(2.0, 3.0)]
end

@testset "zero divisor propagates Inf (not a crash)" begin
    r = leftinverse(*, CI(1.0, 3.0), 0.0)
    @test isinf(first(r)) || isinf(last(r)) || isnan(first(r)) || isnan(last(r))
end

@testset "evaluate sqrt" begin
    @test evaluate(sqrt, 4.0) == 2.0
    @test evaluate(sqrt, [4.0, 9.0, 16.0]) == [2.0, 3.0, 4.0]
end

@testset "inverse sqrt" begin
    # normal case: [l,u] with 0 ≤ l ≤ u
    @test inverse(sqrt, CI(1.0, 3.0)) == CI(1.0, 9.0)
    @test inverse(sqrt, CI(0.0, 2.0)) == CI(0.0, 4.0)

    # l < 0 ≤ u: clamp l to 0
    @test inverse(sqrt, CI(-1.0, 3.0)) == CI(0.0, 9.0)

    # u < 0: unreachable target → invalid sentinel CI(Inf,Inf)
    @test inverse(sqrt, CI(-3.0, -1.0)) == CI(Inf, Inf)

    # degenerate point interval
    @test inverse(sqrt, CI(2.0, 2.0)) == CI(4.0, 4.0)

    # vector broadcasting
    ts = [CI(0.0, 2.0), CI(1.0, 3.0)]
    @test inverse(sqrt, ts) == [CI(0.0, 4.0), CI(1.0, 9.0)]
end

@testset "evaluate exp" begin
    @test evaluate(exp, 0.0) ≈ 1.0
    @test evaluate(exp, [0.0, 1.0]) ≈ [1.0, exp(1.0)]
end

@testset "inverse exp" begin
    # normal: 0 < l ≤ u
    @test inverse(exp, CI(1.0, exp(3.0))) ≈ CI(0.0, 3.0)

    # l ≤ 0 < u: lower bound becomes -Inf
    r = inverse(exp, CI(-1.0, exp(2.0)))
    @test first(r) == -Inf
    @test last(r) ≈ 2.0

    # l = 0: lower bound is -Inf (log(0) = -Inf)
    r2 = inverse(exp, CI(0.0, exp(1.0)))
    @test first(r2) == -Inf
    @test last(r2) ≈ 1.0

    # u ≤ 0: unreachable → invalid sentinel
    @test inverse(exp, CI(-2.0, -1.0)) == CI(Inf, Inf)
    @test inverse(exp, CI(-1.0,  0.0)) == CI(Inf, Inf)

    # vector broadcasting
    ts = [CI(1.0, exp(1.0)), CI(exp(1.0), exp(2.0))]
    rs = inverse(exp, ts)
    @test first(rs[1]) ≈ 0.0 && last(rs[1]) ≈ 1.0
    @test first(rs[2]) ≈ 1.0 && last(rs[2]) ≈ 2.0
end

@testset "evaluate log" begin
    @test evaluate(log, 1.0) ≈ 0.0
    @test evaluate(log, [1.0, exp(1.0)]) ≈ [0.0, 1.0]
end

@testset "inverse log" begin
    # normal case
    @test inverse(log, CI(0.0, 1.0)) ≈ CI(1.0, exp(1.0))
    @test inverse(log, CI(-1.0, 1.0)) ≈ CI(exp(-1.0), exp(1.0))

    # l = -Inf: lower bound becomes exp(-Inf) = 0
    r = inverse(log, CI(-Inf, 1.0))
    @test first(r) == 0.0
    @test last(r) ≈ exp(1.0)

    # always valid — no sentinel for log inverse
    @test inverse(log, CI(-100.0, 100.0)) ≈ CI(exp(-100.0), exp(100.0))

    # vector broadcasting
    ts = [CI(0.0, 1.0), CI(1.0, 2.0)]
    rs = inverse(log, ts)
    @test rs[1] ≈ CI(1.0, exp(1.0))
    @test rs[2] ≈ CI(exp(1.0), exp(2.0))
end
