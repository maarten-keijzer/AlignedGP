using Test
using AlignedGP
using AlignedGP.ReverseIntervals
using IntervalArithmetic: issubset_interval, in_interval, bounds, inf, sup

# utility functions for easier testing
AlignedGP.inverse(op, iv::IntervalType, y::Float64) = inverse(op, IntervalVector([iv]), [y])[1]
AlignedGP.inverse(op, x::Float64, iv::IntervalType) = inverse(op, [x], IntervalVector([iv]))[1]
AlignedGP.inverse(op, iv::IntervalType) = inverse(op, IntervalVector([iv]))[1]

Base.isapprox(iv1::IntervalType, iv2::IntervalType) = iv1.lo ≈ iv2.lo && iv1.hi ≈ iv2.hi

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
    @test inverse(+, intervaltype(3.0, 5.0), 2.0)[1]  ≈ intervaltype(1.0, 3.0)
    @test inverse(+, intervaltype(3.0, 5.0), -1.0)[1] ≈ intervaltype(4.0, 6.0)
    
    # point intervals can work under addition
    @test inverse(+, intervaltype(0.0, 0.0), 1.0)[1] ≈ intervaltype(-1.0, -1.0)

    ts = IntervalVector([intervaltype(3.0, 5.0), intervaltype(4.0, 6.0)])
    @test all(inverse(+, ts, [1.0, 2.0]).intervals .≈ [intervaltype(2.0, 4.0), intervaltype(2.0, 4.0)])
    @test all(inverse(+, ts, [1.0, 1.0]).intervals .≈ [intervaltype(2.0, 4.0), intervaltype(3.0, 5.0)])
end

@testset "rightinverse +" begin
    @test inverse(+, 1.0, intervaltype(3.0, 5.0))[1]  ≈ intervaltype(2.0, 4.0)
    @test inverse(+, -1.0, intervaltype(3.0, 5.0))[1] ≈ intervaltype(4.0, 6.0)

    ts = IntervalVector([intervaltype(3.0, 5.0), intervaltype(4.0, 6.0)])
    @test all(inverse(+, [1.0, 2.0], ts).intervals .≈ [intervaltype(2.0, 4.0), intervaltype(2.0, 4.0)])
end

@testset "leftinverse -" begin
    @test inverse(-, intervaltype(3.0, 5.0), 2.0)[1]  ≈ intervaltype(5.0, 7.0)
    @test inverse(-, intervaltype(3.0, 5.0), -2.0)[1] ≈ intervaltype(1.0, 3.0)

    ts = IntervalVector([intervaltype(1.0, 3.0), intervaltype(2.0, 4.0)])
    @test all(inverse(-, ts, [1.0, 2.0]).intervals .≈ [intervaltype(2.0, 4.0), intervaltype(4.0, 6.0)])
end

@testset "rightinverse -" begin
    @test inverse(-, 7.0, intervaltype(3.0, 5.0))[1] ≈ intervaltype(2.0, 4.0)
    @test inverse(-, 0.0, intervaltype(1.0, 4.0))[1] ≈ intervaltype(-4.0, -1.0)
    @test inverse(-, 5.0, intervaltype(1.0, 3.0))[1] ≈ intervaltype(2.0, 4.0)

    ts = IntervalVector([intervaltype(1.0, 3.0), intervaltype(2.0, 4.0)])
    @test all(inverse(-, [5.0, 6.0], ts).intervals .≈ [intervaltype(2.0, 4.0), intervaltype(2.0, 4.0)])
end

@testset "leftinverse *" begin
    @test inverse(*, intervaltype(2.0, 6.0), 2.0)[1]   ≈ intervaltype(1.0, 3.0)
    @test inverse(*, intervaltype(2.0, 6.0), -2.0)[1]  ≈ intervaltype(-3.0, -1.0)
    @test inverse(*, intervaltype(-6.0, -2.0), -2.0)[1] ≈ intervaltype(1.0, 3.0)

    ts = IntervalVector([intervaltype(2.0, 4.0), intervaltype(6.0, 9.0)])
    @test all(inverse(*, ts, [2.0, 3.0]).intervals .≈ [intervaltype(1.0, 2.0), intervaltype(2.0, 3.0)])
end

@testset "rightinverse *" begin
    @test inverse(*, 2.0, intervaltype(2.0, 6.0))[1]  ≈ intervaltype(1.0, 3.0)
    @test inverse(*, -2.0, intervaltype(2.0, 6.0))[1] ≈ intervaltype(-3.0, -1.0)

    ts = IntervalVector([intervaltype(2.0, 4.0), intervaltype(6.0, 9.0)])
    @test all(inverse(*, [2.0, 3.0], ts).intervals .≈ [intervaltype(1.0, 2.0), intervaltype(2.0, 3.0)])
end


@testset "leftinverse /" begin
    @test inverse(/, intervaltype(1.0, 3.0), 2.0)[1]   ≈ intervaltype(2.0, 6.0)
    @test inverse(/, intervaltype(1.0, 3.0), -2.0)[1]  ≈ intervaltype(-6.0, -2.0)
    @test inverse(/, intervaltype(-3.0, -1.0), -2.0)[1] ≈ intervaltype(2.0, 6.0)

    ts = IntervalVector([intervaltype(1.0, 2.0), intervaltype(2.0, 3.0)])
    @test all(inverse(/, ts, [2.0, 3.0]).intervals .≈ [intervaltype(2.0, 4.0), intervaltype(6.0, 9.0)])

    # y = 0: x/0 is undefined → no intervals produced
    @test isempty(inverse(/, intervaltype(1.0, 3.0), 0.0))
    @test isempty(inverse(/, intervaltype(-2.0, 2.0), 0.0))

    # vectorised: one zero divisor produces an empty view in that slot
    ts2 = IntervalVector([intervaltype(1.0, 2.0), intervaltype(3.0, 4.0)])
    r2 = inverse(/, ts2, [0.0, 2.0])
    @test isempty(r2[1])
    @test first(r2[2]) ≈ intervaltype(6.0, 8.0)
end

@testset "rightinverse /" begin
    @test inverse(/, 6.0, intervaltype(1.0, 3.0))[1]  ≈ intervaltype(2.0, 6.0)
    @test inverse(/, -4.0, intervaltype(1.0, 2.0))[1] ≈ intervaltype(-4.0, -2.0)

    ts = IntervalVector([intervaltype(1.0, 2.0), intervaltype(2.0, 3.0)])
    @test all(inverse(/, [4.0, 6.0], ts).intervals .≈ [intervaltype(2.0, 4.0), intervaltype(2.0, 3.0)])

    # x = 0, target contains 0 → any d works → full real line (except 0)
    @test bounds(inverse(/, 0.0, intervaltype(0.0, 2.0))[1]) == (nextfloat(-Inf), prevfloat(0.0))
    @test bounds(inverse(/, 0.0, intervaltype(0.0, 2.0))[2]) == (nextfloat(0.0), prevfloat(Inf))
    
    # x = 0, target does NOT contain 0 → impossible → sentinel
    @test isempty(inverse(/, 0.0, intervaltype(1.0, 3.0)))
    @test isempty(inverse(/, 0.0, intervaltype(-2.0, -1.0)))

    # sign-crossing target: pre-image is two disjoint rays
    rays = inverse(/, 3.0, intervaltype(-1.0, 2.0))
    @test length(rays) == 2
    # rays: (-∞, 3/(-1)] ∪ [3/2, ∞) = (-∞, -3] ∪ [1.5, ∞)
    # Unbounded ray ends are capped to ±floatmax: the forward-inner surrogate admits only
    # finite child evals (an infinite eval forwards to 0 or NaN and is handled separately).
    @test rays[1].lo ≈ -floatmax(Float64)
    @test rays[1].hi ≈ -3.0
    @test rays[2].lo ≈ 1.5
    @test rays[2].hi ≈ floatmax(Float64)

    # vectorised: mixed cases
    ts2 = IntervalVector([intervaltype(1.0, 2.0), intervaltype(2.0, 4.0)])
    r2 = inverse(/, [4.0, 0.0], ts2)
    @test r2[1][1] ≈ intervaltype(2.0, 4.0)
    @test isempty(r2[2])   # x=0, 0 ∉ [2,4]
end

@testset "evaluate sqrt" begin
    @test evaluate(sqrt, 4.0) == [2.0 for a in 0]
    @test evaluate(sqrt, [4.0, 9.0, 16.0]) == [2.0, 3.0, 4.0]
end

@testset "inverse sqrt" begin
    @test inverse(sqrt, intervaltype(1.0, 3.0))[1]  ≈ intervaltype(1.0, 9.0)
    @test inverse(sqrt, intervaltype(0.0, 2.0))[1]  ≈ intervaltype(0.0, 4.0)  
    @test inverse(sqrt, intervaltype(-1.0, 3.0))[1] ≈ intervaltype(0.0, 9.0) 

    # u < 0: unreachable target → invalid sentinel
    @test isempty(inverse(sqrt, intervaltype(-3.0, -1.0)))

    # Point interval intervaltype(2,2) → narrow(4,4): no float strictly inside → invalid_interval
    @test isempty(inverse(sqrt, intervaltype(2.0, 2.0)))

    ts = IntervalVector([intervaltype(0.0, 2.0), intervaltype(1.0, 3.0)])
    @test all(inverse(sqrt, ts).intervals .≈ [intervaltype(0.0, 4.0), intervaltype(1.0, 9.0)])
end

@testset "evaluate exp" begin
    @test evaluate(exp, 0.0) ≈ 1.0
    @test evaluate(exp, [0.0, 1.0]) ≈ [1.0, exp(1.0)]
end

@testset "inverse exp" begin
    @test inverse(exp, intervaltype(1.0, exp(3.0)))[1] ≈ intervaltype(nextfloat(0.0), 3.0)

    # Unbounded lower preimage is capped to -floatmax (forward-inner surrogate: exp(-floatmax)
    # underflows to 0, still inside the target; an infinite child eval is excluded).
    r = inverse(exp, intervaltype(-1.0, exp(2.0)))[1]
    @test r.lo ≈ nextfloat(-Inf)
    @test r.hi ≈ 2.0

    r2 = inverse(exp, intervaltype(0.0, exp(1.0)))[1]
    @test inf(r2) ≈ nextfloat(-Inf)
    @test sup(r2) ≈ 1.0

    @test isempty(inverse(exp, intervaltype(-2.0, -1.0)))
    @test isempty(inverse(exp, intervaltype(-1.0,  0.0)))

    ts = IntervalVector([intervaltype(1.0, exp(1.0)), intervaltype(exp(1.0), exp(2.0))])
    rs = inverse(exp, ts)
    @test inf(rs[1][1]) ≈ 0.0 atol=1e-10
    @test sup(rs[1][1])  ≈ 1.0
    @test inf(rs[2][1]) ≈ 1.0
    @test sup(rs[2][1])  ≈ 2.0
end

@testset "evaluate log" begin
    @test evaluate(log, 1.0) ≈ [log(v) for v in 1.0]
    @test evaluate(log, [1.0, exp(1.0)]) ≈ [0.0, 1.0]
end

@testset "inverse log" begin
    @test inverse(log, intervaltype(0.0, 1.0))[1] ≈ intervaltype(1.0, exp(1.0))
    @test inverse(log, intervaltype(-1.0, 1.0))[1] ≈ intervaltype(exp(-1.0), exp(1.0))

    r = inverse(log, intervaltype(-Inf, 1.0))[1]
    @test inf(r) ≈ 0.0 atol=1e-10
    @test sup(r) ≈ exp(1.0)

    @test inverse(log, intervaltype(-100.0, 100.0))[1] ≈ intervaltype(exp(-100.0), exp(100.0))

    ts = IntervalVector([intervaltype(0.0, 1.0), intervaltype(1.0, 2.0)])
    rs = inverse(log, ts)
    @test rs[1][1] ≈ intervaltype(1.0, exp(1.0))
    @test rs[2][1] ≈ intervaltype(exp(1.0), exp(2.0))
end

@testset "evaluate sin" begin
    @test evaluate(sin, [0.0, π/2, π]) ≈ [0.0, 1.0, 0.0] atol=1e-10
    @test isnan(evaluate(sin, [Inf])[1])
    @test isnan(evaluate(sin, [NaN])[1])
end

@testset "inverse sin" begin
    # Standard target: two branches
    r = inverse(sin, intervaltype(0.0, 0.5))
    @test r[1] isa IntervalType
    @test length(r) == 2
    # Both branches contain valid values
    for ci in r
        x = (ci.lo + ci.hi) / 2
        @test in_interval(sin(x), intervaltype(0.0, 0.5))
    end

    # Fully out of range → empty
    @test isempty(inverse(sin, intervaltype(1.5, 2.0)))
    @test isempty(inverse(sin, intervaltype(-3.0, -1.5)))

    # Upper-clamped target: peak included → single arc [asin(lo), π − asin(lo)]
    r_upper = inverse(sin, intervaltype(-0.5, 2.0))
    @test length(r_upper) == 1
    @test inf(r_upper[1]) ≈ asin(-0.5) atol=1e-10
    @test sup(r_upper[1]) ≈ π - asin(-0.5) atol=1e-10

    # Lower-clamped target: trough included → single wrapped arc [π − asin(hi), asin(hi) + 2π]
    
    target = intervaltype(-2.0, 0.5)
    r_lower = inverse(sin, target)
    for r in r_lower
        @test issubset_interval(sin(r), target)
    end
    
    # Vectorised
    ts = IntervalVector([intervaltype(0.0, 0.5), intervaltype(-0.5, 0.0)])
    rs = inverse(sin, ts)
    @test length(rs.intervals) == 4
end

