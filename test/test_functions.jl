using Test
using AlignedGP

const CI = CInterval

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
    @test leftinverse(+, CI(3.0, 5.0), 2.0)  ≈ CI(1.0, 3.0)
    @test leftinverse(+, CI(3.0, 5.0), -1.0) ≈ CI(4.0, 6.0)
    # Point interval CI(0,0) - 0.0 = narrow(0,0): no float strictly inside → invalid_interval
    @test _is_invalid(leftinverse(+, CI(0.0, 0.0), 0.0))

    ts = [CI(3.0, 5.0), CI(4.0, 6.0)]
    @test leftinverse(+, ts, [1.0, 2.0]) ≈ [CI(2.0, 4.0), CI(2.0, 4.0)]
    @test leftinverse(+, ts, 1.0)        ≈ [CI(2.0, 4.0), CI(3.0, 5.0)]
end

@testset "rightinverse +" begin
    @test rightinverse(+, CI(3.0, 5.0), 1.0)  ≈ CI(2.0, 4.0)
    @test rightinverse(+, CI(3.0, 5.0), -1.0) ≈ CI(4.0, 6.0)

    ts = [CI(3.0, 5.0), CI(4.0, 6.0)]
    @test rightinverse(+, ts, [1.0, 2.0]) ≈ [CI(2.0, 4.0), CI(2.0, 4.0)]
end

@testset "leftinverse -" begin
    @test leftinverse(-, CI(3.0, 5.0), 2.0)  ≈ CI(5.0, 7.0)
    @test leftinverse(-, CI(3.0, 5.0), -2.0) ≈ CI(1.0, 3.0)

    ts = [CI(1.0, 3.0), CI(2.0, 4.0)]
    @test leftinverse(-, ts, [1.0, 2.0]) ≈ [CI(2.0, 4.0), CI(4.0, 6.0)]
end

@testset "rightinverse -" begin
    @test rightinverse(-, CI(3.0, 5.0), 7.0) ≈ CI(2.0, 4.0)
    @test rightinverse(-, CI(1.0, 4.0), 0.0) ≈ CI(-4.0, -1.0)

    ts = [CI(1.0, 3.0), CI(2.0, 4.0)]
    @test rightinverse(-, ts, [5.0, 6.0]) ≈ [CI(2.0, 4.0), CI(2.0, 4.0)]
end

@testset "leftinverse *" begin
    @test leftinverse(*, CI(2.0, 6.0), 2.0)   ≈ CI(1.0, 3.0)
    @test leftinverse(*, CI(2.0, 6.0), -2.0)  ≈ CI(-3.0, -1.0)
    @test leftinverse(*, CI(-6.0, -2.0), -2.0) ≈ CI(1.0, 3.0)

    ts = [CI(2.0, 4.0), CI(6.0, 9.0)]
    @test leftinverse(*, ts, [2.0, 3.0]) ≈ [CI(1.0, 2.0), CI(2.0, 3.0)]
end

@testset "rightinverse *" begin
    @test rightinverse(*, CI(2.0, 6.0), 2.0)  ≈ CI(1.0, 3.0)
    @test rightinverse(*, CI(2.0, 6.0), -2.0) ≈ CI(-3.0, -1.0)

    ts = [CI(2.0, 4.0), CI(6.0, 9.0)]
    @test rightinverse(*, ts, [2.0, 3.0]) ≈ [CI(1.0, 2.0), CI(2.0, 3.0)]
end

@testset "leftinverse /" begin
    @test leftinverse(/, CI(1.0, 3.0), 2.0)   ≈ CI(2.0, 6.0)
    @test leftinverse(/, CI(1.0, 3.0), -2.0)  ≈ CI(-6.0, -2.0)
    @test leftinverse(/, CI(-3.0, -1.0), -2.0) ≈ CI(2.0, 6.0)

    ts = [CI(1.0, 2.0), CI(2.0, 3.0)]
    @test leftinverse(/, ts, [2.0, 3.0]) ≈ [CI(2.0, 4.0), CI(6.0, 9.0)]

    # y = 0: x/0 is undefined → invalid sentinel
    @test _is_invalid(leftinverse(/, CI(1.0, 3.0), 0.0))
    @test _is_invalid(leftinverse(/, CI(-2.0, 2.0), 0.0))

    # vectorised: one zero divisor produces a sentinel in that slot
    ts2 = [CI(1.0, 2.0), CI(3.0, 4.0)]
    r2 = leftinverse(/, ts2, [0.0, 2.0])
    @test _is_invalid(r2[1])
    @test r2[2] ≈ CI(6.0, 8.0)
end

@testset "rightinverse /" begin
    @test rightinverse(/, CI(1.0, 3.0), 6.0)  ≈ CI(2.0, 6.0)
    @test rightinverse(/, CI(1.0, 2.0), -4.0) ≈ CI(-4.0, -2.0)

    ts = [CI(1.0, 2.0), CI(2.0, 3.0)]
    @test rightinverse(/, ts, [4.0, 6.0]) ≈ [CI(2.0, 4.0), CI(2.0, 3.0)]

    # x = 0, target contains 0 → any d works → full real line
    @test rightinverse(/, CI(0.0, 2.0), 0.0) == CI(-Inf, Inf)
    @test rightinverse(/, CI(-1.0, 1.0), 0.0) == CI(-Inf, Inf)

    # x = 0, target does NOT contain 0 → impossible → sentinel
    @test _is_invalid(rightinverse(/, CI(1.0, 3.0), 0.0))
    @test _is_invalid(rightinverse(/, CI(-2.0, -1.0), 0.0))

    # sign-crossing target: pre-image is two disjoint rays returned as CIntervals
    r_sc = rightinverse(/, CIntervals(CI(-1.0, 2.0)), 3.0)
    @test r_sc isa CIntervals
    @test length(r_sc) == 2
    # rays: (-∞, 3/(-1)] ∪ [3/2, ∞) = (-∞, -3] ∪ [1.5, ∞)
    items = sort(r_sc.items, by = ci -> ci.lo)
    # Unbounded ray ends are capped to ±floatmax: the forward-inner surrogate admits only
    # finite child evals (an infinite eval forwards to 0 or NaN and is handled separately).
    @test items[1].lo ≈ -floatmax(Float64)
    @test items[1].hi ≈ -3.0
    @test items[2].lo ≈ 1.5
    @test items[2].hi ≈ floatmax(Float64)

    # sign-crossing with negative x
    r_sc2 = rightinverse(/, CIntervals(CI(-1.0, 2.0)), -3.0)
    @test r_sc2 isa CIntervals && length(r_sc2) == 2

    # vectorised: mixed cases
    ts2 = [CI(1.0, 2.0), CI(2.0, 4.0)]
    r2 = rightinverse(/, ts2, [4.0, 0.0])
    @test r2[1] ≈ CI(2.0, 4.0)
    @test _is_invalid(r2[2])   # x=0, 0 ∉ [2,4]
end

@testset "zero divisor / unconstrained denominator" begin
    # leftinverse: y = 0 means x/0, never finite → sentinel for any target
    @test _is_invalid(leftinverse(/, CI(1.0, 3.0), 0.0))
    @test _is_invalid(leftinverse(/, CI(-1.0, 1.0), 0.0))

    # rightinverse: x = 0 → result is always 0, so constrained iff target contains 0
    @test rightinverse(/, CI(-1.0, 1.0), 0.0) == CI(-Inf, Inf)  # 0 ∈ target → unconstrained
    @test _is_invalid(rightinverse(/, CI(1.0, 3.0), 0.0))        # 0 ∉ target → impossible

    # rightinverse: single-interval sign-crossing target → two rays (not a sentinel)
    r = rightinverse(/, CIntervals(CI(-1.0, 1.0)), 5.0)
    @test r isa CIntervals && length(r) == 2
    @test !_is_invalid(r._solo)   # result is valid (two rays)
end

@testset "evaluate identity" begin
    @test evaluate(identity, 3.0) == 3.0
    @test evaluate(identity, [1.0, 2.0, 3.0]) == [1.0, 2.0, 3.0]
end

@testset "inverse identity" begin
    @test inverse(identity, CI(1.0, 3.0)) == CI(1.0, 3.0)
    @test inverse(identity, CI(-Inf, Inf)) == CI(-Inf, Inf)
    @test isnan(inverse(identity, invalid_interval).lo) # invalid sentinel passes through unchanged
    @test isnan(inverse(identity, invalid_interval).hi) # invalid sentinel passes through unchanged

    ts = [CI(1.0, 2.0), CI(3.0, 4.0)]
    @test inverse(identity, ts) == ts
end

@testset "evaluate sqrt" begin
    @test evaluate(sqrt, 4.0) == [2.0 for a in 0]
    @test evaluate(sqrt, [4.0, 9.0, 16.0]) == [2.0, 3.0, 4.0]
end

@testset "inverse sqrt" begin
    @test inverse(sqrt, CI(1.0, 3.0))  ≈ CI(1.0, 9.0)
    @test inverse(sqrt, CI(0.0, 2.0))  ≈ CI(0.0, 4.0)  atol=1e-10
    @test inverse(sqrt, CI(-1.0, 3.0)) ≈ CI(0.0, 9.0) atol=1e-10

    # u < 0: unreachable target → invalid sentinel
    @test _is_invalid(inverse(sqrt, CI(-3.0, -1.0)))

    # Point interval CI(2,2) → narrow(4,4): no float strictly inside → invalid_interval
    @test _is_invalid(inverse(sqrt, CI(2.0, 2.0)))

    ts = [CI(0.0, 2.0), CI(1.0, 3.0)]
    @test inverse(sqrt, ts) ≈ [CI(0.0, 4.0), CI(1.0, 9.0)] atol=1e-10
end

@testset "evaluate exp" begin
    @test evaluate(exp, 0.0) ≈ 1.0
    @test evaluate(exp, [0.0, 1.0]) ≈ [1.0, exp(1.0)]
end

@testset "inverse exp" begin
    @test inverse(exp, CI(1.0, exp(3.0))) ≈ CI(0.0, 3.0) atol=1e-10

    # Unbounded lower preimage is capped to -floatmax (forward-inner surrogate: exp(-floatmax)
    # underflows to 0, still inside the target; an infinite child eval is excluded).
    r = inverse(exp, CI(-1.0, exp(2.0)))
    @test first(r) ≈ -floatmax(Float64)
    @test last(r) ≈ 2.0

    r2 = inverse(exp, CI(0.0, exp(1.0)))
    @test first(r2) ≈ -floatmax(Float64)
    @test last(r2) ≈ 1.0

    @test _is_invalid(inverse(exp, CI(-2.0, -1.0)))
    @test _is_invalid(inverse(exp, CI(-1.0,  0.0)))

    ts = [CI(1.0, exp(1.0)), CI(exp(1.0), exp(2.0))]
    rs = inverse(exp, ts)
    @test first(rs[1]) ≈ 0.0 atol=1e-10
    @test last(rs[1])  ≈ 1.0
    @test first(rs[2]) ≈ 1.0
    @test last(rs[2])  ≈ 2.0
end

@testset "evaluate log" begin
    @test evaluate(log, 1.0) ≈ [log(v) for v in 1.0]
    @test evaluate(log, [1.0, exp(1.0)]) ≈ [0.0, 1.0]
end

@testset "inverse log" begin
    @test inverse(log, CI(0.0, 1.0)) ≈ CI(1.0, exp(1.0))
    @test inverse(log, CI(-1.0, 1.0)) ≈ CI(exp(-1.0), exp(1.0))

    r = inverse(log, CI(-Inf, 1.0))
    @test first(r) ≈ 0.0 atol=1e-10
    @test last(r) ≈ exp(1.0)

    @test inverse(log, CI(-100.0, 100.0)) ≈ CI(exp(-100.0), exp(100.0))

    ts = [CI(0.0, 1.0), CI(1.0, 2.0)]
    rs = inverse(log, ts)
    @test rs[1] ≈ CI(1.0, exp(1.0))
    @test rs[2] ≈ CI(exp(1.0), exp(2.0))
end

@testset "invalid sentinel propagates through scalar arithmetic" begin
    sentinel = invalid_interval
    @test _is_invalid(sentinel - 5.0)
    @test _is_invalid(sentinel + 3.0)
    @test _is_invalid(2.0 - sentinel)
    @test _is_invalid(_scale(sentinel, 2.0))
    @test _is_invalid(_div_into(3.0, sentinel))
end

@testset "evaluate cos / sin" begin
    @test evaluate(cos, [0.0, π/2, π]) ≈ [1.0, 0.0, -1.0] atol=1e-10
    @test evaluate(sin, [0.0, π/2, π]) ≈ [0.0, 1.0, 0.0] atol=1e-10
    @test isnan(evaluate(cos, [Inf])[1])
    @test isnan(evaluate(sin, [Inf])[1])
    @test isnan(evaluate(cos, [-Inf])[1])
    @test isnan(evaluate(sin, [NaN])[1])
end

@testset "inverse cos" begin
    # Standard target in (−1, 1): two branches in [−π, π]
    r = inverse(cos, CI(0.0, 0.5))
    @test r isa CIntervals
    @test length(r) == 2
    items = r.items
    # Branch 1 in [0, π]: x ∈ [acos(0.5), acos(0.0)] = [π/3, π/2]
    @test any(ci -> isapprox(ci.lo, acos(0.5), atol=1e-9) || isapprox(ci.hi, acos(0.0), atol=1e-9), items)
    # Both branches contain valid values: cos of a point in each branch hits the target
    for ci in r
        x = (ci.lo + ci.hi) / 2
        @test -1 <= cos(x) <= 1
        @test cos(x) ∈ CI(0.0, 0.5)
    end

    # Fully out of range → empty
    @test isempty(inverse(cos, CI(1.5, 2.0)))
    @test isempty(inverse(cos, CI(-3.0, -1.5)))

    # Partial clamp: [0.5, 2.0] → effective [0.5, 1.0]
    r_clamp = inverse(cos, CI(0.5, 2.0))
    r_exact = inverse(cos, CI(0.5, 1.0))
    @test length(r_clamp) == length(r_exact)
    @test r_clamp ≈ r_exact

    # Full range [−1, 1]: branches cover [0,π] and [−π,0], each touching 0
    r_full = inverse(cos, CI(-1.0, 1.0))
    @test length(r_full) == 2
    for ci in r_full
        @test isapprox(ci.lo, 0.0, atol=1e-10) || isapprox(ci.hi, 0.0, atol=1e-10)
    end

    # Vectorised call returns Vector{CIntervals}
    ts = [CI(0.0, 0.5), CI(-1.0, 0.0)]
    rs = inverse(cos, ts)
    @test rs isa Vector
    @test all(r -> r isa CIntervals && length(r) == 2, rs)
end

@testset "inverse sin" begin
    # Standard target: two branches
    r = inverse(sin, CI(0.0, 0.5))
    @test r isa CIntervals
    @test length(r) == 2
    # Both branches contain valid values
    for ci in r
        x = (ci.lo + ci.hi) / 2
        @test sin(x) ∈ CI(0.0, 0.5)
    end

    # Fully out of range → empty
    @test isempty(inverse(sin, CI(1.5, 2.0)))
    @test isempty(inverse(sin, CI(-3.0, -1.5)))

    # Upper-clamped target: peak included → single arc [asin(lo), π − asin(lo)]
    r_upper = inverse(sin, CI(-0.5, 2.0))
    @test length(r_upper) == 1
    @test r_upper._solo.lo ≈ asin(-0.5) atol=1e-10
    @test r_upper._solo.hi ≈ π - asin(-0.5) atol=1e-10

    # Lower-clamped target: trough included → single wrapped arc [π − asin(hi), asin(hi) + 2π]
    r_lower = inverse(sin, CI(-2.0, 0.5))
    @test length(r_lower) == 1
    @test r_lower._solo.lo ≈ π - asin(0.5) atol=1e-10
    @test r_lower._solo.hi ≈ asin(0.5) + 2π atol=1e-10

    # Vectorised
    ts = [CI(0.0, 0.5), CI(-0.5, 0.0)]
    rs = inverse(sin, ts)
    @test all(r -> r isa CIntervals && length(r) == 2, rs)
end

@testset "inverse cos/sin CIntervals lift" begin
    # CIntervals input with one sub-interval → 2 output branches
    cis = CIntervals(CI(0.0, 0.5))
    r = inverse(cos, cis)
    @test r isa CIntervals
    @test length(r) == 2

    # Empty input → empty output
    @test isempty(inverse(cos, CIntervals()))
    @test isempty(inverse(sin, CIntervals()))

    # Out-of-range sub-interval is dropped
    cis2 = CIntervals([CI(0.0, 0.5), CI(2.0, 3.0)])
    r2 = inverse(cos, cis2)
    @test length(r2) == 2   # only the valid sub-interval contributes branches
end

@testset "preimage_circular sin" begin
    # Whole circle
    r = AlignedGP._preimage_circular(sin, CI(-1.0, 1.0))
    @test length(r) == 1
    @test r[1][1] ≈ 0.0 atol=1e-10
    @test r[1][2] ≈ 2π  atol=1e-10

    # Out of range → empty
    @test isempty(AlignedGP._preimage_circular(sin, CI(1.5, 2.0)))
    @test isempty(AlignedGP._preimage_circular(sin, CI(-3.0, -1.5)))

    # Clamped target (includes peak): two arcs merge into one
    r_peak = AlignedGP._preimage_circular(sin, CI(0.5, 2.0))
    @test length(r_peak) == 1
    x = (r_peak[1][1] + r_peak[1][2]) / 2
    @test sin(x) >= 0.5 - 1e-10

    # Clamped target (includes trough): arcs connect through 3π/2
    r_trough = AlignedGP._preimage_circular(sin, CI(-2.0, -0.5))
    @test length(r_trough) == 1
    x = (r_trough[1][1] + r_trough[1][2]) / 2
    @test sin(x) <= -0.5 + 1e-10

    # Fully inside [-1,1]: rising branch wraps across 0/2π → 3 arcs on [0,2π).
    # Arc layout: [0,π/6] (rising near 0), [5π/6,7π/6] (falling), [11π/6,2π] (rising near 2π).
    r2 = AlignedGP._preimage_circular(sin, CI(-0.5, 0.5))
    @test length(r2) == 3
    for (a, b) in r2
        x = (a + b) / 2
        @test -0.5 <= sin(x) <= 0.5 + 1e-10
    end

    # All arcs are within [0, 2π)
    for target in [CI(-0.5, 0.5), CI(0.0, 1.0), CI(-1.0, 0.0)]
        for (a, b) in AlignedGP._preimage_circular(sin, target)
            @test 0.0 <= a <= 2π + 1e-10
            @test 0.0 <= b <= 2π + 1e-10
            @test a <= b
        end
    end
end

@testset "preimage_circular cos" begin
    # Whole circle
    r = AlignedGP._preimage_circular(cos, CI(-1.0, 1.0))
    @test length(r) == 1
    @test r[1][1] ≈ 0.0 atol=1e-10
    @test r[1][2] ≈ 2π  atol=1e-10

    # Out of range → empty
    @test isempty(AlignedGP._preimage_circular(cos, CI(1.5, 2.0)))

    # Fully inside [-1,1]: two arcs
    r2 = AlignedGP._preimage_circular(cos, CI(-0.5, 0.5))
    @test length(r2) == 2
    for (a, b) in r2
        x = (a + b) / 2
        @test -0.5 <= cos(x) <= 0.5 + 1e-10
    end
end
