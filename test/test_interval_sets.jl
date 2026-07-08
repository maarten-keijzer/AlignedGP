using Test
using AlignedGP

@testset "CInterval construction" begin
    ci = CInterval(1.0, 3.0)
    @test ci.lo == 1.0
    @test ci.hi == 3.0

    # Constructor sorts: hi < lo swaps
    ci2 = CInterval(3.0, 1.0)
    @test ci2.lo == 1.0
    @test ci2.hi == 3.0

    @test CInterval(Inf, Inf) == CInterval(Inf, Inf)
    @test CInterval(-Inf, -Inf) == CInterval(-Inf, -Inf)
end

@testset "first / last" begin
    ci = CInterval(2.0, 5.0)
    @test first(ci) == 2.0
    @test last(ci)  == 5.0
end

@testset "equality and approx" begin
    @test CInterval(1.0, 2.0) == CInterval(1.0, 2.0)
    @test CInterval(1.0, 2.0) != CInterval(1.0, 3.0)
    @test CInterval(1.0, 2.0) ≈ CInterval(1.0, 2.0)
    @test CInterval(nextfloat(1.0), prevfloat(2.0)) ≈ CInterval(1.0, 2.0)
end

@testset "containment" begin
    ci = CInterval(1.0, 3.0)
    @test 1.0 in ci
    @test 2.0 in ci
    @test 3.0 in ci
    @test !(0.9 in ci)
    @test !(3.1 in ci)
end

@testset "invalid sentinel" begin
    @test _is_invalid(invalid_interval)
    @test _is_invalid(CInterval(1.0, NaN))     # any NaN bound is invalid
    @test !_is_invalid(CInterval(Inf, Inf))    # degenerate but valid (was old sentinel)
    @test !_is_invalid(CInterval(-Inf, -Inf))  # degenerate but valid
    @test !_is_invalid(CInterval(-Inf, Inf))   # whole real line is valid
    @test !_is_invalid(CInterval(1.0, 3.0))
end

@testset "narrow" begin
    ci = narrow(1.0, 3.0)
    @test ci.lo > 1.0
    @test ci.hi < 3.0
    @test ci ≈ CInterval(1.0, 3.0)

    # narrow handles hi < lo by swapping
    ci2 = narrow(3.0, 1.0)
    @test ci2.lo > 1.0
    @test ci2.hi < 3.0

    # infinite bounds are preserved (they denote truly unbounded)
    @test narrow(-Inf, 1.0).lo == -Inf
    @test narrow(0.0, Inf).hi  == Inf
    @test narrow(-Inf, Inf)    == CInterval(-Inf, Inf)
end

@testset "scalar arithmetic (narrowing)" begin
    ci = CInterval(2.0, 4.0)

    @test (ci + 1.0) ≈ CInterval(3.0, 5.0)
    @test (1.0 + ci) ≈ CInterval(3.0, 5.0)
    @test (ci - 1.0) ≈ CInterval(1.0, 3.0)
    @test (5.0 - ci) ≈ CInterval(1.0, 3.0)   # bounds flip: 5-4=1, 5-2=3

    # Results are strictly inside (narrowed)
    @test (ci + 1.0).lo > 3.0
    @test (ci + 1.0).hi < 5.0
end

@testset "scalar multiplication (narrowing)" begin
    ci = CInterval(2.0, 4.0)

    # positive scalar: bounds scale uniformly
    @test (ci * 2.0) ≈ CInterval(4.0, 8.0)
    @test (2.0 * ci) ≈ CInterval(4.0, 8.0)

    # negative scalar: bounds flip
    @test (ci * -1.0) ≈ CInterval(-4.0, -2.0)
    @test (-1.0 * ci) ≈ CInterval(-4.0, -2.0)

    # zero: collapses to a zero-width interval; narrow returns invalid_interval (sub-ULP)
    @test _is_invalid(ci * 0.0)

    # Results are strictly inside (narrowed) for non-degenerate cases
    @test (ci * 3.0).lo > 6.0
    @test (ci * 3.0).hi < 12.0

    # interval spanning zero with positive scalar
    ci2 = CInterval(-2.0, 3.0)
    @test (ci2 * 2.0) ≈ CInterval(-4.0, 6.0)

    # interval spanning zero with negative scalar (bounds flip)
    @test (ci2 * -2.0) ≈ CInterval(-6.0, 4.0)
end

@testset "scalar multiplication preserves invalid sentinel" begin
    s = invalid_interval
    @test _is_invalid(s * 2.0)
    @test _is_invalid(2.0 * s)
    @test _is_invalid(s * -1.0)
end

@testset "scalar arithmetic preserves invalid sentinel" begin
    s = invalid_interval
    @test _is_invalid(s + 1.0)
    @test _is_invalid(s - 1.0)
    @test _is_invalid(1.0 + s)
    @test _is_invalid(1.0 - s)
end

@testset "interval-interval arithmetic" begin
    a = CInterval(1.0, 3.0)
    b = CInterval(2.0, 4.0)
    @test (a + b) ≈ CInterval(3.0, 7.0)
    @test (a - b) ≈ CInterval(-3.0, 1.0)   # [1,3] - [2,4] = [1-4, 3-2] = [-3, 1]
end

@testset "CIntervals construction" begin
    empty = CIntervals()
    @test isempty(empty)
    @test length(empty) == 0

    single = CIntervals(CInterval(1.0, 2.0))
    @test length(single) == 1

    from_vec = CIntervals([CInterval(1.0, 2.0), CInterval(3.0, 4.0)])
    @test length(from_vec) == 2

    from_tuple = CIntervals((CInterval(1.0, 2.0), CInterval(3.0, 4.0)))
    @test from_vec == from_tuple
end

@testset "CIntervals.items returns Vector" begin
    cis = CIntervals([CInterval(1.0, 2.0), CInterval(3.0, 4.0)])
    items = cis.items
    @test items isa Vector{CInterval}
    @test items == [CInterval(1.0, 2.0), CInterval(3.0, 4.0)]
end

@testset "CIntervals containment" begin
    cis = CIntervals([CInterval(1.0, 2.0), CInterval(4.0, 5.0)])
    @test 1.5 in cis
    @test 4.5 in cis
    @test !(3.0 in cis)
end

@testset "CIntervals iteration" begin
    cis = CIntervals([CInterval(1.0, 2.0), CInterval(3.0, 4.0)])
    collected = collect(cis)
    @test collected == [CInterval(1.0, 2.0), CInterval(3.0, 4.0)]
end

@testset "flatten" begin
    v = [CIntervals((CInterval(0.1, 1.0), CInterval(2.0, 3.0))), CIntervals(CInterval(4.0, 5.0))]
    @test flatten(v) == [CInterval(0.1, 1.0), CInterval(2.0, 3.0), CInterval(4.0, 5.0)]
end

@testset "convert CInterval → CIntervals" begin
    ci = CInterval(1.0, 2.0)
    cis::CIntervals = ci
    @test length(cis) == 1
    @test first(cis) == ci
end

@testset "CIntervals scalar arithmetic" begin
    cis = CIntervals([CInterval(1.0, 2.0), CInterval(4.0, 5.0)])

    shifted = cis + 3.0
    @test length(shifted) == 2
    @test shifted ≈ CIntervals([CInterval(4.0, 5.0), CInterval(7.0, 8.0)])

    shifted2 = cis - 1.0
    @test shifted2 ≈ CIntervals([CInterval(0.0, 1.0), CInterval(3.0, 4.0)]) atol=1e-10

    flipped = 10.0 - cis   # s - [lo,hi] flips bounds: [10-hi, 10-lo]
    @test flipped ≈ CIntervals([CInterval(8.0, 9.0), CInterval(5.0, 6.0)])

    flipped2 = 3.0 + cis
    @test flipped2 ≈ CIntervals([CInterval(4.0, 5.0), CInterval(7.0, 8.0)])

    # Empty stays empty
    empty = CIntervals()
    @test isempty(empty + 5.0)
    @test isempty(empty - 5.0)
    @test isempty(5.0 - empty)
    @test isempty(5.0 + empty)
end

@testset "CIntervals scalar multiplication" begin
    cis = CIntervals([CInterval(1.0, 2.0), CInterval(4.0, 5.0)])

    scaled = cis * 3.0
    @test length(scaled) == 2
    @test scaled ≈ CIntervals([CInterval(3.0, 6.0), CInterval(12.0, 15.0)])

    # negative scalar flips bounds in each sub-interval
    flipped = cis * -1.0
    @test flipped ≈ CIntervals([CInterval(-2.0, -1.0), CInterval(-5.0, -4.0)])

    # Empty stays empty
    empty = CIntervals()
    @test isempty(empty * 2.0)
end

@testset "CIntervals broadcastable (scalar in broadcast)" begin
    cis = CIntervals([CInterval(2.0, 4.0), CInterval(6.0, 8.0)])
    # CIntervals in a vector-scalar broadcast
    v = [cis, cis]
    result = v .- [1.0, 2.0]
    @test result[1] ≈ CIntervals([CInterval(1.0, 3.0), CInterval(5.0, 7.0)])
    @test result[2] ≈ CIntervals([CInterval(0.0, 2.0), CInterval(4.0, 6.0)]) atol=1e-10
end
