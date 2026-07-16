using Test
using AlignedGP
using AlignedGP.ReverseIntervals
using Random

using IntervalArithmetic: issubset_interval, in_interval

@testset "max_overlap_region" begin

    @testset "empty input" begin
        r = max_overlap_region(IntervalType[])
        @test r.depth == 0
        @test isempty(r.region)
    end

    @testset "single interval" begin
        r = max_overlap_region([intervaltype(1.0, 3.0)])
        @test r.depth == 1
        @test isequal_interval( first(r.region), intervaltype(1.0, 3.0))
    end

    @testset "two disjoint intervals" begin
        r = max_overlap_region([intervaltype(1.0, 2.0), intervaltype(3.0, 4.0)])
        @test r.depth == 1
        @test length(r.region) == 2
    end

    @testset "two fully overlapping intervals" begin
        r = max_overlap_region([intervaltype(1.0, 5.0), intervaltype(2.0, 4.0)])
        @test r.depth == 2
        @test length(r.region) == 1
        @test isequal_interval(first(r.region),  intervaltype(2.0, 4.0))
    end

    @testset "touching boundaries count at shared point" begin
        r = max_overlap_region([intervaltype(1.0, 3.0), intervaltype(3.0, 5.0)])
        @test r.depth == 2
        @test length(r.region) == 1
        @test isequal_interval(first(r.region), intervaltype(3.0, 3.0))
    end

    @testset "degenerate (point) interval participates" begin
        r = max_overlap_region([intervaltype(1.0, 5.0), intervaltype(3.0, 3.0)])
        @test r.depth == 2
        @test length(r.region) == 1
        @test isequal_interval(first(r.region), intervaltype(3.0, 3.0))
    end

    @testset "three staggered intervals — single max region" begin
        r = max_overlap_region([intervaltype(1.0, 5.0), intervaltype(2.0, 6.0), intervaltype(3.0, 7.0)])
        @test r.depth == 3
        @test isequal_interval(first(r.region), intervaltype(3.0, 5.0))
    end

    @testset "two disjoint max-depth components" begin
        r = max_overlap_region([intervaltype(1.0, 3.0), intervaltype(1.0, 3.0), intervaltype(5.0, 7.0), intervaltype(5.0, 7.0)])
        @test r.depth == 2
        @test length(r.region) == 2
        @test any(item -> isequal_interval(intervaltype(1.0, 3.0), item), r.region)
        @test any(item -> isequal_interval(intervaltype(5.0, 7.0), item), r.region)
    end

    @testset "all intervals identical" begin
        r = max_overlap_region([intervaltype(2.0, 4.0), intervaltype(2.0, 4.0), intervaltype(2.0, 4.0)])
        @test r.depth == 3
        @test any(item -> isequal_interval(intervaltype(2.0, 4.0), item), r.region)
    end
end;

@testset "select_constant" begin

    @testset "empty region returns zero" begin
        r = max_overlap_region(IntervalType[])
        @test select_constant(r.region) == 0.0
    end

    @testset "zero in region → returns 0" begin
        r = max_overlap_region([intervaltype(-1.0, 2.0), intervaltype(-2.0, 1.0)])
        @test select_constant(r.region) == 0.0
    end

    @testset "region without zero → continuous sample within overlap" begin
        r = max_overlap_region([intervaltype(1.0, 5.0), intervaltype(3.0, 7.0)])  # overlap [3, 5]
        rng = Random.MersenneTwister(42)
        for _ in 1:20
            c = select_constant(r.region, rng)
            @test 3.0 <= c <= 5.0
        end
    end

    @testset "two components → sample lands in a component, weighted by Cauchy prior mass" begin
        r = max_overlap_region([intervaltype(1.0, 3.0), intervaltype(1.0, 3.0), intervaltype(5.0, 9.0), intervaltype(5.0, 9.0)])
        rng = Random.MersenneTwister(42)
        results = [select_constant(r.region, rng) for _ in 1:4000]
        @test all(c -> (1.0 <= c <= 3.0 || 5.0 <= c <= 9.0), results)
        hits_first = count(c -> 1.0 <= c <= 3.0, results)
        # Selection is by Cauchy(0, γ=1) prior mass tempered by α=0.5, NOT by width.
        # w ∝ (atan(hi) - atan(lo))^0.5: [1,3] → 0.681, [5,9] → 0.295, so the
        # closer-to-zero component [1,3] is chosen ~0.70 of the time despite [5,9]
        # being twice as wide (the Cauchy tail downweights large magnitudes).
        @test 0.66 < hits_first / 4000 < 0.74
    end

    @testset "no integers in region → continuous sample within region" begin
        r = max_overlap_region([intervaltype(1.2, 1.8), intervaltype(1.2, 1.8)])
        rng = Random.MersenneTwister(42)
        for _ in 1:20
            c = select_constant(r.region, rng)
            @test 1.2 <= c <= 1.8
        end
    end

    @testset "zero at boundary counts as in region" begin
        r = max_overlap_region([intervaltype(-1.0, 0.0), intervaltype(0.0, 1.0)])
        @test select_constant(r.region) == 0.0
    end

    # Half-infinite intervals: select_constant must return a value IN the region.
    # Before the fix, total=Inf caused a silent fallback to 0.0 regardless of region.

    @testset "half-infinite region CI(lo, Inf) — result must be in region" begin
        region = [intervaltype(2.0, Inf)]  # does not contain 0
        rng = Random.MersenneTwister(42)
        for _ in 1:20
            c = select_constant(region, rng)
            @test in_interval(c, region[1])      # was: c == 0.0 which is NOT in CI(2, Inf)
        end
    end

    @testset "half-infinite region CI(-Inf, hi) — result must be in region" begin
        region = [intervaltype(-Inf, -3.0)]  # does not contain 0
        rng = Random.MersenneTwister(42)
        for _ in 1:20
            c = select_constant(region, rng)
            @test in_interval(c, region[1])      # was: c == 0.0 which is NOT in CI(-Inf, -3)
        end
    end

    @testset "half-infinite region CI(lo, Inf) → returns truncated Cauchy" begin
        region = [intervaltype(2.3, Inf)]  # width Inf → fall back to the finite bound
        c = select_constant(region)
        @test c >= 2.3
    end

    @testset "half-infinite region CI(-Inf, hi) → returns truncated Cauchy" begin
        region = [intervaltype(-Inf, -1.7)]
        c = select_constant(region)
        @test c <= -1.7
    end

    @testset "two-ray region from division inverse — result in one of the rays" begin
        # _div_into produces (-∞, a] ∪ [b, +∞) when target straddles zero
        region = [intervaltype(-Inf, -2.0), intervaltype(3.0, Inf)]
        rng = Random.MersenneTwister(42)
        for _ in 1:50
            c = select_constant(region, rng)
            @test any(item -> in_interval(c, item), region)
        end
    end

    @testset "very large finite intervals whose widths sum to Inf — no crash" begin
        # Two fully-finite intervals that don't contain 0, both with abs(bound) >= 9e18.
        # Their widths each approach prevfloat(Inf), so sum(widths) overflows to Inf.
        # Before the fix: finite_bounds stayed empty → rand(rng, []) crashed.
        big = prevfloat(Inf)
        region = [intervaltype(1.0e19, big), intervaltype(2.0e19, big)]
        rng = Random.MersenneTwister(42)
        for _ in 1:10
            c = select_constant(region, rng)
            @test any(in_interval.(Ref(c), region))
        end
    end

end
