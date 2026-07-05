using Test
using AlignedGP
using Random

const CI = CInterval
cis(lo, hi) = CIntervals(CI(lo, hi))   # convenience: single-interval CIntervals

@testset "max_overlap_region" begin

    @testset "empty input" begin
        r = max_overlap_region(CIntervals[])
        @test r.depth == 0
        @test isempty(r.region)
    end

    @testset "all invalid sentinels" begin
        r = max_overlap_region([cis(Inf, Inf), cis(Inf, Inf)])
        @test r.depth == 0
        @test isempty(r.region)
    end

    @testset "all empty CIntervals (error state)" begin
        r = max_overlap_region([CIntervals(), CIntervals()])
        @test r.depth == 0
        @test isempty(r.region)
    end

    @testset "mixed valid and invalid sentinels" begin
        r = max_overlap_region([cis(1.0, 3.0), cis(Inf, Inf)])
        @test r.depth == 1
        @test r.region.items == [CI(1.0, 3.0)]
    end

    @testset "-Inf sentinel filtered" begin
        r = max_overlap_region([cis(-Inf, -Inf), cis(1.0, 3.0)])
        @test r.depth == 1
        @test r.region.items == [CI(1.0, 3.0)]
    end

    @testset "single interval" begin
        r = max_overlap_region([cis(1.0, 3.0)])
        @test r.depth == 1
        @test r.region.items == [CI(1.0, 3.0)]
    end

    @testset "two disjoint intervals" begin
        r = max_overlap_region([cis(1.0, 2.0), cis(3.0, 4.0)])
        @test r.depth == 1
        @test length(r.region.items) == 2
    end

    @testset "two fully overlapping intervals" begin
        r = max_overlap_region([cis(1.0, 5.0), cis(2.0, 4.0)])
        @test r.depth == 2
        @test r.region.items == [CI(2.0, 4.0)]
    end

    @testset "touching boundaries count at shared point" begin
        r = max_overlap_region([cis(1.0, 3.0), cis(3.0, 5.0)])
        @test r.depth == 2
        @test r.region.items == [CI(3.0, 3.0)]
    end

    @testset "degenerate (point) interval participates" begin
        r = max_overlap_region([cis(1.0, 5.0), cis(3.0, 3.0)])
        @test r.depth == 2
        @test r.region.items == [CI(3.0, 3.0)]
    end

    @testset "three staggered intervals — single max region" begin
        r = max_overlap_region([cis(1.0, 5.0), cis(2.0, 6.0), cis(3.0, 7.0)])
        @test r.depth == 3
        @test r.region.items == [CI(3.0, 5.0)]
    end

    @testset "two disjoint max-depth components" begin
        r = max_overlap_region([cis(1.0, 3.0), cis(1.0, 3.0), cis(5.0, 7.0), cis(5.0, 7.0)])
        @test r.depth == 2
        @test length(r.region.items) == 2
        @test CI(1.0, 3.0) in r.region.items
        @test CI(5.0, 7.0) in r.region.items
    end

    @testset "all intervals identical" begin
        r = max_overlap_region([cis(2.0, 4.0), cis(2.0, 4.0), cis(2.0, 4.0)])
        @test r.depth == 3
        @test r.region.items == [CI(2.0, 4.0)]
    end

    @testset "CIntervals with multiple sub-intervals" begin
        # Each data point provides two disjoint target ranges
        multi = CIntervals([CI(1.0, 3.0), CI(5.0, 7.0)])
        r = max_overlap_region([multi, cis(2.0, 4.0)])
        @test r.depth == 2
        @test r.region.items == [CI(2.0, 3.0)]
    end

end

@testset "select_constant" begin

    @testset "empty region returns zero" begin
        r = max_overlap_region(CIntervals[])
        @test select_constant(r.region) == 0.0
    end

    @testset "zero in region → returns 0" begin
        r = max_overlap_region([cis(-1.0, 2.0), cis(-2.0, 1.0)])
        @test select_constant(r.region) == 0.0
    end

    @testset "integers in region → result is an integer" begin
        r = max_overlap_region([cis(1.0, 5.0), cis(3.0, 7.0)])
        rng = Random.MersenneTwister(42)
        for _ in 1:20
            c = select_constant(r.region, rng)
            @test c in [3.0, 4.0, 5.0]
        end
    end

    @testset "two components with integers → sample comes from one of the integer ranges" begin
        r = max_overlap_region([cis(1.0, 3.0), cis(1.0, 3.0), cis(5.0, 9.0), cis(5.0, 9.0)])
        rng = Random.MersenneTwister(42)
        results = [select_constant(r.region, rng) for _ in 1:4000]
        @test all(c -> (1.0 <= c <= 3.0 || 5.0 <= c <= 9.0) && c == floor(c), results)
        hits_first = count(c -> 1.0 <= c <= 3.0, results)
        @test 0.40 < hits_first / 4000 < 0.60
    end

    @testset "no integers in region → continuous sample within region" begin
        r = max_overlap_region([cis(1.2, 1.8), cis(1.2, 1.8)])
        rng = Random.MersenneTwister(42)
        for _ in 1:20
            c = select_constant(r.region, rng)
            @test 1.2 <= c <= 1.8
        end
    end

    @testset "wider component sampled more often (continuous path)" begin
        region = CIntervals([CI(0.1, 0.2), CI(0.4, 0.9)])
        rng = Random.MersenneTwister(42)
        hits_wide = count(1:10000) do _
            c = select_constant(region, rng)
            0.4 <= c <= 0.9
        end
        @test 0.78 < hits_wide / 10000 < 0.88
    end

    @testset "zero at boundary counts as in region" begin
        r = max_overlap_region([cis(-1.0, 0.0), cis(0.0, 1.0)])
        @test select_constant(r.region) == 0.0
    end

    # Half-infinite intervals: select_constant must return a value IN the region.
    # Before the fix, total=Inf caused a silent fallback to 0.0 regardless of region.

    @testset "half-infinite region CI(lo, Inf) — result must be in region" begin
        region = CIntervals(CI(2.0, Inf))  # does not contain 0
        rng = Random.MersenneTwister(42)
        for _ in 1:20
            c = select_constant(region, rng)
            @test c in region      # was: c == 0.0 which is NOT in CI(2, Inf)
        end
    end

    @testset "half-infinite region CI(-Inf, hi) — result must be in region" begin
        region = CIntervals(CI(-Inf, -3.0))  # does not contain 0
        rng = Random.MersenneTwister(42)
        for _ in 1:20
            c = select_constant(region, rng)
            @test c in region      # was: c == 0.0 which is NOT in CI(-Inf, -3)
        end
    end

    @testset "half-infinite region CI(lo, Inf) prefers integer at ceil(lo)" begin
        region = CIntervals(CI(2.3, Inf))
        c = select_constant(region)
        @test c == 3.0   # ceil(2.3) = 3, smallest integer in CI(2.3, Inf)
    end

    @testset "half-infinite region CI(-Inf, hi) prefers integer at floor(hi)" begin
        region = CIntervals(CI(-Inf, -1.7))
        c = select_constant(region)
        @test c == -2.0  # floor(-1.7) = -2, largest integer in CI(-Inf, -1.7)
    end

    @testset "two-ray region from division inverse — result in one of the rays" begin
        # _div_into produces (-∞, a] ∪ [b, +∞) when target straddles zero
        region = CIntervals([CI(-Inf, -2.0), CI(3.0, Inf)])
        rng = Random.MersenneTwister(42)
        for _ in 1:50
            c = select_constant(region, rng)
            @test c in region    # was: c == 0.0, outside both rays
        end
    end

end
