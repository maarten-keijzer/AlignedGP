using Test
using AlignedGP
using Intervals

const CI = Interval{Float64,Closed,Closed}

@testset "max_overlap_region" begin

    @testset "empty input" begin
        r = max_overlap_region(CI[])
        @test r.depth == 0
        @test isempty(r.region.items)
    end

    @testset "all invalid sentinels" begin
        r = max_overlap_region([CI(Inf, Inf), CI(Inf, Inf)])
        @test r.depth == 0
        @test r.region.items == [CI(Inf, Inf)]
    end

    @testset "mixed valid and invalid sentinels" begin
        # Invalid sentinel should be filtered; valid interval contributes normally
        r = max_overlap_region([CI(1.0, 3.0), CI(Inf, Inf)])
        @test r.depth == 1
        @test r.region.items == [CI(1.0, 3.0)]
    end

    @testset "-Inf sentinel filtered" begin
        r = max_overlap_region([CI(-Inf, -Inf), CI(1.0, 3.0)])
        @test r.depth == 1
        @test r.region.items == [CI(1.0, 3.0)]
    end

    @testset "single interval" begin
        r = max_overlap_region([CI(1.0, 3.0)])
        @test r.depth == 1
        @test r.region.items == [CI(1.0, 3.0)]
    end

    @testset "two disjoint intervals" begin
        # max depth is 1; two separate components
        r = max_overlap_region([CI(1.0, 2.0), CI(3.0, 4.0)])
        @test r.depth == 1
        @test length(r.region.items) == 2
    end

    @testset "two fully overlapping intervals" begin
        r = max_overlap_region([CI(1.0, 5.0), CI(2.0, 4.0)])
        @test r.depth == 2
        @test r.region.items == [CI(2.0, 4.0)]
    end

    @testset "touching boundaries count at shared point" begin
        # [1,3] and [3,5] share only point 3 — depth 2 there
        r = max_overlap_region([CI(1.0, 3.0), CI(3.0, 5.0)])
        @test r.depth == 2
        @test r.region.items == [CI(3.0, 3.0)]
    end

    @testset "degenerate (point) interval participates" begin
        r = max_overlap_region([CI(1.0, 5.0), CI(3.0, 3.0)])
        @test r.depth == 2
        @test r.region.items == [CI(3.0, 3.0)]
    end

    @testset "three staggered intervals — single max region" begin
        # [1,5],[2,6],[3,7]: depth 3 only in [3,5]
        r = max_overlap_region([CI(1.0, 5.0), CI(2.0, 6.0), CI(3.0, 7.0)])
        @test r.depth == 3
        @test r.region.items == [CI(3.0, 5.0)]
    end

    @testset "two disjoint max-depth components" begin
        # [1,3],[1,3],[5,7],[5,7]: depth 2 in [1,3] and [5,7]
        r = max_overlap_region([CI(1.0, 3.0), CI(1.0, 3.0), CI(5.0, 7.0), CI(5.0, 7.0)])
        @test r.depth == 2
        @test length(r.region.items) == 2
        @test CI(1.0, 3.0) in r.region.items
        @test CI(5.0, 7.0) in r.region.items
    end

    @testset "all intervals identical" begin
        r = max_overlap_region([CI(2.0, 4.0), CI(2.0, 4.0), CI(2.0, 4.0)])
        @test r.depth == 3
        @test r.region.items == [CI(2.0, 4.0)]
    end

end

@testset "select_constant" begin

    @testset "empty region returns zero" begin
        r = max_overlap_region(CI[])
        @test select_constant(r.region) == 0.0
    end

    @testset "zero in region → returns 0" begin
        r = max_overlap_region([CI(-1.0, 2.0), CI(-2.0, 1.0)])
        @test select_constant(r.region) == 0.0
    end

    @testset "zero not in region → midpoint of widest component" begin
        # depth-2 region is [3,5], midpoint 4.0
        r = max_overlap_region([CI(1.0, 5.0), CI(3.0, 7.0)])
        @test select_constant(r.region) == 4.0
    end

    @testset "two components, widest chosen" begin
        # [1,3] width 2, [5,9] width 4 → midpoint of [5,9] = 7.0
        r = max_overlap_region([CI(1.0, 3.0), CI(1.0, 3.0), CI(5.0, 9.0), CI(5.0, 9.0)])
        @test select_constant(r.region) == 7.0
    end

    @testset "zero at boundary counts as in region" begin
        r = max_overlap_region([CI(-1.0, 0.0), CI(0.0, 1.0)])
        @test select_constant(r.region) == 0.0
    end

end
