using AlignedGP.ReverseIntervals
using Test

@testset "sqrt_rev" begin 
    intervals = [intervaltype(0,1), intervaltype(-4, -1), intervaltype(-0.5, 0.8)]
    iv = IntervalVector(intervals)

    @test eachindex(iv) == 1:3
    for i in eachindex(iv)
        @test length(iv[i]) == 1
        @test isequal_interval(first(iv[i]), intervals[i])
    end

    iv2 = invert(iv, sqrt_rev)
    @test length(iv2) == 3

    for i in eachindex(iv2)
        for intv in iv2[i]
            @test issubset_interval(sqrt(intv), intervals[i])
        end
    end
end;

@testset "inv_rev" begin 
    intervals = [intervaltype(0,1), intervaltype(-4, -1), intervaltype(-0.5, 0.8)]
    iv = IntervalVector(intervals)
    iv2 = invert(iv, inv_rev)
    @test length(iv2) == 3

    for i in eachindex(iv2)
        for intv in iv2[i]
            @test issubset_interval(inv(intv), intervals[i])
        end
    end
end;

@testset "simple addition" begin
    intervals = [intervaltype(0,1), intervaltype(-4, -1), intervaltype(-0.5, 0.8)]
    arg = [0.1, -0.3, 0.5]
    iv = IntervalVector(intervals)
    iv2 = invert(iv, arg, add_rev)

    @test length(iv2) == 3

    for i in eachindex(iv2)
        for intv in iv2[i] 
            @test issubset_interval(intv + arg[i], intervals[i])
        end
    end
end

@testset "Addition with NaN arguments" begin
    intervals = [intervaltype(0,1), intervaltype(-4, -1), intervaltype(-0.5, 0.8)]
    arg = [Inf, -0.3, NaN]
    iv = IntervalVector(intervals)
    iv2 = invert(iv, arg, add_rev)

    @test length(iv2) == 3

    for i in eachindex(iv2)
        for intv in iv2[i] 
            @test issubset_interval(intv + arg[i], intervals[i])
        end
    end
end;

@testset "simple multiplication" begin
    intervals = [intervaltype(0,1), intervaltype(-4, -1), intervaltype(-0.5, 0.8)]
    arg = [0.1, -0.3, 0.5]
    iv = IntervalVector(intervals)
    iv2 = invert(iv, arg, mul_rev)

    @test length(iv2) == 3

    for i in eachindex(iv2)
        for intv in iv2[i] 
            @test issubset_interval(intv * arg[i], intervals[i])
        end
    end
end;

@testset "Multiplication with non-finite arguments" begin
    intervals = [intervaltype(0,1), intervaltype(-4, -1), intervaltype(-0.5, 0.8)]
    arg = [Inf, -0.3, NaN]
    iv = IntervalVector(intervals)
    iv2 = invert(iv, arg, mul_rev)

    @test length(iv2) == 3

    for i in eachindex(iv2)
        for intv in iv2[i] 
            @test issubset_interval(intv * arg[i], intervals[i])
        end
    end
end;

@testset "sin_rev" begin 
    intervals = [intervaltype(0,1), intervaltype(-4, -1), intervaltype(-0.5, 0.8)]
    iv = IntervalVector(intervals)
    iv2 = invert(iv, sin_rev)
    @test length(iv2) == 3

    for i in eachindex(iv2)
        for intv in iv2[i]
            @test issubset_interval(sin(intv), intervals[i])
        end
    end
end;

@testset "exp_rev" begin 
    intervals = [intervaltype(0,1), intervaltype(-4, -1), intervaltype(-0.5, 0.8)]
    iv = IntervalVector(intervals)
    iv2 = invert(iv, exp_rev)
    @test length(iv2) == 3

    for i in eachindex(iv2)
        for intv in iv2[i]
            @test issubset_interval(exp(intv), intervals[i])
        end
    end
end;

@testset "log_rev" begin 
    intervals = [intervaltype(0,1), intervaltype(-4, -1), intervaltype(-0.5, 0.8)]
    iv = IntervalVector(intervals)
    iv2 = invert(iv, log_rev)
    @test length(iv2) == 3

    for i in eachindex(iv2)
        for intv in iv2[i]
            @test issubset_interval(log(intv), intervals[i])
        end
    end
end;

