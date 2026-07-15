# --- Tests -----------------------------------------------------------------
using AlignedGP.ReverseIntervals
using Test 
using IntervalArithmetic: issubset_interval, in_interval, sup, inf

@testset "sqrt_rev" begin
    @testset "normal range" begin
        y = intervaltype(3.00100, 101238874.11232)
        x = sqrt_rev(y)
        @test issubset_interval(sqrt(x), y)
    end

    @testset "negative lower bound" begin
        y = intervaltype(-3, 101233.033)
        x = sqrt_rev(y)
        @test issubset_interval(sqrt(x), y)
    end

    @testset "point target has no inner approx" begin
        y = intervaltype(1000.01)
        x = sqrt_rev(y)
        @test isnothing(x)
        # @test issubset_interval(sqrt(x), y)
        # @test isempty_interval(x)
    end

    @testset "entirely negative is empty" begin
        y = intervaltype(-2, -1)
        x = sqrt_rev(y)
        @test isnothing(x)
        #@test issubset_interval(sqrt(x), y)
        #@test isempty_interval(x)
    end

    @testset "unbounded upper (Inf)" begin
        y = intervaltype(3, Inf)
        x = sqrt_rev(y)
        @test issubset_interval(sqrt(x), y)
    end

    @testset "yh==0 keeps x==0 solution" begin
        for y in (intervaltype(-2.0, 0.0), intervaltype(0.0, 0.0))
            x = sqrt_rev(y)
            @test !isnothing(x)
            @test issubset_interval(sqrt(x), y)
        end
    end
end;

@testset "add_rev" begin
    @testset "positive shift" begin
        z = intervaltype(2, 5)
        x = add_rev(z, 3.0)
        @test issubset_interval(x + 3.0, z)
    end

    @testset "negative shift" begin
        z = intervaltype(-2, 12)
        x = add_rev(z, -3.0)
        @test issubset_interval(x + -3.0, z)
    end

    @testset "non-finite y is empty" begin
        z = intervaltype(-2, 12)
        x = add_rev(z, Inf)
        @test isnothing(x)
    end

    @testset "unbounded z" begin
        z = intervaltype(-Inf, 12)
        x = add_rev(z, 3.0)
        @test issubset_interval(x + 3.0, z)
    end

    # RED: subtraction by a scalar is exact when no rounding occurs, so an
    # exactly-representable preimage must survive rather than be narrowed away.
    @testset "exact shift keeps its endpoints" begin
        z = intervaltype(2.0, 5.0)          # true preimage is [-1, 2]
        x = add_rev(z, 3.0)
        @test issubset_interval(x + 3.0, z)
        @test in_interval(-1.0, x)       # -1 + 3 = 2 ∈ z
        @test in_interval(2.0, x)        #  2 + 3 = 5 ∈ z
    end

    @testset "point target has an exact preimage" begin
        z = intervaltype(5.0)               # x = 2.0 solves 2 + 3 = 5 exactly
        x = add_rev(z, 3.0)
        @test !isnothing(x)
        @test issubset_interval(x + 3.0, z)
    end
end;

@testset "mul_rev" begin
    @testset "positive multiplier" begin
        z = intervaltype(2, 5)
        x = mul_rev(z, 3.0)
        @test issubset_interval(x * 3.0, z)
    end

    @testset "negative multiplier" begin
        z = intervaltype(-2, 12)
        x = mul_rev(z, -3.0)
        @test issubset_interval(x * -3.0, z)
    end

    @testset "y==0 with 0 ∈ z gives all reals" begin
        z = intervaltype(-2, 12)
        x = mul_rev(z, 0.0)
        @test issubset_interval(x * 0.0, z)
    end

    @testset "y==0 with 0 ∉ z is empty" begin
        for z in (intervaltype(2, 5), intervaltype(-5, -2))
            x = mul_rev(z, 0.0)
            @test isnothing(x)
        end
    end

    @testset "non-finite y is empty" begin
        z = intervaltype(-2, 12)
        x = mul_rev(z, Inf)
        @test isnothing(x)
    end

    @testset "unbounded z" begin
        z = intervaltype(-Inf, 12)
        x = mul_rev(z, 3.0)
        @test issubset_interval(x * 3.0, z)
    end

    @testset "point target has no inner approx" begin
        z = intervaltype(5.0)
        x = mul_rev(z, 3.0)
        @test isnothing(x)
    end
end;

@testset "inv_rev" begin
    @testset "positive interval" begin
        y = intervaltype(1, 3)
        x = inv_rev(y)
        @test issubset_interval(inv(x), y)
    end

    @testset "negative interval" begin
        y = intervaltype(-3, -1)
        x = inv_rev(y)
        @test issubset_interval(inv(x), y)
    end

    @testset "lower bound zero" begin
        y = intervaltype(0, 3)
        x = inv_rev(y)
        @test issubset_interval(inv(x), y)
    end

    @testset "upper bound zero" begin
        y = intervaltype(-3, 0)
        x = inv_rev(y)
        @test issubset_interval(inv(x), y)
    end

    @testset "straddling zero splits into two branches" begin
        y = intervaltype(-1, 5)
        x1, x2 = inv_rev(y)
        @test issubset_interval(inv(x1), y)
        @test issubset_interval(inv(x2), y)
    end

    @testset "y==[0,0] has no finite preimage" begin
        y = intervaltype(0.0, 0.0)
        x = inv_rev(y)
        @test isnothing(x)
    end
end;

@testset "sin_rev" begin
    # Normalise the Union{Nothing,Interval,Tuple} return into a vector of arcs.
    arcs(r) = r === nothing ? IntervalType[] : (r isa Tuple ? collect(r) : [r])

    @testset "interior target: two disjoint arcs, both inner approximations" begin
        y = intervaltype(0.2, 0.8)
        x = sin_rev(y)
        @test x isa Tuple                       # rising + falling arc
        for a in arcs(x)
            @test issubset_interval(sin(a), y)  # inner: sin(arc) ⊆ y
        end
    end

    @testset "negative interior target stays inner" begin
        y = intervaltype(-0.8, -0.2)
        for a in arcs(sin_rev(y))
            @test issubset_interval(sin(a), y)
        end
    end

    @testset "peak covered (yh > 1): single merged arc through π/2" begin
        y = intervaltype(0.5, 1.5)
        x = sin_rev(y)
        @test !(x isa Tuple)                    # one connected arc, no seam at peak
        @test issubset_interval(sin(only(arcs(x))), y)
    end

    @testset "yh == 1 exactly merges at the peak (>= boundary)" begin
        y = intervaltype(0.99, 1.0)
        x = sin_rev(y)
        @test !(x isa Tuple)
        @test issubset_interval(sin(only(arcs(x))), y)
    end

    @testset "full range [-1,1] is a single covering arc" begin
        y = intervaltype(-1.0, 1.0)
        x = sin_rev(y)
        @test !(x isa Tuple)
        @test issubset_interval(sin(only(arcs(x))), y)
    end

    @testset "target entirely outside [-1,1] is empty" begin
        for y in (intervaltype(-2.0, -1.5), intervaltype(1.5, 3.0))
            @test isnothing(sin_rev(y))
        end
    end

    @testset "interior point target has no inner approximation" begin
        # Mirrors sqrt_rev/mul_rev: a non-representable point crossing yields nothing.
        y = intervaltype(0.5, 0.5)
        @test isnothing(sin_rev(y))
    end

    @testset "peak point [1,1] has no inner approximation" begin
        y = intervaltype(1.0, 1.0)
        @test isnothing(sin_rev(y))
    end

    @testset "trough covered splits across the seam (two arcs)" begin
        y = intervaltype(-1.5, -0.5)          # trough included, peak not
        x = sin_rev(y)
        @test x isa Tuple                     # two arcs at the window edges,
                                            # adjacent through the seam, not the interior
        a, b = x
        @test sup(a) < inf(b)                 # a ≈ [-π/2, -π/6], b ≈ [7π/6, 3π/2]
        # soundness spot-check: both arcs map back inside the (clamped) target
        @test issubset_interval(sin(a), intervaltype(-1.0, -0.5))
        @test issubset_interval(sin(b), intervaltype(-1.0, -0.5))
    end;

    @testset "Large ranges" begin
        y = intervaltype(nextfloat(-Inf), nextfloat(nextfloat(-Inf)))
        @test isnothing(sin_rev(y))
    end
end;

@testset "exp" begin
    y = intervaltype(1, 2)
    x = exp_rev(y)
    @test issubset_interval(exp(x), y)

    y = intervaltype(-1, 2)
    x = exp_rev(y)
    @test issubset_interval(exp(x), y)

    y = intervaltype(-2, -1)
    x = exp_rev(y)
    @test isnothing(x)
end;

@testset "log" begin
    y = intervaltype(1, 2)
    x = log_rev(y)
    @test issubset_interval(log(x), y)

    y = intervaltype(-100, 2)
    x = log_rev(y)
    @test issubset_interval(log(x), y)

    y = intervaltype(-Inf, 2)
    x = log_rev(y)
    @test issubset_interval(log(x), y)

    y = intervaltype(-Inf, Inf)
    x = log_rev(y)
    @test issubset_interval(log(x), y)

    y = intervaltype(3, 3)
    x = log_rev(y)
    @test isnothing(x)

    y = intervaltype(prevfloat(3.0), nextfloat(3.0))
    x = log_rev(y)
    @test issubset_interval(log(x), y)
end;

