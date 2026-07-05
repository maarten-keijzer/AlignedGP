using Test
using AlignedGP
using Statistics

@testset "linear_scale_l1" begin

    @testset "exact affine fit converges in few iterations" begin
        x = Float64.(1:20)
        t = 3.0 .* x .+ 7.0
        slope, intercept, mae, iters = linear_scale_l1(x, t)
        @test slope     ≈ 3.0  atol=1e-6
        @test intercept ≈ 7.0  atol=1e-6
        @test mae       ≈ 0.0  atol=1e-6
        @test iters     >= 1
    end

    @testset "returns mae not mse" begin
        x = [0.0, 1.0, 2.0, 3.0]
        t = [0.0, 1.0, 2.0, 10.0]   # one large outlier
        slope, intercept, mae, _ = linear_scale_l1(x, t)
        residuals = abs.(slope .* x .+ intercept .- t)
        @test mae ≈ mean(residuals) atol=1e-10
    end

    @testset "robust to outlier: slope closer to true than OLS" begin
        # true line is y = x; one large outlier at the last point
        x = Float64.(1:10)
        t = copy(x); t[10] = 100.0
        slope_l1, _, _, _ = linear_scale_l1(x, t)
        slope_ols, _, _   = linear_scale(x, t)
        @test abs(slope_l1 - 1.0) < abs(slope_ols - 1.0)
    end

    @testset "degenerate constant output uses median intercept" begin
        x = fill(5.0, 10)
        t = [1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 100.0]
        slope, intercept, mae, iters = linear_scale_l1(x, t)
        @test slope     == 0.0
        @test intercept ≈ median(t)
        @test mae       ≈ mean(abs.(t .- median(t))) atol=1e-10
        @test iters     == 0
    end

    @testset "iteration count within bounds" begin
        x = randn(50)
        t = 2.0 .* x .+ 1.0 .+ 0.1 .* randn(50)
        _, _, _, iters = linear_scale_l1(x, t)
        @test 1 <= iters <= 50
    end

    @testset "result length is 4" begin
        result = linear_scale_l1([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
        @test length(result) == 4
    end

end
