using Random
const CI = CInterval

@testset "circular_compute_added_value — minimal failing case" begin
    # From specs/sine-inversion.md: target {1.0}, current arguments {π/2, π/2+2π, 0}.
    # Optimal c=0: sin(π/2+0)=1 ✓, sin(π/2+2π+0)=1 ✓, sin(0+0)=0 ✗  → 2 hits.
    # Linear stab would miss because shifts land in different periods on the real line.
    rng = MersenneTwister(42)
    raw_evals = [π/2, π/2 + 2π, 0.0]
    targets = [CIntervals(CI(1.0, 1.0)), CIntervals(CI(1.0, 1.0)), CIntervals(CI(1.0, 1.0))]

    av, new_evals = AlignedGP.circular_compute_added_value(raw_evals, targets, sin, rng)
    @test av.hits == 2
    @test av.value ≈ 0.0 atol=1e-10
end

@testset "circular_compute_added_value — all points in same period" begin
    # When all points are in [0,2π), circular and linear stab should agree.
    rng = MersenneTwister(1)
    raw_evals = [0.1, 0.5, 1.0, 2.0]
    target = CI(0.0, 1.0)   # sin hits for x around 0 to π
    targets = fill(CIntervals(target), length(raw_evals))

    av, _ = AlignedGP.circular_compute_added_value(raw_evals, targets, sin, rng)
    # Verify hit count is achievable
    actual = sum(sin(raw_evals[i] + av.value) ∈ target for i in eachindex(raw_evals))
    @test actual == av.hits
end

@testset "_optimize with sin node over multi-period inputs" begin
    # Node: sin(x), with x = asin(0.7) + k*2π for k=0..3 (4 points spread over 4 periods).
    # All 4 points hit target [0.6,0.8] at c=0.  Linear stab misses 3 of the 4 points
    # because their shifted preimage lands at ~[-k*2π, -k*2π+ε], disjoint from c=0.
    rng = MersenneTwister(7)
    target = CI(0.6, 0.8)
    xs = [asin(0.7) + k*2π for k in 0:3]
    inputs = [xs]
    targets = fill(CIntervals(target), length(xs))

    root = UnaryNode(sin, Var(1))
    optimized, updated_evals = AlignedGP.optimize(root, inputs, targets, rng)

    c = optimized.child.addition.value
    @test c ≈ 0.0 atol=1e-8

    # All 4 points should be hits.
    actual_hits = sum(sin(xs[i] + c) ∈ target for i in eachindex(xs))
    @test actual_hits == 4
    @test optimized.addition.hits >= 4
end
