using Test
using AlignedGP

import AlignedGP: compute_added_value, align_node

const CI = CInterval

# Per-point hit vector: true where evals[j] lands inside targets[j]
_hits(evals, targets) = [evals[j] ∈ targets[j] for j in eachindex(targets)]

# ──────────────────────────────────────────────────────────────────────────────
@testset "compute_added_value" begin

    @testset "constant shift needed, all hits" begin
        evals   = [1.0, 2.0, 3.0]
        targets = CIntervals.([CI(5.0, 7.0), CI(6.0, 8.0), CI(7.0, 9.0)])
        av, updated = compute_added_value(evals, targets)
        # shifted intervals all CI(4,6); zero not in region → sample within [4, 6]
        @test 4.0 <= av.value <= 6.0
        @test updated ≈ evals .+ av.value
        @test _hits(updated, targets) == [true, true, true]
    end

    @testset "zero is the optimal constant" begin
        evals   = [3.0, 4.0, 5.0]
        targets = CIntervals.([CI(2.0, 5.0), CI(3.0, 6.0), CI(4.0, 7.0)])
        av, updated = compute_added_value(evals, targets)
        # shifted intervals all CI(-1,2); zero lies in that region
        @test av.value == 0.0
        @test updated ≈ [3.0, 4.0, 5.0]
        @test _hits(updated, targets) == [true, true, true]
    end

    @testset "partial hits — oracle hit vector" begin
        evals   = [0.0, 5.0, 10.0]
        targets = CIntervals.([CI(1.0, 3.0), CI(6.0, 8.0), CI(20.0, 22.0)])
        av, updated = compute_added_value(evals, targets)
        # shifted: CI(1,3), CI(1,3), CI(10,12) → depth 2 at CI(1,3), sample within [1,3]
        @test 1.0 <= av.value <= 3.0
        @test _hits(updated, targets) == [true, true, false]
    end

    @testset "Inf in evals → graceful miss for that point" begin
        evals   = [1.0, Inf, 3.0]
        targets = CIntervals.([CI(2.0, 4.0), CI(5.0, 7.0), CI(4.0, 6.0)])
        av, updated = compute_added_value(evals, targets)
        # point 2: CI(5-Inf, 7-Inf) = CI(-Inf,-Inf) is filtered as invalid sentinel
        # remaining depth-2 region from points 1 & 3: CI(1,3) → sample within [1,3]
        @test 1.0 <= av.value <= 3.0
        @test isinf(updated[2])
        @test _hits(updated, targets) == [true, false, true]
    end

end

# ──────────────────────────────────────────────────────────────────────────────
@testset "align_node" begin

    inputs  = [[1.0, 2.0, 3.0]]
    targets = CIntervals.([CI(5.0, 7.0), CI(6.0, 8.0), CI(7.0, 9.0)])

    @testset "fresh Var node" begin
        node, new_outputs = align_node(Var(1), targets, inputs)
        @test node isa Var
        @test 4.0 <= node.addition.value <= 6.0   # overlap of shifted targets is [4,6]
        @test new_outputs ≈ evaluate(node, inputs)
        @test _hits(new_outputs, targets) == [true, true, true]
    end

    @testset "existing AddedValue is stripped before refitting" begin
        # Var(1) + 10 should produce the same addition as a fresh Var(1), not +15
        node_old = Var(1, AddedValue(10.0))
        node, new_outputs = align_node(node_old, targets, inputs)
        @test 4.0 <= node.addition.value <= 6.0   # refit from scratch, not old+shift
        @test node.addition.value ≉ 15.0
        @test new_outputs ≈ evaluate(node, inputs)
        @test _hits(new_outputs, targets) == [true, true, true]
    end

    @testset "donation with large AddedValue is also stripped correctly" begin
        donated = Var(1, AddedValue(100.0))
        node, new_outputs = align_node(donated, targets, inputs)
        @test 4.0 <= node.addition.value <= 6.0   # large donated AddedValue is stripped
        @test new_outputs ≈ evaluate(node, inputs)
        @test _hits(new_outputs, targets) == [true, true, true]
    end

    @testset "node evaluating to Inf → graceful miss, no error" begin
        inf_inputs = [[Inf, 2.0, 3.0]]
        targets2   = CIntervals.([CI(4.0, 6.0), CI(4.0, 6.0), CI(4.0, 6.0)])
        node, new_outputs = align_node(Var(1), targets2, inf_inputs)
        @test new_outputs ≈ evaluate(node, inf_inputs)
        h = _hits(new_outputs, targets2)
        @test h[1] == false     # Inf cannot land in a finite band
        @test h[2] == true
        @test h[3] == true
    end

end
