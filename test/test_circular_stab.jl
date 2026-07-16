using Test
using AlignedGP
using AlignedGP.ReverseIntervals
using Random

using IntervalArithmetic: in_interval

const C = 2π

# sin_rev arcs (u-space) for a per-case band target, as an IntervalVector
sinrev_targets(band, n) = invert(IntervalVector(fill(band, n)), sin_rev)

# ---------------------------------------------------------------------------
# fold_stab — the circular stab for the constant directly below sin (§3.3)
# ---------------------------------------------------------------------------
@testset "fold_stab" begin

    @testset "empty" begin
        r = AlignedGP.fold_stab(IntervalType[])
        @test r.depth == 0
        @test isempty(r.region)
    end

    @testset "wrapping region found whole via the -2π copy" begin
        # Two arcs that both sit just below 2π and overlap near the seam.
        # A plain linear sweep on the un-duplicated arcs would not see them
        # stack across 0; the -C copy makes the overlap visible.
        r = AlignedGP.fold_stab([intervaltype(6.0, 6.4), intervaltype(6.1, 6.5)])
        @test r.depth == 2
    end

    @testset "two cases 2π apart stack (multi-period agreement, §3.4)" begin
        # Canonical arcs from two different cases that are 2π apart are the *same*
        # constant mod 2π, so they must stack — this is the whole point of folding.
        r = AlignedGP.fold_stab([intervaltype(0.1, 0.3), intervaltype(0.1 + C, 0.3 + C)])
        @test r.depth == 2
    end

    @testset "a single arc never self-stacks (its two copies are C apart)" begin
        # One case's arc is emitted twice, C apart; disjoint since width < C, so no
        # point is ever counted twice for one case (disjointness invariant, §1).
        r = AlignedGP.fold_stab([intervaltype(0.1, 0.3)])
        @test r.depth == 1
    end

    @testset "nfull folded into depth" begin
        # One full-circle arc (hits every constant) + one ordinary arc.
        r = AlignedGP.fold_stab([intervaltype(0.0, C), intervaltype(1.0, 1.2)])
        @test r.depth == 2                       # 1 (full) + 1 (ordinary) where they overlap
    end

    @testset "all-full → empty region, depth = nfull" begin
        r = AlignedGP.fold_stab([intervaltype(0.0, C), intervaltype(0.0, C)])
        @test r.depth == 2
        @test isempty(r.region)                  # every constant equally good ⇒ select_constant → 0
    end
end

@testset "fold_stab per-case disjointness (IntervalVector)" begin

    @testset "two touching arcs in one case count once" begin
        # A single case whose two arcs meet at ≈0 — as arises when a sin node sits
        # below deeper structure and sin_rev maps several incoming arcs, two of which
        # touch after inner-rounding. Pooling would count the shared point twice;
        # the case must contribute depth 1 (disjointness invariant, §1 / open-Q #4).
        carcs = IntervalVector([intervaltype(-0.02, 0.0), intervaltype(0.0, 0.03)], [1, 3])
        r = AlignedGP.fold_stab(carcs)
        @test r.depth == 1
    end

    @testset "genuinely disjoint arcs in one case still count once" begin
        # rising + falling of one sin_rev case: disjoint mod 2π, so at any c at most
        # one covers ⇒ the case contributes 1, never 2.
        carcs = IntervalVector([intervaltype(0.1, 0.4), intervaltype(2.7, 3.0)], [1, 3])
        r = AlignedGP.fold_stab(carcs)
        @test r.depth == 1
    end

    @testset "two different cases agreeing still stack" begin
        # case 1 = {[0.0,0.2]}, case 2 = {[0.1,0.3]} overlap at [0.1,0.2] ⇒ depth 2.
        carcs = IntervalVector([intervaltype(0.0, 0.2), intervaltype(0.1, 0.3)], [1, 2, 3])
        r = AlignedGP.fold_stab(carcs)
        @test r.depth == 2
    end
end

# ---------------------------------------------------------------------------
# circular_hits — truthful hit count mod 2π, must agree with fold depth
# ---------------------------------------------------------------------------
@testset "circular_hits" begin

    @testset "agrees with fold_stab depth at the chosen constant" begin
        evals   = [π/2, π/2 + C, 0.0]                 # spread over periods
        targets = sinrev_targets(intervaltype(0.9, 1.0), 3)
        carcs   = targets - evals
        res, depth = AlignedGP.fold_stab(carcs.intervals)
        value = select_constant(res)
        @test AlignedGP.circular_hits(carcs, value) == depth
    end

    @testset "full-circle arc always hits" begin
        # target ⊇ [-1,1] ⇒ sin_rev returns the full circle ⇒ every point hits
        carcs = sinrev_targets(intervaltype(-2.0, 2.0), 2) - [123.4, -987.6]
        @test AlignedGP.circular_hits(carcs, 0.37) == 2
    end
end

# ---------------------------------------------------------------------------
# compute_added_value(...; circular=true) — the folded additive constant
# ---------------------------------------------------------------------------
@testset "circular compute_added_value" begin

    @testset "multi-period agreement picks the shared constant" begin
        # From specs/sine-inversion.md: arguments span several periods but agree that
        # c=0 places cases 1 & 2 at the crest. Case 3 (sin(0)=0) misses. Linear stab
        # would miss cases in other periods; the fold sees them stack mod 2π.
        rng     = MersenneTwister(42)
        evals   = [π/2, π/2 + C, 0.0]
        band    = intervaltype(0.9, 1.0)
        targets = sinrev_targets(band, 3)

        av, _ = AlignedGP.compute_added_value(evals, targets, rng; circular=true)

        @test av.value ≈ 0.0 atol = 1e-8
        @test av.hits == 2
        # hits are truthful against the real forward model
        fwd = sum(in_interval(sin(evals[i] + av.value), band) for i in eachindex(evals))
        @test fwd == av.hits
    end

    @testset "value canonicalised into (-π, π]" begin
        rng     = MersenneTwister(3)
        evals   = [π/2, π/2 + C, 0.0]
        targets = sinrev_targets(intervaltype(0.9, 1.0), 3)
        av, _ = AlignedGP.compute_added_value(evals, targets, rng; circular=true)
        @test -π < av.value <= π
    end

    @testset "same period ⇒ hits stay truthful" begin
        rng     = MersenneTwister(1)
        evals   = [0.1, 0.5, 1.0, 2.0]               # all within one period
        targets = sinrev_targets(intervaltype(0.6, 0.8), length(evals))
        av, _ = AlignedGP.compute_added_value(copy(evals), targets, rng; circular=true)
        fwd = sum(in_interval(sin(evals[i] + av.value), intervaltype(0.6, 0.8)) for i in eachindex(evals))
        @test fwd == av.hits
        @test av.hits >= 1
    end
end

# ---------------------------------------------------------------------------
# integration through insert_with_alignment — the production path
# ---------------------------------------------------------------------------
@testset "insert_with_alignment folds the constant below sin" begin

    @testset "four points over four periods all become hits" begin
        # x = asin(0.7) + k·2π, k=0..3.  At c=0 every point hits [0.6,0.8].
        # Linear stab would keep only one; the circular fold keeps all four.
        band    = intervaltype(0.6, 0.8)
        xs      = [asin(0.7) + k*C for k in 0:3]
        inputs  = [xs]
        targets = IntervalVector(fill(band, length(xs)))

        root = UnaryNode(sin, Var(1))
        # Insert Var(1) at the child slot (index 2); recursive_stabbing = true.
        new_root, _, _ = insert_with_alignment(root, Var(1), 2, 1, inputs, targets, true)

        c = new_root.child.addition.value
        @test -π < c <= π
        forward_hits = sum(in_interval(sin(xs[i] + c), band) for i in eachindex(xs))
        @test forward_hits == 4
        @test new_root.child.addition.hits == 4
    end

    @testset "assertion holds on a sin tree over multi-period data (no throw)" begin
        rng     = MersenneTwister(11)
        band    = intervaltype(0.2, 0.9)
        xs      = [0.3 + k*C + 0.7*randn(rng) for k in 0:9]
        inputs  = [xs]
        targets = IntervalVector(fill(band, length(xs)))
        root = UnaryNode(sin, Var(1))
        # Must not trip the `hits >= depth` assertion in compute_added_value.
        @test_nowarn insert_with_alignment(root, Var(1), 2, 1, inputs, targets, true)
    end
end
