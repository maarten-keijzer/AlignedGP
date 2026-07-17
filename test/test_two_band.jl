using Test
using AlignedGP
using AlignedGP.ReverseIntervals
using Random

# ---------------------------------------------------------------------------
# Module 1 — two_band_hits: (evals, inner interval_targets, noisy, τ_outer)
#            → (primary::BitVector, secondary::BitVector)
# Primary = interval membership against the inner band.
# Secondary = membership against the symmetric outer band noisy ± τ_outer.
# ---------------------------------------------------------------------------
@testset "two_band_hits" begin

    # inner band is ±0.1 around noisy, outer band is ±0.5.
    noisy = zeros(5)
    tau_inner = 0.1
    tau_outer = 0.5
    inner = IntervalVector(intervaltype.(noisy .- tau_inner, noisy .+ tau_inner))

    #        inside inner | outer only | outside both | on inner edge | on outer edge
    evals = [0.05,          0.3,         1.0,           0.1,            0.5]

    @testset "primary and secondary per case" begin
        primary, secondary = AlignedGP.two_band_hits(evals, inner, noisy, tau_outer)
        @test primary   == BitVector([1, 0, 0, 1, 0])
        @test secondary == BitVector([1, 1, 0, 1, 1])
    end

    @testset "containment invariant: primary ⟹ secondary" begin
        primary, secondary = AlignedGP.two_band_hits(evals, inner, noisy, tau_outer)
        @test all(secondary[i] for i in eachindex(primary) if primary[i])
    end

    @testset "on the inner edge counts as a primary (closed) hit" begin
        primary, _ = AlignedGP.two_band_hits(evals, inner, noisy, tau_outer)
        @test primary[4]                       # eval == inner hi bound
    end

    @testset "inside outer only is a secondary-not-primary hit" begin
        primary, secondary = AlignedGP.two_band_hits(evals, inner, noisy, tau_outer)
        @test !primary[2] && secondary[2]
    end

    @testset "degenerate τ_outer == τ_inner ⇒ identical bitvectors" begin
        primary, secondary = AlignedGP.two_band_hits(evals, inner, noisy, tau_inner)
        @test primary == secondary
    end

    @testset "containment holds on random symmetric bands" begin
        rng = MersenneTwister(7)
        for _ in 1:50
            n = rand(rng, 1:12)
            nz = randn(rng, n)
            ti = 0.05 + rand(rng)
            to = ti + rand(rng)                 # outer ≥ inner
            iv = IntervalVector(intervaltype.(nz .- ti, nz .+ ti))
            ev = nz .+ (2rand(rng, n) .- 1) .* (to + 0.3)
            p, s = AlignedGP.two_band_hits(ev, iv, nz, to)
            @test all(s[i] for i in eachindex(p) if p[i])
        end
    end
end

# ---------------------------------------------------------------------------
# Module 2 — two_band_loss(primary, secondary)
#   loss = sum(primary_misses) + mean(secondary_misses)·(n−1)/n
#   secondary term ∈ [0,1); loss == 0 exactly on a full inner solve.
# ---------------------------------------------------------------------------
@testset "two_band_loss" begin

    @testset "exact value" begin
        p = BitVector([1, 1, 1, 0])
        s = BitVector([1, 1, 1, 0])
        @test AlignedGP.two_band_loss(p, s) ≈ 1 + (1 / 4) * (3 / 4)   # 1.1875
    end

    @testset "loss == 0 iff zero primary misses (n-invariant)" begin
        for n in (1, 3, 10)
            allhit = trues(n)
            @test AlignedGP.two_band_loss(allhit, allhit) == 0.0
        end
        # one primary miss ⇒ strictly positive regardless of secondary
        @test AlignedGP.two_band_loss(BitVector([1, 0, 1]), BitVector([1, 1, 1])) > 0
    end

    @testset "secondary term always in [0, 1)" begin
        rng = MersenneTwister(11)
        for _ in 1:200
            n = rand(rng, 1:15)
            p = bitrand(rng, n)
            s = p .| bitrand(rng, n)             # containment: primary ⟹ secondary
            loss = AlignedGP.two_band_loss(p, s)
            frac = loss - count(!, p)            # subtract integer primary-miss part
            @test 0.0 <= frac < 1.0
        end
    end

    @testset "primary misses dominate the ordering" begin
        # A has more primary hits but the worst possible secondary;
        # B has fewer primary hits but perfect secondary. A must still win.
        A_p = BitVector([1, 1, 0]); A_s = BitVector([0, 0, 0])
        B_p = BitVector([1, 0, 0]); B_s = BitVector([1, 1, 1])
        @test AlignedGP.two_band_loss(A_p, A_s) < AlignedGP.two_band_loss(B_p, B_s)
    end

    @testset "among equal primary, fewer secondary misses wins" begin
        p = BitVector([1, 1, 0])
        A_s = BitVector([1, 1, 0])               # 1 secondary miss
        B_s = BitVector([1, 1, 1])               # 0 secondary misses
        @test AlignedGP.two_band_loss(p, B_s) < AlignedGP.two_band_loss(p, A_s)
    end

    @testset "ordering matches the maximize (score) form" begin
        rng = MersenneTwister(3)
        for _ in 1:300
            n = rand(rng, 1:15)
            pa = bitrand(rng, n); sa = pa .| bitrand(rng, n)
            pb = bitrand(rng, n); sb = pb .| bitrand(rng, n)
            la = AlignedGP.two_band_loss(pa, sa); lb = AlignedGP.two_band_loss(pb, sb)
            ga = AlignedGP.two_band_score(pa, sa); gb = AlignedGP.two_band_score(pb, sb)
            @test sign(la - lb) == -sign(ga - gb)
        end
    end
end

# ---------------------------------------------------------------------------
# Module 3 — two_band lexicase fallback.
#   two_band_lexicase_pool(primary, secondary, case_order) filters, per case:
#   to inner-band hitters; else outer-band hitters; else no filter.
# ---------------------------------------------------------------------------
@testset "two_band_lexicase" begin

    # independent single-band reference (filter-or-skip) for the degenerate case
    function single_band_pool(primary, case_order)
        remaining = collect(1:length(primary))
        for case in case_order
            length(remaining) <= 1 && break
            hitters = [j for j in remaining if primary[j][case]]
            isempty(hitters) || (remaining = hitters)
        end
        return remaining
    end

    @testset "a case with an inner hitter filters to inner hitters" begin
        primary   = [BitVector([1]), BitVector([0]), BitVector([0])]
        secondary = [BitVector([1]), BitVector([1]), BitVector([0])]
        @test AlignedGP.two_band_lexicase_pool(primary, secondary, [1]) == [1]
    end

    @testset "no inner hitter falls back to outer hitters" begin
        primary   = [BitVector([0]), BitVector([0]), BitVector([0])]
        secondary = [BitVector([1]), BitVector([1]), BitVector([0])]
        pool = AlignedGP.two_band_lexicase_pool(primary, secondary, [1])
        @test sort(pool) == [1, 2]               # ind 3 (misses outer) dropped
    end

    @testset "no hitter of either band does not filter" begin
        primary   = [BitVector([0]), BitVector([0]), BitVector([0])]
        secondary = [BitVector([0]), BitVector([0]), BitVector([0])]
        pool = AlignedGP.two_band_lexicase_pool(primary, secondary, [1])
        @test sort(pool) == [1, 2, 3]            # pool untouched
    end

    @testset "multi-case: inner pressure then outer fallback" begin
        # case 1: only ind1 hits inner → collapses to {1}
        # (ind1 alone, so later cases can't change it)
        primary   = [BitVector([1, 0]), BitVector([0, 0]), BitVector([0, 0])]
        secondary = [BitVector([1, 1]), BitVector([0, 1]), BitVector([0, 0])]
        @test AlignedGP.two_band_lexicase_pool(primary, secondary, [1, 2]) == [1]

        # case 1: nobody hits inner, {1,2} hit outer → {1,2};
        # case 2: ind1 hits inner among the survivors → {1}
        primary2   = [BitVector([0, 1]), BitVector([0, 0]), BitVector([0, 0])]
        secondary2 = [BitVector([1, 1]), BitVector([1, 1]), BitVector([0, 0])]
        @test AlignedGP.two_band_lexicase_pool(primary2, secondary2, [1, 2]) == [1]
    end

    @testset "secondary ≡ primary reduces to single-band lexicase" begin
        rng = MersenneTwister(99)
        for _ in 1:100
            npop = rand(rng, 2:6)
            m = rand(rng, 1:6)
            primary = [bitrand(rng, m) for _ in 1:npop]
            order = randperm(rng, m)
            @test AlignedGP.two_band_lexicase_pool(primary, primary, order) ==
                  single_band_pool(primary, order)
        end
    end

    @testset "two_band_lexicase returns a valid pool member" begin
        rng = MersenneTwister(5)
        primary   = [BitVector([1, 0]), BitVector([0, 1]), BitVector([0, 0])]
        secondary = [BitVector([1, 1]), BitVector([0, 1]), BitVector([1, 0])]
        for _ in 1:20
            sel = AlignedGP.two_band_lexicase(primary, secondary, 2, rng)
            @test sel in 1:3
        end
    end
end

# ---------------------------------------------------------------------------
# Integration — Tree fields, retarget, and the ratchet band transition through
# a real ProblemSetup (no re-evaluation on retarget).
# ---------------------------------------------------------------------------
@testset "two-band Tree integration" begin
    setup = keijzer1(tol=0.5, noise=0.0, rng=MersenneTwister(2024))
    n = length(setup.interval_targets)
    node = BinaryNode(+, Var(1), Constant(0.3))
    tree = AlignedGP.evaluate_to_tree(node, setup)

    @testset "fields are populated and self-consistent" begin
        @test length(tree.evals) == n
        @test length(tree.secondary_hits) == n
        @test tree.loss ≈ two_band_loss(tree.hits, tree.secondary_hits)
        @test all(tree.secondary_hits[i] for i in 1:n if tree.hits[i])   # containment
    end

    @testset "single-band start (τ_outer == tol): secondary ≡ primary" begin
        @test setup.params.tau_outer == 0.5
        @test tree.hits == tree.secondary_hits
    end

    @testset "ratchet: retarget tightens inner band without re-evaluating" begin
        setup.params.tau_outer = 0.5                       # freeze floor
        new_tol = 0.1
        newint = IntervalVector(intervaltype.(setup.noisy_targets .- new_tol,
                                              setup.noisy_targets .+ new_tol))
        setup.interval_targets.intervals .= newint.intervals

        rt = retarget(tree, setup)
        @test isequal(rt.evals, tree.evals)                # evals reused, not recomputed
        @test rt.slope == tree.slope && rt.intercept == tree.intercept && rt.mse == tree.mse
        @test all(rt.secondary_hits[i] for i in 1:n if rt.hits[i])   # containment holds
        @test sum(rt.hits) <= sum(rt.secondary_hits)       # tighter inner ⊆ outer floor
        @test rt.loss ≈ two_band_loss(rt.hits, rt.secondary_hits)

        # retarget agrees with a full re-evaluation against the same bands
        fresh = AlignedGP.evaluate_to_tree(tree.root, setup)
        @test rt.hits == fresh.hits
        @test rt.secondary_hits == fresh.secondary_hits
        @test rt.loss ≈ fresh.loss
    end
end
