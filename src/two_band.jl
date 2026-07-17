# Two-band tolerance ratchet — pure primitives (GitHub issue #2).
#
# The search maintains two symmetric tolerance bands around the noisy targets:
#   * the INNER (primary) band  — the tight band the search is actively pushing into,
#     carried as the `interval_targets` IntervalVector;
#   * the OUTER (secondary) band — a looser floor pinned at the last solved tolerance.
# Because both bands are symmetric around the same `noisy_targets` and τ_outer ≥ τ_inner,
# the inner interval nests inside the outer one, so a primary hit is always a secondary
# hit (the containment invariant). With τ_outer == τ_inner the two bands coincide and
# every two-band quantity below reduces exactly to its single-band counterpart.

using IntervalArithmetic: in_interval

# --- Module 1: two-band hit computation ------------------------------------------------
# Returns (primary, secondary) BitVectors. `primary` is membership in the inner band
# (`interval_targets`); `secondary` is membership in the symmetric outer band
# noisy ± τ_outer. Building the outer band the same way as the inner band (interval
# membership rather than a raw |eval − noisy| ≤ τ test) keeps the containment invariant
# exact under floating-point rounding: the inner interval nests inside the outer one.
function two_band_hits(evals::AbstractVector{<:Real},
                       interval_targets::IntervalVector,
                       noisy_targets::AbstractVector{<:Real},
                       tau_outer::Real)
    primary = hitvector(evals, interval_targets)
    n = length(interval_targets)
    secondary = BitVector(undef, n)
    for i in 1:n
        outer = intervaltype(noisy_targets[i] - tau_outer, noisy_targets[i] + tau_outer)
        secondary[i] = in_interval(evals[i], outer)
    end
    return primary, secondary
end

# --- Module 2: loss / score scalars ----------------------------------------------------
# loss = sum(primary_misses) + mean(secondary_misses)·(n−1)/n.
# The secondary term lies in [0, 1) so the integer primary-miss count always dominates
# ordering; under containment it is exactly 0 on a full inner solve, giving an
# n-invariant "0 = done" reported value.
function two_band_loss(primary::AbstractVector{Bool}, secondary::AbstractVector{Bool})
    n = length(primary)
    n == 0 && return 0.0
    primary_misses = count(!, primary)
    secondary_misses = count(!, secondary)
    return primary_misses + (secondary_misses / n) * (n - 1) / n
end

# Maximize-form sibling of `two_band_loss`: score = sum(primary_hits) +
# mean(secondary_hits)·(n−1)/n. For a fixed n, loss + score is constant, so the two
# induce opposite orderings — `two_band_loss` is preferred only for its n-invariant zero.
function two_band_score(primary::AbstractVector{Bool}, secondary::AbstractVector{Bool})
    n = length(primary)
    n == 0 && return 0.0
    primary_hits = count(identity, primary)
    secondary_hits = count(identity, secondary)
    return primary_hits + (secondary_hits / n) * (n - 1) / n
end

# --- Module 3: two-band lexicase -------------------------------------------------------
# `primary`/`secondary` are per-individual BitVectors indexed by case. For each case in
# `case_order`: filter the surviving pool to individuals that hit the inner band; if none
# do, fall back to those that hit the outer band; if neither, leave the pool unfiltered.
# Standard pool-collapse termination (stop once a single individual remains). With
# secondary ≡ primary this reduces to single-band lexicase (filter-or-skip).
function two_band_lexicase_pool(primary::AbstractVector{<:AbstractVector{Bool}},
                                secondary::AbstractVector{<:AbstractVector{Bool}},
                                case_order)
    remaining = collect(1:length(primary))
    for case in case_order
        length(remaining) <= 1 && break
        inner = [j for j in remaining if primary[j][case]]
        if !isempty(inner)
            remaining = inner
        else
            outer = [j for j in remaining if secondary[j][case]]
            isempty(outer) || (remaining = outer)   # else: no filter on this case
        end
    end
    return remaining
end

# Production wrapper: shuffle case order, filter, and break remaining ties at random.
function two_band_lexicase(primary::AbstractVector{<:AbstractVector{Bool}},
                           secondary::AbstractVector{<:AbstractVector{Bool}},
                           max_lexicase::Int, rng)
    m = length(first(primary))
    order = randperm(rng, m)
    ncases = min(max_lexicase, m)
    pool = two_band_lexicase_pool(primary, secondary, view(order, 1:ncases))
    return length(pool) == 1 ? pool[1] : rand(rng, pool)
end
