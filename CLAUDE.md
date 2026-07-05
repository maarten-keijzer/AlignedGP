# AlignedGP — development notes

## Invalid interval sentinel: `invalid_interval` = `CInterval(NaN, NaN)`

When an inverse function has no valid solution (e.g. `sqrt(x) < 0` is impossible),
it must return the **invalid sentinel** `invalid_interval` rather than `nothing` or
a misleading valid interval.

**Why not `nothing`?** Returning `nothing` forces every caller to check for it,
breaking the uniform interval arithmetic pipeline.

**Why not a clamped valid interval?** Returning e.g. `CI(0, 0)` implies `x = 0` is
a solution when it isn't. This silently inflates hit counts.

**Why `CInterval(NaN, NaN)`?** NaN propagates naturally through arithmetic, and any
interval with a NaN bound is unambiguously detectable via `_is_invalid(ci) =
isnan(ci.lo) || isnan(ci.hi)`. The constructor treats any NaN input as producing the
invalid sentinel (since `NaN <= NaN` is false, both fields end up NaN). `CI(Inf, Inf)`
and `CI(-Inf, -Inf)` are **valid** degenerate intervals and must not be treated as
sentinels.

## `CIntervals` filters invalid sentinels; `max_overlap_region` filters degenerate-infinite

Invalid sentinels produced at the `CInterval` level are filtered out by the `CIntervals`
lifting functions (`_lift_inverse`, `_scale`, `_div_into`), so `CIntervals` never
contains invalid sentinels. `max_overlap_region` applies two filters before the sweep:

- `_is_invalid`: removes any NaN-bearing intervals (defensive, should already be absent)
- `_is_useless(ci) = !isfinite(ci.lo) && ci.lo == ci.hi`: removes degenerate intervals
  like `CI(-Inf, -Inf)` produced by subtracting `Inf` evals from finite targets — these
  carry no information about a valid finite constant and would corrupt the depth count.

If all input intervals are filtered, returns `(; region=CIntervals(), depth=0)`.
Callers use `depth == 0` to detect that no valid constant exists.

## NaN signals undefined evaluation output

Node-level `evaluate` for unsafe functions returns `NaN` for inputs outside their
domain: `sqrt(x) = NaN` for `x < 0`, `log(x) = NaN` for `x ≤ 0`. This is the
canonical signal for "this data point produced no meaningful output."

**Why NaN (not Inf or 0)?**
- `Inf` is a valid output (e.g. `log` of a tiny positive number approaches `-Inf`;
  arithmetic can legitimately produce large values). Using `Inf` for "undefined" would
  conflate two distinct cases.
- `0` is a valid output and would silently inflate hit counts.
- `NaN` propagates through all arithmetic (`NaN op x = NaN`), so a single undefined
  node poisons the whole subtree evaluation for that point. `any(isnan, evals)` is a
  reliable, cheap check.

**Consequence for `compute_added_value`:** NaN evals produce `CInterval(NaN, NaN)` =
`invalid_interval` through normal CInterval subtraction arithmetic, which
`max_overlap_region` already filters via `_is_invalid`. NaN points simply don't vote
on the constant. The returned `hits` is the true optimal count over the finite points;
NaN evals stay NaN in `updated_evals` (NaN + c = NaN) and miss naturally.

**Consequence for `_optimize(BinaryNode)`:** the assertion is
`parent_hits >= new_right.addition.hits` (T3 ≥ T2, right-child Gauss-Seidel guarantee),
not `>= max_child`. See the `narrow` section below for why left.hits (T1) is not
directly comparable.

## `narrow` is an inner approximation — T3 ≥ T2 holds, T3 ≥ T1 does not

`narrow(lo, hi)` shrinks the interval by one ULP on each side, making surrogates
**inner approximations** of the true preimage. This guarantees: if
`right_eval ∈ narrow(lo, hi)`, then `left op right_eval ∈ target` after floating-point
rounding (one ULP of error is absorbed by the narrowing). Inner approximation is
correct in the forward direction: surrogate hit → parent hit.

**Why T3 ≥ T1 fails with inner approximation.** The left child's hit count
(`new_left.addition.hits`) is computed against left-surrogates built from
`right_evals_pre`. These count points where `left_eval - right_eval_pre ∈ target`
(T1 hits). For `_optimize(right)` to preserve all T1 hits, it needs `right_eval_pre`
to be inside every corresponding right-surrogate. But because right-surrogates are
INNER approximations, boundary values are excluded — a T1 hit where `right_eval_pre`
is at the exact boundary of `rightinverse(op, target, left_eval)` will not be in the
narrowed surrogate. That point is invisible to the right optimization.

Additionally, nested inverse computations (e.g. `inverse(exp, narrow(...))`) can
collapse intervals to `invalid_interval` at sub-ULP scales, further breaking T1
comparability. The cascading effect through deep trees amplifies this.

**`narrow` returning `invalid_interval` for sub-ULP intervals** (when `nextfloat(lo) >
prevfloat(hi)` after sorting) is correct and necessary: no float is strictly inside
such an interval, and allowing the `CInterval` constructor to silently swap
`(nextfloat(lo), prevfloat(lo))` into a spurious 3-float-wide interval produced false
surrogate hits that violated inner approximation (the original assertion failure with
`sqrt`).

**The assertion in `_optimize(BinaryNode)`** is therefore `parent_hits >=
new_right.addition.hits` (T3 ≥ T2). The right child is the last to be optimized, uses
the Gauss-Seidel context (updated left_evals), and its inner-approximated surrogates
directly bound parent hits.

## DomainErrors: no Tree-level catch

There is no `try/catch` at the Tree level. Unsafe node evaluates (`sqrt`, `log`)
return `NaN` explicitly for out-of-domain inputs, so no `DomainError` is thrown.
Do not rely on a Tree-level safety net — handle undefined outputs at the node level.

## Intervals.jl cannot represent empty bounded-closed intervals

`Interval{T,Closed,Closed}` validates and sorts its constructor arguments, so a
traditional "empty by construction" trick (`CI(high, low)`) does not work. Use
`invalid_interval` for unreachable cases instead.

## Two distinct size measures: `length` vs `complexity`

Nodes carry two separate size measures that must not be conflated.

**`Base.length(node)`** — structural node count. Counts only indexable nodes
(`BinaryNode`, `UnaryNode`, `Var`, `Constant`). Additive constants (`AddedValue`)
never contribute. This guarantees `node[rand(1:length(node))]` always returns a
valid subtree. Used by `getindex`, `insert`, `insert_with_alignment`, and `targetsize`
in `init`.

**`complexity(node)`** — effective expression size. Adds 1 for each non-zero
`AddedValue` on top of structural length. Used by `SizeStratum` and any size/fit
balance logic. `BinaryNode` and `UnaryNode` cache both `size` (structural) and
`complexity` as `Int` fields; `Tree` also caches `complexity`.

**`Constant` is special:** its `addition.value` IS the constant's value, not an
offset on top of some other computation. Therefore `complexity(::Constant) = 1`
always, even when `addition.value == 0`. A zero `Constant` is still a concrete
numeric node in the expression.

**`complexity(::AddedValue)`** encodes the conditional: `value == 0 ? 0 : 1`. This
is the only place that tests whether an additive constant contributes to complexity.
Do not use `Base.length` on `AddedValue`.
