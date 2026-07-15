# AlignedGP — development notes

## Interval representation: `ReverseIntervals` over `IntervalArithmetic.BareInterval`

Intervals live in the `ReverseIntervals` module (`src/Intervals/`). The concrete type is
`IntervalType = BareInterval{Float64}`, built with `intervaltype(x)` (a point) or
`intervaltype(lo, hi)`. There is no bespoke `CInterval` type anymore — the earlier
hand-rolled `CInterval`/`CIntervals`/`IntervalSets.jl` machinery has been removed.

**`IntervalVector`** is a ragged (flattened) array of interval slices, one *case* per
data point:

- `intervals::Vector{IntervalType}` — all sub-intervals, concatenated.
- `offsets::Vector{Int}` — case `i` spans `offsets[i]:(offsets[i+1]-1)`.

`getindex(iv, i)` returns a **view** into that slice (no allocation). A case can hold
**zero, one, or several disjoint** intervals: several because a preimage can split
(e.g. `inv_rev` and `sin_rev` return tuples of arcs), and **zero** when the case has no
valid preimage. An empty slice therefore means "this data point is unreachable under the
current surrogate target."

## `nothing` signals "no valid preimage" (there is no invalid sentinel)

When an inverse function has no valid solution (e.g. `sqrt(x) < 0` is impossible), the
reverse functions in `rev_functions.jl` return **`nothing`** (or, for the point/degenerate
cases, `maybeinterval` returns `nothing`). The `invert` lifting functions in `Intervals.jl`
**skip `nothing` results** — they simply are not pushed into `intervals`, so the case ends
up with an empty slice. Consequences:

- An `IntervalVector` never contains a NaN-bearing or "invalid" interval; unreachable cases
  are represented by *absence* (empty slice), not a sentinel value.
- Callers detect "no valid constant" via an empty region / `depth == 0` from
  `max_overlap_region`, not by inspecting interval bounds.

`maybeinterval(lo, hi) = (lo != Inf && lo <= hi) ? intervaltype(lo, hi) : nothing` is the
single gate that collapses empty/unreachable results to `nothing`. This replaces the old
`invalid_interval = CInterval(NaN, NaN)` design: we no longer need NaN to propagate a
sentinel because `nothing` is dropped at the lifting boundary instead.

## Reverse functions are inner approximations via `nextfloat`/`prevfloat`

Each reverse function narrows its computed bounds inward by one ULP (`nextfloat` on the
low side, `prevfloat` on the high side) — e.g. `sqrt_rev` uses `xl = nextfloat(yl*yl)`,
`xh = prevfloat(yh*yh)`; `mul_rev`, `inv_rev`, `exp_rev`, `log_rev` do the same. This makes
surrogate targets **inner approximations** of the true preimage, which is correct in the
forward direction: if `arg ∈ rev_interval`, then `f(arg) ∈ target` even after
floating-point rounding (the one-ULP margin absorbs the rounding error). Surrogate hit ⇒
parent hit. There is no separate `narrow()` helper — narrowing is baked into each reverse
function, and `maybeinterval` returns `nothing` when the inward narrowing makes `lo > hi`.

## `max_overlap_region` — depth sweep over a flat interval vector

`max_overlap_region(intervals::Vector{IntervalType})` performs an enter/leave sweep and
returns `(; region::Vector{IntervalType}, depth)`: the maximum overlap depth `k` and every
closed range covered by exactly `k` intervals. It does **no** pre-filtering — because
`invert` already dropped unreachable/non-finite cases upstream, the input carries no invalid
sentinels. `depth == 0` (empty input) means no valid constant exists; `select_constant`
returns `0.0` for an empty region. Half-infinite intervals (from `mul_rev` with `y == 0`,
`inv_rev`, `exp_rev`) can legitimately appear, and `select_constant` handles them via its
finite-bound fallback.

## NaN signals undefined evaluation output

Node-level `evaluate` for unsafe functions returns `NaN` for inputs outside their domain:
`sqrt(x) = NaN` for `x < 0`, `log(x) = NaN` for `x ≤ 0`, and `sin`/`cos` return `NaN` for
non-finite inputs. This is the canonical signal for "this data point produced no meaningful
output."

**Why NaN (not Inf or 0)?**
- `Inf` is a valid output (e.g. `log` of a tiny positive number approaches `-Inf`;
  arithmetic can legitimately produce large values). Using `Inf` for "undefined" would
  conflate two distinct cases.
- `0` is a valid output and would silently inflate hit counts.
- `NaN` propagates through all arithmetic (`NaN op x = NaN`), so a single undefined node
  poisons the whole subtree evaluation for that point. `any(isnan, evals)` is a reliable,
  cheap check.

**Consequence for `compute_added_value`:** it forms surrogate targets via
`targets - evals` (which lifts through `add_rev`). For a NaN eval, `add_rev(z, NaN)` sees
`!isfinite(NaN)` and returns `nothing`, so `invert` drops that case → empty slice → the
point contributes no interval and does not vote on the constant. The returned `hits` is the
true optimal count over the finite points; NaN evals stay NaN in `updated_evals`
(`NaN + c = NaN`) and miss naturally.

**Consequence for `_optimize(BinaryNode)`:** the assertion is
`parent_hits >= new_right.addition.hits` (T3 ≥ T2, right-child Gauss-Seidel guarantee),
not `>= max_child`. See the inner-approximation section below for why left.hits (T1) is not
directly comparable.

## Inner approximation ⇒ T3 ≥ T2 holds, T3 ≥ T1 does not

Because reverse functions produce **inner approximations** (one-ULP inward narrowing),
`_optimize` can assert the right-child bound but not the left-child bound:

**Why T3 ≥ T1 fails.** The left child's hit count (`new_left.addition.hits`) is computed
against left-surrogates built from `right_evals_pre` (T1 hits). For `_optimize(right)` to
preserve all T1 hits, `right_eval_pre` would need to be inside every corresponding
right-surrogate — but inner approximations exclude boundary values, so a T1 hit where
`right_eval_pre` sits at the exact boundary of the right-inverse is invisible to the right
optimization. Nested inverse computations (e.g. `inverse(exp, …)`) can also collapse
sub-ULP intervals to `nothing`, further breaking T1 comparability, and the effect cascades
through deep trees.

**The assertion in `_optimize(BinaryNode)`** is therefore `parent_hits >=
new_right.addition.hits` (T3 ≥ T2). The right child is optimized last, uses the
Gauss-Seidel context (updated `left_evals`), and its inner-approximated surrogates directly
bound the parent's hits.

## DomainErrors: no Tree-level catch

There is no `try/catch` at the Tree level. Unsafe node evaluates (`sqrt`, `log`, `sin`,
`cos`) return `NaN` explicitly for out-of-domain / non-finite inputs, so no `DomainError`
is thrown. Do not rely on a Tree-level safety net — handle undefined outputs at the node
level.

## Three size measures: `length`, `complexity`, `pathlen_complexity`

Nodes carry three separate size measures that must not be conflated.

**`Base.length(node)`** — structural node count. Counts only indexable nodes (`BinaryNode`,
`UnaryNode`, `Var`, `Constant`). Additive constants (`AddedValue`) never contribute. This
guarantees `node[rand(1:length(node))]` always returns a valid subtree. Used by `getindex`,
`insert`, `insert_with_alignment`, and `targetsize` in `init`. `BinaryNode`/`UnaryNode`
cache it as the `size` field.

**`complexity(node)`** — effective expression size. Adds 1 for each non-zero `AddedValue`
on top of structural length. Used to bucket individuals into strata (`strata` is indexed by
complexity) and for any size/fit balance logic. `BinaryNode` and `UnaryNode` cache
`complexity` as an `Int` field; `Tree` also caches it.

- **`Constant` is special:** its `addition.value` IS the constant's value, not an offset on
  top of some other computation. Therefore `complexity(::Constant) = 1` always, even when
  `addition.value == 0`. A zero `Constant` is still a concrete numeric node.
- **`complexity(::AddedValue)`** encodes the conditional `value == 0 ? 0 : 1` — the only
  place that tests whether an additive constant contributes to complexity. Do not use
  `Base.length` on `AddedValue`.

**`pathlen_complexity(node)`** — complexity weighted by depth (each subtree's complexity
also counts toward its ancestors), i.e. `cmplx + pathlen_complexity(children...)`. Cached on
`BinaryNode`/`UnaryNode` and on `Tree`. Used as a tie-break / structural-balance signal in
reporting and selection; leaves reduce to `complexity`.
