# AlignedGP — development notes

## Invalid interval sentinel: `CI(typemax(T), typemax(T))`

When an inverse function has no valid solution (e.g. `sqrt(x) < 0` is impossible),
it must return an **invalid sentinel** rather than `nothing` or a misleading valid
interval. The sentinel is `CInterval{T}(typemax(T), typemax(T))` — for `Float64`
this is `CI(Inf, Inf)`.

**Why not `nothing`?** Returning `nothing` forces every caller to check for it,
breaking the uniform interval arithmetic pipeline.

**Why not a clamped valid interval?** Returning e.g. `CI(0, 0)` implies `x = 0` is
a solution when it isn't. This silently inflates hit counts.

**Why `CI(Inf, Inf)`?** It propagates correctly through interval arithmetic (stays
at infinity under addition/subtraction) and is unambiguously detectable via
`isinf(first(t)) && first(t) == last(t)`. Note that `CI(-Inf, Inf)` (the whole
real line) is a **valid** interval and must not be treated as a sentinel.

## `max_overlap_region` filters invalid sentinels automatically

`max_overlap_region` detects and removes invalid sentinels before the sweep using
`_is_invalid(t) = isinf(first(t)) && first(t) == last(t)`. If **all** input
intervals are invalid, it returns `(; region=CIntervalSet([CI(Inf,Inf)]), depth=0)`.
Callers use `depth == 0` to detect that no valid constant exists.

## DomainErrors are caught at the Tree level

Node-level `evaluate` (e.g. `sqrt.(x)`) may throw `DomainError` for invalid inputs
such as negative values. These are **not** caught at the node level — they propagate
up to `evaluate(tree::Tree, x)`, which catches `DomainError` and returns
`fill(Inf, length(x[1]))`. This produces `CI(Inf, Inf)` sentinels in the shifted
interval vector, which `max_overlap_region` then filters, yielding `depth = 0`.

## Intervals.jl cannot represent empty bounded-closed intervals

`Interval{T,Closed,Closed}` validates and sorts its constructor arguments, so a
traditional "empty by construction" trick (`CI(high, low)`) does not work. Use the
`CI(typemax(T), typemax(T))` sentinel for unreachable cases instead.
