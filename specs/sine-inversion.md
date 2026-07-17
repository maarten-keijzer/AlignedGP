# Constant fitting under `sin` nodes via interval stabbing — problem statement & current approach

## 0. Purpose of this document

This is a seed for a fresh discussion. It states the problem, the facts already established
(with worked counterexamples), the plan currently favoured, and the questions still open.
Code is Julia. Snippets are the actual working code, not pseudocode.

---

## 1. Domain

Genetic programming applied to symbolic regression. Alongside the usual forward evaluation
there is a **backward (inverse) pass** used to fit numeric constants inside candidate
expressions.

Setup:

- There are `n` fitness cases. Case `i` has input `X_i` (scalar) and a **target interval**
  `t_i :: IntervalType` — the band of output values that counts as a hit for that case.
- A **hit** is: the expression's output interval for case `i` intersects `t_i`. Fitness is
  the hit count, so it is coarse and integer-valued — ties are the common case, not the
  exception.
- To fit a constant `c` inside an expression, the target is pushed *backwards* through the
  expression tree to produce a **surrogate target** per case: the set of `c` values that
  would make case `i` a hit. Then find the `c` stabbed by the most surrogate targets.

**Stabbing** = given a set of intervals (one or more per fitness case), find the point
covered by the most of them, and return the maximal region(s) achieving that depth.
Max **depth**, not max clique — we want a *point*, so the Helly-property complications of
circular-arc graphs do not apply.

Constraint that keeps fitness meaningful:

> **Disjointness invariant.** The regions contributed by a *single* fitness case must be
> pairwise disjoint. Otherwise a case is counted twice at some point and `depth ≠ hits`.

This invariant is easy to violate subtly; §3.2 covers the case that shapes `sin_rev`'s design.

---

## 2. The specific problem

Node `sin(f(X))`, target `t`. We want to find `c` such that

```
sin(f(X) + c)
```

is optimal — i.e. `c` maximises the number of cases `i` where `sin(f(X_i) + c) ∈ t_i`.

Procedure: compute `t_s = sin_rev(t)` (arcs in u-space, where `u = f(X) + c`), then for
case `i` subtract the known forward value `u_i = f(X_i)` to get arcs in c-space, then stab.

Afterwards we want to recurse: feed `t_s - c` back through `f` to fit a deeper constant,
e.g. `c2` in `sin(f(X + c2) + c)`.

The whole difficulty is **periodicity**: `sin`'s true preimage is the principal arcs
`+ 2πk` for all integers `k`, and `sin_rev` returns only one representative per class.

---

## 3. Established findings

### 3.1 The affine test — which constants live on a circle

The path from the `sin` node down to the constant determines everything:

- Path is **affine** in the constant (`c ↦ αc + β`, i.e. only `+`/`-` and constant `*`):
  the constant is only defined **mod `2π/|α|`**. For `sin(f(X) + c)`, `α = 1`, period `2π`.
  `c` and `c + 2π` are *the same constant*. Folding mod 2π loses nothing.
- Path contains **anything nonlinear** (`sqrt`, `exp`, `log`, nested `sin`): the period is
  dead. Distinct `k` give genuinely distinct constants. Must choose one — inherently
  arbitrary, accepted as a limitation.

For `sin(sqrt(X + c1) + c2)`: `c2` is affine (circle, period 2π); `c1` is not (line).
The two constants need opposite treatment even though they sit under the same `sin`.

### 3.2 Why `sin_rev` merges at the crest and the dip

`sin_rev` (§6.2) returns **canonical arcs**: at most two, and merged wherever the target
reaches an extremum. This is a deliberate design property that the rest of the scheme
depends on, so it is worth stating why rather than leaving it as an implementation detail.

Consider `y = [-1.5, -0.5]` — dip covered, crest not. The naive treatment sends this to the
generic branch, giving two arcs:

```
rising  ≈ [-π/2, -π/6]
falling ≈ [ 7π/6, 3π/2]
```

These look disjoint on the line. Mod 2π, `[-π/2, -π/6]` becomes `[3π/2, 11π/6]` — so they
are really **one contiguous arc `[7π/6, 11π/6]` through the trough, touching at `3π/2`**.
The point `3π/2` would be in both → that case counted **twice** → depth ≠ hits, violating
the disjointness invariant of §1. `sin_rev` therefore merges the dip around `3π/2`,
exactly mirroring the long-standing crest merge around `π/2`, and returns the single arc.

Same reasoning at the top of the range: when `yl <= -1 && yh >= 1` every `x` is a preimage.
Deriving that arc from the branch formulas yields an inner-rounded width just under 2π,
which leaves a hairline gap once a consumer reduces mod 2π — so the condition is tested on
the *target* and a genuinely full arc is returned.

Resulting arc counts. Merging never costs arcs, and its translates never touch:

| target | arcs | note |
|---|---|---|
| generic (`-1 < cl`, `ch < 1`) | 2 | separated by both crest and dip; disjoint mod 2π |
| crest (`yh ≥ 1`) | 1 | merged around `π/2` |
| dip (`yl ≤ -1`, `yh < 1`) | 1 | merged around `3π/2`; extends past 2π when `ch > 0` |
| full (`yl ≤ -1 ≤ 1 ≤ yh`) | 1 | full circle, exact coverage |

**Consequence for consumers:** output now ranges over roughly `[-π/2, 5π/2]`, not
`[-π/2, 3π/2]`. Anything assuming the tighter window must reduce or translate (§4).

### 3.3 Folding mod 2π is exact for the direct constant — and needs no new stab routine

For `c` in `sin(f(X) + c)`: reduce every arc mod 2π, then emit each arc **twice** — once
where it is, once translated by `-2π`. Feed the result to the *existing linear* stab
routine. The duplicate is what handles the seam: a region wrapping past 0 appears as one
unbroken interval in the shifted copy, so the linear sweep finds it whole. No merge pass,
no second stabbing routine, no `nfull` bookkeeping in the sweep itself (full circles are
counted separately since both copies of a ≥2π arc would cover the same `c`).

### 3.4 Why a fixed window like `[-π, 3π]` degrades — the moving-window argument

Tempting alternative: have `sin_rev` return all arcs in a fixed window (say `[-π, 3π]`,
max 3 arcs), backpropagate those, and use one plain linear stab everywhere. One rev
function, one stab routine, no coordination.

This is **sound but lossy**, and the loss is not where you'd expect. The window is fixed in
*u-space*, but each case then subtracts its own `u_i = f(X_i)`. So in **c-space** case `i`
can only vote inside:

```
[-π - u_i,  3π - u_i]
```

**Every case gets its own window, in a different place.** Two cases whose `u_i` differ by
more than 4π ≈ 12.6 have disjoint voting windows and can *never* be counted together, no
matter how much they agree about the right constant. Depths at different `c` values are
therefore computed from *different subsets of the data*, which breaks the argmax.

Worked counterexample (13 cases, all targets `[0.99, 1.0]`, single arc `[1.43, 1.71]`/period):

- Bin A — 5 cases, `u ≈ 0`     → hit when `c ≡ 1.55` (mod 2π)
- Bin B — 4 cases, `u ≈ 47.33` → hit when `c ≡ 4.50`
- Bin C — 4 cases, `u ≈ 97.60` → hit when `c ≡ 4.50`

Truth: `c = 4.50` hits **8** cases (B and C agree). `c = 1.55` hits **5**.
Fixed window: A votes at ≈1.55; B at ≈-45.75; C at ≈-96.2 (all ≡ 4.50, never stacked).
Sweep reports depth 5 / 4 / 4 → picks `c = 1.55` → **5 true hits when 8 were available**.

So the answer to "can a constant have fewer true hits but more window hits?" is **yes**:
at `c = 1.55`, 5 window-hits / 5 true hits; at `c = 4.50`, 4 window-hits / 8 true hits.
It happens whenever cases that agree with each other have `u_i` values far apart.

What still holds for the fixed-window scheme:

1. **Reported depth is never an overcount** — whatever is picked has at least that many
   true hits. GP re-evaluates the full expression, so fitness is never inflated; the only
   loss is picking a worse constant than was available.
2. **Exact if `max(u) - min(u) ≤ 2π ≈ 6.28`** across the fitness cases.
3. Degradation is smooth, roughly: you choose based on ~1 bin out of `M = spread/2π` bins.
   Spread 6 → exact. Spread 50 → ~8 bins → choosing on ~1/8 of the data. Spread 500 → noise.

**This is measurable, and measuring it is a concrete open action:** log
`min/max` of `u_i = f(X_i)` at each `sin` node across fitness cases. `sqrt` squashes hard
(`X ∈ [0, 40]` → spread 6.3, fine); `exp` blows the window instantly.

### 3.5 The branch-explosion concern, and where it does *not* apply

Concern: materializing all `k` over a wide inner range × 1000 fitness cases → 50,000
surrogate targets.

This **does not** touch the direct-constant stabbing path, for two reasons:
(a) `u_i = f(X_i)` is a *scalar* per case, not an interval — there is no range to sweep `k`
across; (b) `c` is only defined mod 2π anyway, so `k`-translates are literally the same
constant written many times. ≤2n arcs, one sort.

It *does* apply when recursing under `f` to fit a deeper constant, where the inner value
becomes an interval (driven by the prior on that constant). Note the sensitivity:
`∂/∂c1 sqrt(X + c1) = 1/(2·sqrt(X+c1))` — **sqrt is contractive**, so the range stays
narrow (`X = 100`, `c1 ∈ [-10, 10]` → `u ∈ [9.49, 10.49]`, 0.16 periods → 1 branch).
The explosion needs an *expansive* path: `exp(X + c)` with `X = 5`, `c ∈ [-10, 10]` →
`exp([4,6]) = [54.6, 403]` → **55 periods**. That is the case where a `kmax` bail is the
honest answer, because across 55 periods `sin` genuinely does not localise `c` and any
finite branch set is arbitrary.

Also note: 50,000 arcs is not by itself fatal *if* they are **pooled into one stab** rather
than **branched into separate subtree explorations**. One `sort!` of 50k intervals is
milliseconds. Pooling is only legitimate while the disjointness invariant (§1) holds.

### 3.6 A correction on the record

An earlier line of argument in the previous discussion claimed a split (principal-branch)
rev should be retained for *domain contraction*, because merging "relocates" arcs by 2π.
That argument was **wrong** and is withdrawn; noted here so it does not get re-derived.

Both forms return one representative per class mod 2π and both under-approximate the true
preimage — neither window is privileged. For `sin(sqrt(x))` with `t = [-1.5, -0.5]` the
split form is in fact *worse*: its rising arc `[-π/2, -π/6]` is entirely negative, gets
killed by `sqrt`'s range, and loses `x ∈ [22.21, 33.18]` — e.g. `x = 30` is a real solution
(`sin(sqrt(30)) = -0.717 ∈ t`) that the split form drops and the merged form finds, since
the merged dip arc `[7π/6, 11π/6]` is entirely positive. The real distinction is not
split-vs-merged but **whether the consumer reduces mod 2π**.

---

## 4. Current plan

**One rev function** — `sin_rev` (§6.2), returning canonical arcs: crest → 1 arc merged
around π/2; dip → 1 arc merged around 3π/2; generic → 2 arcs; full → 1 full-circle arc.
Source of truth for all consumers. Output ranges over roughly `[-π/2, 5π/2]`; consumers
reduce or translate as needed. **This is implemented.**

**One stab routine** — the existing linear `max_overlap_region`, unchanged.

Three situations, differing only in preprocessing before the stab:

| situation | treatment |
|---|---|
| **1. Fit `c`, direct child of `sin`** (affine path) | Fold arcs mod 2π, duplicate each at `-2π`, feed to linear stab. **Exact.** |
| **2. Backpropagate `t_s` down through `f`** | Must pick one translate per case. Do **not** use a fixed window — **anchor to the data**: pick the translate nearest the known forward value `u_i = f(X_i) + c`. One arc per case, no explosion, least-arbitrary choice available. |
| **3. Fit `c2`, under `f`** (nonlinear path) | Plain linear stab on the anchored arcs. No folding — `f` killed the periodicity. Accept the loss. |

Coordination cost is **one scalar threaded down the backward pass**: set period `= 2π` at
the `sin` node; unchanged through `+ c` / `- X`; halved through `* 2`; set to `nothing`
through `sqrt`/`exp`/`log`/nested `sin`. At the leaf: `nothing` → anchor+translate,
otherwise → fold.

**Known coupling:** `c` is found first, then `t_s - c` is backpropagated; the anchoring in
(2) uses `c`, so the passes are coupled — arcs sent down depend on the constant just
picked. Probably desirable (local refinement around the current expression), but it means
re-running the `c` pass after `c2` moves would shift the anchors. Iterate or don't; one
pass is a reasonable local search step.

**Tie-breaking:** prefer the **widest** deepest region (fitness is coarse → ties are the
norm; a wide plateau gives a `c` robust to later refinement). Filter zero-width components
— a constant good at a single exact real is worthless. Canonicalise the chosen `c` into
`(-π, π]` for a smaller printed literal and better conditioning.

**Sibling rev functions:** if there is a dispatch table mapping node type → rev function,
`cos_rev` / `tan_rev` must agree with `sin_rev` about which convention they are in, or the
walker cannot tell which entries return what.

---

## 5. Open questions

1. Measure the spread of `u_i = f(X_i)` at real `sin` nodes across fitness cases. If it is
   consistently ≤ ~6, the fold in situation 1 buys little and the simplest scheme wins. If
   routinely ≥ 12, the fold is free accuracy.
2. Is anchoring (situation 2) actually better than a fixed window in practice, or does the
   coupling to `c` cause thrash across generations?
3. `kmax` bail policy for expansive paths (`exp`): bail to `hull(arcs) ∩ U`? Or skip
   constant fitting at that node entirely?
4. Does the disjointness invariant survive every rev function on the path? Worth a debug
   assertion (sort per-case regions, check no two touch) rather than an assumption. Note
   `v ↦ v²` folds at `v = 0` and needs a `v ≥ 0` clip to stay monotone when inverting
   `sqrt`.
5. Nested `sin` — currently intractable, presumably just kill the period and accept.

---

## 6. Code

Julia. Interval type is `IntervalType`, constructed via `intervaltype(lo, hi)` — **note the
constructor requires `lo <= hi` and throws otherwise**. Available operations:
`bounds(y) -> (lo, hi)`, `inf(x)`, `sup(x)`, `mid(x)`, `issubset_interval(a, b)`,
`intersect_interval(a, b)`, `isempty_interval(a)`, `in_interval(x, a)`. Rounding is
rigorous; `intervaltype(π)` is a genuine π enclosure.

Note the codebase currently mixes `bounds(y)` and direct field access `iv.lo` / `iv.hi`.
On recent IntervalArithmetic.jl the fields are not public API — `inf`/`sup`/`bounds` are.
Worth unifying.

`sup`/`inf` are used deliberately to **inner-round** the arcs (`sup` on the lower end, `inf`
on the upper end), so the returned arcs are guaranteed-sound subsets rather than
outward-rounded supersets.

### 6.1 The stabbing routine (linear, current, works)

```julia
function max_overlap_region(intervals::Vector{IntervalType})

    isempty(intervals) && return (; region=IntervalType[], depth=0)

    # Encode: (coordinate, 0) = enter, (coordinate, 1) = leave
    events = Vector{Tuple{Float64,Int}}(undef, 2 * length(intervals))
    for (i, iv) in enumerate(intervals)
        events[2i-1] = (iv.lo, 0)
        events[2i]   = (iv.hi, 1)
    end
    sort!(events)

    depth     = 0
    max_depth = 0
    seg_start = 0.0
    components = IntervalType[]

    for (coord, kind) in events
        if kind == 0        # enter
            depth += 1
            if depth > max_depth
                max_depth = depth
                empty!(components)
                seg_start = coord
            elseif depth == max_depth
                seg_start = coord
            end
        else                # leave
            if depth == max_depth
                push!(components, intervaltype(seg_start, coord))
            end
            depth -= 1
        end
    end

    return (; region=components, depth=max_depth)
end
```

The core loop is correct and is intended to survive unchanged. Note `sort!` on
`Tuple{Float64,Int}` sorts enters (0) before leaves (1) at equal coordinates, which is the
right convention for closed intervals.

### 6.2 `sin_rev` — canonical arcs (current, implemented)

Returns canonical arcs per §3.2: merged at crest and dip, full circle caught up front.
Output ranges over roughly `[-π/2, 5π/2]`. Single source of truth for all three situations
in §4.

```julia
const TWOPI       = 2 * intervaltype(π)
const FULL_CIRCLE = intervaltype(0.0, sup(TWOPI))

function sin_rev(y::IntervalType) :: Union{Nothing, IntervalType, Tuple{IntervalType, IntervalType}}
    yl, yh = bounds(y)

    cl = max(yl, -1.0)          # clamp target to sin's range [-1, 1]
    ch = min(yh,  1.0)
    cl > ch && return nothing   # target unreachable by sin

    peak   = yh >= 1.0          # target reaches the crest
    trough = yl <= -1.0         # target reaches the dip

    # Both: y ⊇ [-1, 1], so every x is a preimage. Return a genuinely full arc —
    # deriving it from the branches below gives an inner-rounded width just under
    # 2π, which leaves a hairline gap once the consumer reduces mod 2π.
    peak && trough && return FULL_CIRCLE

    PI = intervaltype(π)

    if peak
        # Crest covered: rising+falling merge into one arc spanning x = π/2,
        # so no spurious hole at the peak.
        A = asin(intervaltype(cl))
        lo = sup(A)                  # asin(cl) inward (up)
        hi = inf(PI - A)             # π−asin(cl) inward (down)
        return lo <= hi ? intervaltype(lo, hi) : nothing

    elseif trough
        # Dip covered: mirror of the peak case. Rising+falling merge into one arc
        # spanning x = 3π/2. Emitted as the two separate arcs they would otherwise
        # be, they touch at 3π/2 once reduced mod 2π, double-counting the case there.
        # May extend past 2π (whenever ch > 0); that is correct and expected.
        B = asin(intervaltype(ch))
        lo = sup(PI - B)             # π−asin(ch) inward (up)
        hi = inf(TWOPI + B)          # 2π+asin(ch) inward (down)
        return lo <= hi ? intervaltype(lo, hi) : nothing
    end

    # Neither crest nor dip: -1 < cl and ch < 1, so the two arcs are separated by
    # both the peak and the trough region and stay disjoint mod 2π.
    Acl, Ach = asin(intervaltype(cl)), asin(intervaltype(ch))

    rlo, rhi = sup(Acl), inf(Ach)             # rising arc  ⊆ [−π/2, π/2]
    flo, fhi = sup(PI - Ach), inf(PI - Acl)   # falling arc ⊆ [π/2, 3π/2]
    rising, falling = rlo <= rhi, flo <= fhi

    if rising && falling
        return (intervaltype(rlo, rhi), intervaltype(flo, fhi))
    elseif rising
        return intervaltype(rlo, rhi)
    elseif falling
        return intervaltype(flo, fhi)
    else
        return nothing
    end
end
```

Notes for anyone modifying this:

- `Acl` / `Ach` are computed inside their branches rather than hoisted — the crest branch
  needs only `asin(cl)`, the dip branch only `asin(ch)`; hoisting either wastes an `asin`
  on the other path.
- The branch guards keep `asin` off its domain edges. Inside the dip branch `peak` is
  false, so `yh < 1` and `ch = yh`, putting `asin(ch)` strictly inside `(-π/2, π/2)`.
  Symmetrically for the crest branch.
- The two independent validity checks in the generic branch (`rlo <= rhi` and
  `flo <= fhi`) are worth keeping: they are separate inner-rounded computations, so on a
  very narrow target one can collapse while the other survives.

Sanity checks: `y = [-1.5, -0.5]` → dip branch → `B = asin(-0.5) = -π/6` →
`[π + π/6, 2π - π/6] = [7π/6, 11π/6]`, width `2π/3`, `3π/2` interior and covered once.
`y = [-1.5, 0.5]` → `[5π/6, 13π/6]`, which wraps past `2π` as intended. (The dip arc
exceeds `2π` only when `ch > 0`.)

### 6.3 The fold, for situation 1 (proposed, §3.3)

```julia
# arcs :: Vector{Tuple{Float64,Float64}} of (start, width) in c-space, pre-reduction
function fold(arcs, C = 2π)
    out, nfull = IntervalType[], 0
    for (s, w) in arcs
        if w >= C
            nfull += 1               # every c works for this case
        else
            s = mod(s, C)
            push!(out, intervaltype(s, s + w))
            push!(out, intervaltype(s - C, s - C + w))
        end
    end
    out, nfull
end

folded, nfull = fold(arcs)
res = max_overlap_region(folded)     # existing routine, untouched
depth = res.depth + nfull
```

Width must be computed **before** reduction and carried; do not recover it from reduced
endpoints, and do not `mod` the two endpoints separately (that splits the arc).

### 6.4 Test for the dip merge (current behaviour)

Asserts the canonical single arc and single coverage at `3π/2`. (Supersedes an older test
that asserted two arcs `[-π/2, -π/6]` and `[7π/6, 3π/2]` with `sup(a) < inf(b)` — a linear
picture of a circular object; if that test is still around anywhere, it is stale.)

```julia
@testset "dip covered merges through 3π/2 (one arc)" begin
    y = intervaltype(-1.5, -0.5)          # dip included, crest not
    x = sin_rev(y)
    @test x isa IntervalType               # merged, not split
    lo, hi = bounds(x)
    @test lo ≈ 7π/6  atol=1e-9             # π − asin(−0.5)
    @test hi ≈ 11π/6 atol=1e-9             # 2π + asin(−0.5)
    @test lo <= hi
    @test in_interval(3π/2, x)             # dip is interior, covered once
    # soundness: the arc maps back inside the target
    @test issubset_interval(sin(x), y)
end;
```

Caution on the soundness assertion: the arcs are inner-rounded but `sin` rounds outward,
and at `11π/6` the derivative is `cos(11π/6) = √3/2 ≈ 0.87`, so the two effects are the
same order — whether `sup(sin(x)) <= -0.5` holds exactly comes down to the last ulp.
Asserting against the unclamped target `y` (as above) has room on both sides; asserting
against the clamped `[-1.0, -0.5]` needs explicit slack. The `-1.0` end is safe regardless,
because `3π/2` is interior and a decent interval `sin` returns exactly `-1` there.

Worth adding alongside: a wrapping dip (`y = intervaltype(-1.5, 0.5)` → `[5π/6, 13π/6]`,
`hi > 2π`), and the full circle (`y = intervaltype(-2, 2)` → `FULL_CIRCLE`, width ≥ 2π,
no hairline gap after `mod`).