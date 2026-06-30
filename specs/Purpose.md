# Additive Alignment for Interval-Target Symbolic Regression

A reference for the **additive-alignment** mechanism. Scope is deliberately narrow:
the hits objective, the closed-form constant fit, the max-overlap region computation,
the constant-selection policy, and the structural consequence (a single arbitrary-constant
terminal). Search strategy, target/surrogate construction, and population design are out of
scope here.

---

## 1. Problem setup

Targets are **intervals**, not point values. For each data point `i` we are given a band

```
band_i = [l_i, u_i],   W_i = u_i - l_i   (the band width)
```

A prediction `f(x_i)` either **hits** band `i` (lands inside `[l_i, u_i]`, edges inclusive)
or **misses** it. Bands are *hard*: every value inside is equally acceptable, with a cliff at
the edge. This is a tolerance/spec, not a blurred point estimate.

### Hits loss (objective)

The model fitness is a **count of hits** (equivalently, a weighted miss-count to minimize):

```
L(f) = sum_i  rho_i * 1[ f(x_i) NOT in band_i ]
```

The per-point weight comes from a uniform-inlier / uniform-outlier mixture likelihood:

```
rho_i = log(1 + kappa / W_i)      kappa = a single width-scale parameter (units of y)
```

Key behaviors:
- `W_i << kappa` → `rho_i ≈ log(kappa/W_i)`: tight bands are strong evidence (high weight).
- `W_i >> kappa` → `rho_i ≈ kappa/W_i → 0`: inflated bands are nearly worthless; degenerate
  cases self-neutralize.

**Homoscedastic case (the default).** When all widths are equal (`W_i ≡ W`), `rho_i` is
constant, `kappa` drops out of the argmin, and the loss is just `rho * (#misses)`. In that
regime **treat the weight as 1 and score by raw hit count**. Use the weighted form only if
the target bands have differing widths; otherwise it is parameter-free.

### Why hits, not squared error

The hit relation is preserved exactly under monotone bijective transforms of the prediction
(the kind introduced when reasoning about a node inside a larger expression):

```
f(x_i) in [l_i, u_i]   <=>   g(f(x_i)) in g([l_i, u_i])      for g strictly monotone per point
```

This is a *biconditional on the hit set*, so the hit count is invariant. Squared error is not:
its argmin moves under the transform's Jacobian. This invariance is the property the whole
approach relies on — hit counts computed in a transformed coordinate equal real hit counts in
the original coordinate.

---

## 2. Additive alignment

Each subtree `S` is evaluated at the data points to produce **point values** `S(x_i)`
(subtrees are point-based, not interval-based). When placing `S` against a target, we fit a
single **additive constant** `c` so that `S + c` maximizes hits:

```
maximize over c:   #{ i :  S(x_i) + c  in  [l_i, u_i] }
```

This decouples **shape** (the subtree's x-dependence) from **offset** (the DC level):
the subtree only has to match the *shape* of the target; the constant absorbs placement.

### Reduction to 1-D max overlap

`S(x_i) + c in [l_i, u_i]` is equivalent to:

```
c in [ l_i - S(x_i),  u_i - S(x_i) ]
```

So fitting `c` is exactly: **find the value `c` contained in the maximum number of shifted
intervals**

```
shifted_i = [ l_i - S(x_i),  u_i - S(x_i) ]
```

This is a standard 1-D max-overlap (interval stabbing) problem, solved in `O(n log n)`.

> **Cost note:** only the *additive* constant is fit in closed form. Fitting a multiplicative
> scale (`a*S + c`) is **not** done — it is `O(n^2)` and is explicitly excluded for
> scalability. Leading coefficients must be expressed structurally, not fit. See §6.

---

## 3. Max-overlap algorithm (full region, not just a point)

A sweep computes not just the best `c` but the **entire set of optimal `c` values** — every
constant achieving the maximum hit depth. This full region is needed for the selection policy
in §4 and is produced at no extra cost.

### Sweep

1. Build `2n` events from the shifted intervals: a `+1` (enter) at each left endpoint, a `-1`
   (leave) at each right endpoint.
2. Sort events by coordinate. **Tie-ordering at equal coordinates matters** (see edge cases).
3. Sweep left→right tracking running depth; record the maximum depth `k`.
4. The optimal set is the union of all coordinate ranges where depth == `k`.

### Output: a union of closed components

The max-depth region is **not necessarily a single interval**. Several disjoint
constant-ranges can each achieve the same max depth `k`:

```
optimal_c_set = [a_1, b_1] ∪ [a_2, b_2] ∪ ...
```

Emit **all** maximal components, each as a closed interval `[a_j, b_j]`, together with the
depth `k`. Do **not** stop at the first witness, and do **not** collapse to the convex hull —
the gap between two components can hit *fewer* targets than the components themselves.

### Edge cases (must handle for correctness)

- **Closed bands / boundary hits.** Bands are closed, so `c` exactly equal to a shifted
  endpoint still counts as inside. Order events so that at equal coordinates **enter (+1)
  events are processed before leave (-1) events**; otherwise touching bands give off-by-one
  depth. This is not hypothetical — many equal endpoints arise when several points share
  values.
- **Disjoint components.** Collect every maximal-depth run, not just the first.
- **Empty / degenerate intervals.** If `l_i == u_i` the shifted interval is a single point;
  it still participates (enter and leave at the same coordinate, enter first).

### Interval width as a free diagnostic

The width `b_j - a_j` of the chosen component is an **identifiability / robustness** readout,
available for free:
- **Wide** region → the constant is weakly constrained (many values equally good).
- **Narrow** region → the constant is sharply pinned.

This can be kept as a tiebreak between equal-hit placements (prefer the one that pins its
constant tighter) but is not required.

---

## 4. Constant selection policy

Given the optimal set (union of closed components) and depth `k`, pick the concrete `c`:

1. **If `0` lies in any component → choose `c = 0`.**
   The additive constant vanishes (no offset node needed; size saving).
2. **Otherwise → choose the midpoint of the widest component.**
   The midpoint is the most robust witness: it stays optimal under the largest symmetric
   perturbation, which matters if this node is later embedded in a larger expression.

**Do not** take the midpoint across disjoint components (the convex hull): that point can fall
in a hit-losing gap. Midpoint is correct *within* a component, prefer-zero / widest-then-mid
*across* components.

These are the correct **uninformed defaults**. The choice of `c` within its allowable interval
genuinely couples to the rest of the tree (a constant below an open slot co-determines what
that slot must do), so a context-aware resolution is possible later — see §7.

---

## 5. Point commitment — what alignment must NOT become

The output of an aligned subtree is always a **single committed point** `S(x_i) + c` per data
point. It is tempting to instead emit the optimal `c`-interval as an **interval-valued output**
and score a hit as "output interval overlaps band." **Do not do this.**

Overlap-scoring an interval output lets *each point pick a different witness* from the interval,
i.e. it silently treats `c` as a per-point variable rather than one committed constant. That
inflates the hit count above what any single constant achieves, and the inflated hits are not
realizable by any concrete expression. It breaks the core property that a score equals the real
loss of a committed function.

If deferring the constant is desired, carry the `[a, b]` interval as **deferred state** and
**commit to a single point at scoring time** (§7). Never score by interval overlap. (If an
interval semantics is ever genuinely wanted, the *sound* predicate is **containment** — output
range ⊆ band — not overlap; but the default and intended path is point commitment.)

---

## 6. The arbitrary-constant terminal

Because every node carries an additive alignment constant, a **separate numeric constant
terminal / ephemeral-random-constant subsystem is unnecessary**. Replace it with a single
terminal:

```
0   (size 0)   — "the optimal constant for this location"
```

Its value is determined by additive alignment at scoring time, not stored in the genome and
not mutated. This removes ERC generation, constant mutation, and the special-casing of numeric
leaves in crossover.

### Consequences

- **Additive-only boundary.** A bare `0` terminal gives a free *additive* constant. It does
  **not** give a free *multiplicative* coefficient (e.g. the `a` in `a*x + b`). Non-unit
  coefficients deep in the tree must be built structurally. This is consistent with the
  additive-only decision in §2 — be aware of it for targets that need genuine scale factors.
- **Constant folding is now a needed canonicalization.** The `0` terminal makes "is this
  subtree secretly just a constant?" come up often (`(0 + 0)`, `(0 * x)`, etc.). Run a
  **simplify-on-evaluate pass**: any subtree that evaluates to a constant collapses to the
  size-0 `0` terminal. Otherwise redundant constant-expressions accumulate size without
  improving hits.
- **`0` is the limiting near-constant subtree.** A subtree whose whole contribution is the
  alignment constant is legitimate and desirable as a leaf, but selection that rewards raw
  aligned hits can drift toward near-constant shapes that free-ride on alignment. If this
  becomes an issue, select on **lift over the bare-constant baseline** (score of `S + c*`
  minus score of the best plain constant) rather than raw aligned hits. Optional.

### Size accounting

A committed alignment constant is a real node (`+ c`), so a placed subtree's size is
`size(S) + (constant node)`. Convention: a `c = 0` constant costs nothing. If non-zero
constants are counted uniformly, that merely scales all sizes by a constant and is irrelevant
for comparisons — decide once and be consistent.

---

## 7. Deferred constant resolution (future, optional)

The selection policy in §4 is uninformed (prefer 0, else midpoint). A real optimization exists:
the best point in `[a, b]` is the one that best serves the *rest* of the tree, because a
committed constant below an open slot co-defines that slot's requirements.

The cheap version is **not** joint optimization (combinatorial), but a **resolution order**:

1. Keep alignment constants as their `[a, b]` intervals (the deferred state from §3).
2. Commit them outside-in, or at final scoring, picking each point against the constraints its
   context has by then accumulated.

Each individual commitment stays a closed-form point-pick (`O(n log n)`, no joint solve); only
the *timing* of collapse changes. The same deferred-interval machinery from §3 supports this,
so it is cheap to add onto a working aligned implementation later. Not needed for a first
version.

---

## Summary of invariants

- Subtrees are **point-valued**; targets are **interval-valued**.
- Fit **one additive constant** per placement, in closed form, via 1-D max overlap on shifted
  bands. `O(n log n)`. No multiplicative scale.
- Emit the **full max-depth region** as a union of closed components, plus depth `k`.
- Select `c`: **0 if available, else midpoint of widest component**. Never midpoint across
  components.
- Always **commit to a point** at scoring; never score interval outputs by overlap.
- A single **`0` terminal (size 0)** replaces all explicit constants; pair it with a
  **constant-folding** canonicalization pass.
- Score by **hit count** (uniform weight) in the homoscedastic case; use `rho_i = log(1 +
  kappa/W_i)` only when band widths differ.