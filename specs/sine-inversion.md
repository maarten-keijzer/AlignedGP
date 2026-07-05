# Sine Inversion for Offset Stabbing — Implementer's Spec

## TL;DR

To stab the offset $c$ on the argument of $\sin(s + c)$, do the max-overlap
**on the circle $\mathbb{R}/2\pi\mathbb{Z}$, not on the real line.** One period of
search is provably sufficient. The only correctness requirement is that each
point's surrogate is *aligned by its current argument mod $2\pi$* before
overlapping. Clipping each surrogate to one fixed period of the line is the bug
that breaks the invariant.

---

## What we're inverting

We are stabbing an additive offset $c$ at the argument of a sin node:

$$\sin(s_i + c) \in \tau_i,$$

where for point $i$:
- $s_i = v_i$ is the **current value of the argument subtree** (a known scalar),
- $\tau_i$ is the **target for the sin's output**: the node's target interval(s)
  with any already-fixed outer constant removed. $\tau_i \subseteq \mathbb{R}$,
  but only its intersection with $[-1,1]$ can ever be satisfied.

Goal: find $c$ hitting the maximum number of points.

---

## Step 1 — One-period preimage $P_i$

Compute the set of arguments in $[0, 2\pi)$ whose sine lands in $\tau_i$:

$$P_i = \{\, u \in [0,2\pi) : \sin u \in \tau_i \,\}.$$

For a target interval $\tau_i = [\ell, h]$, clamp to the sine's range,
$[\ell', h'] = [\max(\ell,-1), \min(h,1)]$:

- If $\ell' > h'$ (target misses $[-1,1]$ entirely): $P_i = \varnothing$.
  **This point is unwinnable by this subtree — drop it from the overlap.**
  (Correct behavior, not an error.)
- Otherwise $P_i$ is the set with $\sin u \in [\ell', h']$. Using
  $\theta_\ell = \arcsin(\ell')$, $\theta_h = \arcsin(h') \in [-\tfrac\pi2,\tfrac\pi2]$,
  the preimage over one period is **up to two arcs**:

  - rising branch: $u \in [\theta_\ell,\ \theta_h]$ (near $\tfrac\pi2$),
  - falling branch: $u \in [\pi - \theta_h,\ \pi - \theta_\ell]$ (near $\tfrac{3\pi}{2}$),

  then reduced into $[0,2\pi)$. Edge cases: $h' = 1$ merges the two arcs at the
  peak $u=\tfrac\pi2$; $\ell' = -1$ merges them at the trough $u=\tfrac{3\pi}2$;
  a degenerate target ($\ell'=h'$) gives arcs of zero width (two points) — widen
  by the target tolerance so the stabber has interior to sample.

If $\tau_i$ is itself a union of intervals/rays, compute $P_i$ per component and
union the arcs.

---

## Step 2 — Align onto the circle (the critical step)

The offset that hits point $i$ is any $c$ with $s_i + c \in P_i + 2\pi\mathbb{Z}$,
i.e.

$$c \in A_i := (P_i - s_i) \bmod 2\pi \quad\subseteq\ [0, 2\pi).$$

**Subtract $s_i$ and reduce mod $2\pi$.** This is what places points from
different periods onto a common fundamental domain. Two points whose current
arguments differ by a multiple of $2\pi$ but want the same offset will now have
*overlapping* arcs; on the un-reduced line they would look disjoint.

Each $A_i$ is one or two arcs on the circle; an arc may wrap across $0 \equiv
2\pi$.

---

## Step 3 — Circular max-overlap stab

Find $c^\ast \in [0,2\pi)$ covered by the most arcs $A_i$:

1. Turn each arc into `(start, +1)` and `(end, -1)` events.
2. For an arc that wraps past $2\pi$, either split it into
   $[\text{start}, 2\pi) \cup [0, \text{end}]$, **or** duplicate all arcs shifted
   by $+2\pi$ and sweep over $[0, 4\pi)$ — standard circular-stabbing unrolling.
3. Sweep, tracking a running cover count; take the max.

Complexity $O(n \log n)$. This is your existing 1-D stabber plus wraparound.

---

## Step 4 — Pick $c^\ast$

Sample $c^\ast$ from the **interior** of a max-overlap arc (never an endpoint) to
avoid open/closed-boundary and floating-point issues. Random selection within the
max-overlap region is fine — every point of it attains the maximum. Then apply
$c^\ast$ to the argument subtree as the additive offset.

---

## Why one period is enough (and multiple periods are not a problem)

The count function $N(c) = |\{i : c \in A_i\}|$ is $2\pi$-periodic, because every
$A_i$ is. Hence

$$\max_{c \in \mathbb{R}} N(c) = \max_{c \in [0,2\pi)} N(c).$$

Searching one period of the **offset** loses nothing — the constant can shift the
argument into whichever period is best. Arguments spanning many periods are
handled automatically by the mod-$2\pi$ alignment in Step 2; they are **not** a
reason to invert more than one period.

---

## The bug to avoid

Do **not** compute the literal real-line interval $P_i - s_i$ and stab without
wrapping. That clips each point to one fixed period of the line, so points in
different periods that share an optimal offset appear disjoint. The surrogate
becomes a strict *subset* of the true preimage → false negatives → the stabber
undervalues the good offset (including $c=0$) and can pick a worse one, dropping
actual hits below the starting count and **breaking the monotonicity invariant.**

Minimal failing case: target $\{1\}$, current arguments
$\{\tfrac\pi2,\ \tfrac\pi2 + 2\pi,\ 0\}$. Circular stab keeps both currently-hit
points (hits stays 2); line-clipped stab sees three disjoint singletons, and a
tie-break to $c^\ast=\tfrac\pi2$ drops to 1 hit.

---

## Self-check (drop-in oracle)

After the stab, assert that the count the stabber reported ($M$) equals the hits
you get by **direct evaluation** of $\sin(s_i + c^\ast)$ against $\tau_i$:

```
assert reported_M == sum(sin(s_i + c_star) in tau_i for all i)
```

If these disagree, the surrogate is not the true preimage. For sin, the
disagreeing points will be ones whose current argument $s_i$ sits a nonzero
multiple of $2\pi$ from the base period — i.e. Step 2 alignment is missing or
wrong.

---

## Checklist

- [ ] $P_i$ computed over exactly one period $[0,2\pi)$; up to two arcs.
- [ ] Target-out-of-$[-1,1]$ → empty preimage → point dropped, not forced.
- [ ] Peak/trough arc-merging handled ($h'=1$, $\ell'=-1$).
- [ ] Degenerate targets widened by tolerance.
- [ ] Surrogate aligned as $(P_i - s_i)\bmod 2\pi$ **before** overlap.
- [ ] Stab is circular (wraparound), not linear.
- [ ] $c^\ast$ sampled from arc interior.
- [ ] Post-stab oracle assertion in place.
