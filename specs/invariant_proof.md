# The Monotone-Chain Invariant for Sequential Offset Stabs

**Claim.** In a symbolic-regression system where subtree constants are found by
sequential (Gauss–Seidel) offset stabbing against exact per-point surrogate
targets, the number of hit intervals is non-decreasing at every stab. Division
does **not** weaken this invariant; the proof never references the shape of the
surrogate set, so unions of rays are handled with no modification.

---

## 1. Setup and notation

Fix data points $\mathbf{x}_1,\dots,\mathbf{x}_n$ and target sets
$t_1,\dots,t_n \subseteq \mathbb{R}$. The targets may be intervals, unions of
intervals, or unions of rays — this is **irrelevant** to everything below.

Adjoin a symbol $\bot$ for "undefined" (e.g. division by zero) and decree

$$\bot \notin t_i \quad \text{for all } i.$$

For an expression $E$, write $E(\mathbf{x}_i) \in \mathbb{R}\cup\{\bot\}$ for its
evaluation on point $i$, and define the hit count

$$\text{hits}(E) := \bigl|\{\, i : E(\mathbf{x}_i) \in t_i \,\}\bigr|.$$

**Context function.** Let $T$ be the current tree and $p$ a node position in it.
For each point $i$, define

$$g_i^p : \mathbb{R} \to \mathbb{R}\cup\{\bot\}$$

as the map sending *the value produced at node $p$* to *the value produced at the
root*, with every other input and subtree frozen at its value for point $i$.
Concretely: evaluate $T$ at $\mathbf{x}_i$, but replace the subtree at $p$ by a
free variable $u$. The result is a fixed partial function of $u$, with
$g_i^p(u)=\bot$ wherever downstream evaluation is undefined.

Let $v_i^p$ be the actual value of the subtree at $p$ on point $i$, so that

$$g_i^p(v_i^p) = T(\mathbf{x}_i). \tag{1}$$

**Exact preimage surrogate.** Define

$$S_i^p := (g_i^p)^{-1}(t_i) = \{\, u \in \mathbb{R} : g_i^p(u) \in t_i \,\}. \tag{2}$$

This is a **definition**, not an approximation: $S_i^p$ is exactly the set of
node-$p$ outputs that make point $i$ a hit. Because $\bot \notin t_i$, every $u$
producing an undefined result is automatically excluded from $S_i^p$.

**The stab.** We stab an additive offset $c$ at node $p$: replace the subtree
$s$ at $p$ by $s + c$. The node's output on point $i$ becomes $v_i^p + c$, so
point $i$ is hit iff $v_i^p + c \in S_i^p$, i.e. iff $c \in S_i^p - v_i^p$. The
max-overlap routine returns

$$c^\ast \in \arg\max_{c\in\mathbb{R}}
\bigl|\{\, i : c \in S_i^p - v_i^p \,\}\bigr|,
\qquad M := \text{that maximum count.} \tag{3}$$

---

## 2. Assumptions

- **(A1) Exactness.** The surrogate the implementation feeds to the stabber
  equals $S_i^p$ in (2) — the true preimage of $t_i$ under $g_i^p$.
- **(A2) Correct stabbing.** The routine returns a point $c^\ast$ attaining the
  true maximum overlap $M$ of the sets $\{S_i^p - v_i^p\}_i$, and reports that
  count. Random selection within the argmax region is fine: every point of the
  argmax region attains $M$ by definition.

These are the **only** assumptions. (A1) is the one division can violate in an
implementation. Notice that the proof below does not touch it — which is exactly
why any observed violation is a bug in the surrogate computation, not a gap in
the theorem.

---

## 3. Lemma 1 — Stab count equals resulting hits

**Statement.** Let $T'$ be $T$ with node $p$ replaced by $s + c^\ast$. Then
$\text{hits}(T') = M$.

**Proof.** In $T'$, node $p$ outputs $v_i^p + c^\ast$ on point $i$, so the root
outputs $g_i^p(v_i^p + c^\ast)$. Hence

$$
\text{hits}(T')
= \bigl|\{\, i : g_i^p(v_i^p + c^\ast)\in t_i \,\}\bigr|
\overset{(2)}{=} \bigl|\{\, i : v_i^p + c^\ast \in S_i^p \,\}\bigr|
= \bigl|\{\, i : c^\ast \in S_i^p - v_i^p \,\}\bigr|
\overset{(3)}{=} M.
$$

The middle equality is where exactness (A1) is used, and it is a **pure
membership equivalence** — it holds verbatim whether $S_i^p$ is an interval, a
union of two rays, a punctured line, or the empty set. $\blacksquare$

---

## 4. Lemma 2 — Zero is feasible

**Statement.** $M \ge \text{hits}(T)$.

**Proof.** $c = 0$ is a candidate in (3), and

$$
\bigl|\{\, i : 0 \in S_i^p - v_i^p \,\}\bigr|
= \bigl|\{\, i : v_i^p \in S_i^p \,\}\bigr|
\overset{(2)}{=} \bigl|\{\, i : g_i^p(v_i^p)\in t_i \,\}\bigr|
\overset{(1)}{=} \text{hits}(T).
$$

Since $M$ is the maximum over all $c$, we get $M \ge \text{hits}(T)$.
$\blacksquare$

---

## 5. Theorem — Monotone chain

**Statement.** Let $T_0$ be any tree and let
$T_0 \to T_1 \to \cdots \to T_m$ be a sequence in which each $T_k$ is obtained
from $T_{k-1}$ by stabbing an additive offset at one node position $p_k$, using
the exact preimage surrogate computed **against $T_{k-1}$** (A1) and a correct
stabber (A2). Then

$$\text{hits}(T_m) \ge \text{hits}(T_{m-1}) \ge \cdots
\ge \text{hits}(T_1) \ge \text{hits}(T_0).$$

**Proof.** Fix a step $k$. Apply Lemma 1 with $T=T_{k-1}$, $p=p_k$, $T'=T_k$:
this gives $\text{hits}(T_k)=M_k$. Apply Lemma 2 with the same data:
$M_k \ge \text{hits}(T_{k-1})$. Hence $\text{hits}(T_k)\ge\text{hits}(T_{k-1})$.
Chain over $k=1,\dots,m$. $\blacksquare$

The hypothesis "surrogate computed against $T_{k-1}$" is the Gauss–Seidel
condition. It is essential: it is what makes $c = 0$ reproduce $T_{k-1}$ exactly
in Lemma 2. Computing a step's surrogate against a stale version of the tree
(Jacobi) removes the feasible zero and voids the chain.

---

## 6. Corollary — the three-constant operator

Take $T_0 = X \circ Y$ (with $\circ = \times$, or $/$, or any operator whose
per-point inverse you compute exactly), and the sequence

$$
T_0=X\circ Y
\;\to\; T_1=(X{+}c_X)\circ Y
\;\to\; T_2=(X{+}c_X)\circ(Y{+}c_Y)
\;\to\; T_3=\bigl((X{+}c_X)\circ(Y{+}c_Y)\bigr)+c_\ast,
$$

stabbing at the left child ($p_1$), the right child ($p_2$), then the root
($p_3$; here $g_i^{p_3}=\mathrm{id}$, so $S_i^{p_3}=t_i$ and it is ordinary
residual stabbing). The Theorem gives

$$\text{hits}(T_3)\ge \text{hits}(T_2)\ge \text{hits}(T_1)\ge \text{hits}(X\circ Y).$$

This is strictly stronger than $\text{hits}(\text{parent})\ge
\max(\text{hits}(\text{left}),\text{hits}(\text{right}))$: the whole operator is
a non-decreasing mutation, hit-by-hit at every intermediate stage.
$\blacksquare$

---

## 7. Why division cannot touch the proof

Division appears **only** inside $g_i^p$, and therefore only inside the
*definition* (2) of $S_i^p$. It plays two roles:

- **Denominator stab** (parent $= X/Y$, stabbing $Y$): $g_i^{Y}(u) = x_i/u$ with
  $g_i^Y(0)=\bot$, so
  $$S_i^Y = \{u : x_i/u \in t_i\},$$
  generically a **union of two rays**, and $\mathbb{R}\setminus\{0\}$ or
  $\varnothing$ in the degenerate $x_i=0$ cases.
- **Numerator stab**: $g_i^X(u)=u/y_i$, so
  $$S_i^X = t_i\cdot y_i \;\text{ if } y_i\neq 0, \qquad S_i^X=\varnothing \;\text{ if } y_i=0.$$

Every one of these is still just "some subset of $\mathbb{R}$." Lemmas 1–2 only
ever manipulate $u \in S_i^p$ as a **membership predicate** — connectedness,
boundedness, number of components, and sign-flip behavior never appear. So the
monotone chain holds through division *with equality of stab count and hits at
each step*. There is no weakened, division-adjusted invariant to fall back to.

### The two-ray sign-flip is not a counterexample

A common (incorrect) argument runs: an optimal $c_Y$ may place $Y_{\text{post}}$
on the *opposite ray* from $Y_{\text{pre}}$, flipping $X/Y$ in sign and magnitude
so that points hit by $T_1$ become misses in $T_2$; therefore
$\text{hits}(T_2)\ge\text{hits}(T_1)$ breaks.

The flaw is the jump from "these individual points are lost" to "the count
drops." **Set churn is not count decrease.** By Lemma 2, $c_Y = 0$ is always
feasible in the right stab and reproduces $T_1$ exactly, hitting precisely the
points $T_1$ hit. The stabber maximizes the total count, so it selects a
sign-flipping $c_Y$ **only when the points gained are at least as many as the
points lost**. The right stab therefore floors at $\text{hits}(T_1)$ regardless
of how violently individual points move across rays.

Two further notes on the same argument:

1. **"$+,-,\times$ degrade gracefully because their inverses are single
   intervals" is false.** A single-interval surrogate does not bound how far the
   chosen offset sits from zero; multiplication's stab can move the sibling
   arbitrarily far and across zero just as freely. If "large discontinuous move
   in the sibling" were the failure mechanism, $\times$ would fail too — it
   doesn't. The number of components of the surrogate is simply not a variable
   in the proof.
2. The two-ray observation *does* correctly identify division as the operator
   where the implementation is most likely to get the surrogate wrong — the
   right operator, for the wrong reason.

---

## 8. Consequence: the invariant is its own test oracle

Because Lemma 1 asserts $\text{hits}(T_k)=M_k$ **exactly**, you can assert it in
code. After each stab, recompute hits by direct evaluation and compare against
the claimed $M_k$. If you observe a genuine count drop
$\text{hits}(T_k) < \text{hits}(T_{k-1})$, the theorem guarantees the cause is a
violated assumption — never the mathematics. For a two-ray (division) step the
diagnosis narrows to:

| Observation | Diagnosis |
|---|---|
| claimed $M_k \ne$ evaluated $\text{hits}(T_k)$ | Counting/union bug in the stabber (Lemma 1 broken): a point can be covered by *either* ray; a point whose two rays both cover $c$ must not be double-counted. |
| claimed $M_k =$ evaluated $\text{hits}(T_k)$ but $M_k < \text{hits}(T_{k-1})$ | **Impossible** under a correct stabber, since $c=0$ is feasible and yields $\text{hits}(T_{k-1})$. The argmax search is not considering $c=0$ / is incomplete over ray-unions (A2 broken). |
| $\text{hits}(T_{k-1})$ logged now $\ne$ what the previous stab claimed | Tree state changed between steps — Jacobi/caching bug (A1's "against $T_{k-1}$" broken). |
| Surrogate silently became $(-\infty,\infty)$ | A union of rays was convex-hulled to a single interval: strict superset, inflates the intermediate count, drops at the next step. |
| Protected vs. true division mismatch | $g_i^p$ used to build the surrogate differs from $g_i^p$ used to evaluate, at $u=0$. Make both consistent. |
| $[0,0]$ returned when $y_i=0$, or line including $0$ when $x_i=0$ | A $\bot$-producing point was admitted into the preimage; it should be excluded ($S_i^p$ empty or punctured). |

Every branch lands on an implementation bug; none lands on "division inherently
violates the invariant."

**Forcing the question empirically.** Construct data where the optimal $c_Y$
flips several points across rays. Log, for the right stab: the claimed $M$, the
directly-evaluated $\text{hits}(T_2)$, and the directly-evaluated
$\text{hits}(T_1)$. You will see the lost points replaced by **at least as many
gained** ones, with the count holding. If instead the count drops on that
constructed case, the three-way log names the broken assumption immediately.

