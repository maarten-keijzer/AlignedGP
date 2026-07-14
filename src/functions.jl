


# evaluate: real × real → real (broadcasting over vectors)
evaluate(::typeof(+), x, y) = x .+ y
evaluate(::typeof(-), x, y) = x .- y
evaluate(::typeof(*), x, y) = x .* y
evaluate(::typeof(/), x, y) = x ./ y

# evaluate: real → real , prevent DomainErrors from being thrown
evaluate(::typeof(sqrt), x) = [v < 0 ? NaN : sqrt(v) for v in x]
evaluate(::typeof(log), x) = [v <= 0 ? NaN : log(v) for v in x]
evaluate(::typeof(exp), x) = exp.(x)
evaluate(::typeof(sin), x) = [isfinite(v) ? sin(v) : NaN for v in x]
evaluate(::typeof(cos), x) = [isfinite(v) ? cos(v) : NaN for v in x]


leftinverse(::typeof(+), t::IntervalType, y::Real) = add_rev(t, y)
rightinverse(::typeof(+), t::IntervalType, x::Real) = add_rev(t, x)

leftinverse(::typeof(/), t::IntervalType, y::Real) = mul_rev(t, inv(y))
# t = x / y -- t/x = inv(y) -- x/t = y -- y = x * inv(t)
function rightinverse(::typeof(/), t::IntervalType, y::Real)
    invt = inv_rev(t)
    isnothing(invt) && return nothing 
    isa(invt, IntervalType) && return mul_rev(invt, y)
    return mul_rev(first(invt), y), mul_rev(last(invt), y)
end


export unsafe 
unsafe(_) = false 
unsafe(::typeof(sqrt)) = true 
unsafe(::typeof(log)) = true


# leftinverse(op, t, y)  → interval of valid left  args given target t and right arg y
# rightinverse(op, t, x) → interval of valid right args given target t and left  arg x
#
# The surrogate must be a FORWARD-inner approximation: every u in the returned interval
# must satisfy `op(·) ∈ t` after floating-point rounding. Narrowing by one ULP in the
# surrogate's own coordinates is NOT enough — under scale mismatch (e.g. a near-zero
# multiplicative sibling, or a huge additive sibling) the forward op rounds at a coarser
# scale than the surrogate, so a "surrogate hit" can miss the target. That breaks the
# monotone-chain invariant (see specs/invariant_proof.md §2, assumption A1). We therefore
# round-trip-verify every surrogate against the forward op and drop points we cannot honor.

# `_guard(a, b, t, fwd)` — narrow [a,b] and keep it only if the forward image of BOTH
# endpoints lands in `t`. `fwd` is monotone on each candidate piece, so endpoint
# containment implies the whole (narrowed) interval maps into `t`. Overflow / extreme
# scaling that pushes an endpoint's image outside `t` yields the invalid sentinel — the
# point simply does not vote, which is the correct inner-approximation behavior.
function _guard(a::Float64, b::Float64, t::CInterval, fwd)
    (isnan(a) || isnan(b)) && return invalid_interval
    # Cap infinite bounds to the largest finite magnitude (by sign, via clamp). An infinite
    # child eval forwards to ±Inf or NaN and can never land in a finite target, so it must
    # not be admitted; the true preimage of a finite target is a set of finite reals, and
    # floatmax is the most extreme finite child that can occur. Capping (rather than changing
    # hit-counting) leaves the parent side untouched, and a guarded child hit then provably
    # forwards into the target. A degenerate [Inf,Inf] preimage (e.g. log of a target above
    # ~709) clamps to [floatmax,floatmax], which narrow rejects as invalid. clamp also keeps
    # `fwd` from being evaluated at an infinite argument outside its domain (e.g. log(-Inf)).
    fm = floatmax(Float64)
    s = narrow(clamp(min(a, b), -fm, fm), clamp(max(a, b), -fm, fm))
    _is_invalid(s) && return s
    (fwd(s.lo) in t && fwd(s.hi) in t) ? s : invalid_interval
end

# Per-piece guarded preimages for the multiplicative operators. Division of the target
# bounds by the sibling is used directly (no `inv(x)` reciprocal), avoiding the double
# rounding and the overflow-to-(-Inf,Inf) that produced phantom hits.
_mul_pre(t::CInterval, x::Real) =                       # u with x*u ∈ t   (* is commutative)
    # No iszero(x) special case: with x==0 the bounds become ±Inf and the guard rejects
    # them because 0*(±Inf) = NaN ∉ t. A near-/exactly-zero sibling therefore certifies no
    # child value (correct: 0*u is a constant, an additive child stab cannot change it, and
    # an infinite child eval would forward to NaN). Undercounting here is invariant-safe.
    _guard(t.lo / x, t.hi / x, t, u -> x * u)

_divl_pre(t::CInterval, y::Real) =                      # L with L/y ∈ t
    iszero(y) ? invalid_interval :                      # L/0 = ⊥ ∉ t
    _guard(t.lo * y, t.hi * y, t, L -> L / y)

_divr_pre(t::CInterval, x::Real) =                      # R with x/R ∈ t  (non-straddling t)
    iszero(x) ? (0.0 in t ? CInterval(-Inf, Inf) : invalid_interval) :
    t.lo < 0.0 < t.hi ? invalid_interval :              # two rays; handled at the CIntervals level
    _guard(x / t.lo, x / t.hi, t, R -> x / R)

# Per-piece left/right inverse of a single target interval, forward-guarded.
_linv(op, t::CInterval, y::Real) =
    _is_invalid(t)  ? invalid_interval :
    op === (+)      ? _guard(t.lo - y, t.hi - y, t, L -> L + y) :
    op === (-)      ? _guard(t.lo + y, t.hi + y, t, L -> L - y) :
    op === (*)      ? _mul_pre(t, y)  :
                      _divl_pre(t, y)                   # op === (/)

_rinv(op, t::CInterval, x::Real) =
    _is_invalid(t)  ? invalid_interval :
    op === (+)      ? _guard(t.lo - x, t.hi - x, t, R -> x + R) :
    op === (-)      ? _guard(x - t.hi, x - t.lo, t, R -> x - R) :
    op === (*)      ? _mul_pre(t, x)  :
                      _divr_pre(t, x)                   # op === (/)

# Lift over a multi-piece target. `x/R ∈ t` with a sign-crossing piece splits into two
# rays, so `/` right-inverse gets its own CIntervals method below.
_linv(op, t::CIntervals, y::Real) = CIntervals(filter(!_is_invalid, [_linv(op, ci, y) for ci in t]))
_rinv(op, t::CIntervals, x::Real) = CIntervals(filter(!_is_invalid, [_rinv(op, ci, x) for ci in t]))

function _rinv(::typeof(/), t::CIntervals, x::Real)     # R with x/R ∈ t, ray-aware
    t.n == 0 && return t
    iszero(x) && return 0.0 in t ? CIntervals([narrow(-Inf, 0.0), narrow(0.0, Inf)]) : CIntervals()
    results = CInterval[]
    for ci in t
        _is_invalid(ci) && continue
        if ci.lo < 0.0 < ci.hi
            # x/R ∈ ci ⟺ R ∈ (-∞, x/lo] ∪ [x/hi, +∞); guard the finite ray endpoints.
            a, b = x / ci.lo, x / ci.hi
            lo_end, hi_end = min(a, b), max(a, b)
            r1 = _guard(-Inf, lo_end, ci, R -> x / R); _is_invalid(r1) || push!(results, r1)
            r2 = _guard(hi_end, Inf, ci, R -> x / R);  _is_invalid(r2) || push!(results, r2)
        else
            r = _divr_pre(ci, x)
            _is_invalid(r) || push!(results, r)
        end
    end
    CIntervals(results)
end

# Public API — broadcast the per-piece builders, preserving the previous polymorphism
# over CInterval, CIntervals, and vectors thereof.
leftinverse(op, t, y)  = _linv.(op, t, y)
rightinverse(op, t, x) = _rinv.(op, t, x)

# * and / target-bound scaling helpers. Retained for their invalid-sentinel contract
# (see test_functions.jl); the inverse path above uses the forward-guarded builders.
_scale(t::CInterval, s::Real) =
    _is_invalid(t) || iszero(s) ? invalid_interval :
    isinf(s) ? (t.lo <= 0.0 <= t.hi ? CInterval(-Inf, Inf) : invalid_interval) :
    let a = t.lo * s, b = t.hi * s
        (isnan(a) || isnan(b)) ? invalid_interval : narrow(min(a, b), max(a, b))
    end

_div_into(x::Real, t::CInterval) =
    _is_invalid(t) ? invalid_interval :
    iszero(x) ? (0.0 in t ? CInterval(-Inf, Inf) : invalid_interval) :
    t.lo < 0.0 < t.hi ? invalid_interval :     # two disjoint rays; handled correctly at CIntervals level
    let a = x / t.lo, b = x / t.hi
        (isnan(a) || isnan(b)) ? invalid_interval : narrow(min(a, b), max(a, b))
    end

# Unary functions: evaluate(op, x) and inverse(op, t)


# sqrt(x) = z ∈ [l,u]  →  x = z² ∈ [l²,u²] for z≥0. Forward-guarded so an infinite or
# overflowing preimage bound cannot admit a child eval whose sqrt escapes the target.
function inverse(::typeof(sqrt), t::CInterval)
    t.hi < 0.0 && return invalid_interval
    l = max(0.0, t.lo)
    _guard(l^2, t.hi^2, t, sqrt)
end

inverse(::typeof(sqrt), t) = inverse.(sqrt, t)

# exp(x) = z ∈ [l,u]  →  x = log(z) ∈ [log(l), log(u)] for z>0.
# exp is always positive, so u≤0 is unreachable. l≤0 imposes no lower bound (→ -Inf).
function inverse(::typeof(exp), t::CInterval)
    t.hi <= 0.0 && return invalid_interval
    l_log = t.lo <= 0.0 ? -Inf : log(t.lo)
    _guard(l_log, log(t.hi), t, exp)
end

inverse(::typeof(exp), t) = inverse.(exp, t)

# log(x) = z ∈ [l,u]  →  x = exp(z) ∈ [exp(l), exp(u)], always valid. Forward-guarded so an
# overflowing exp(u)=Inf bound cannot admit an infinite child eval (log(Inf)=Inf ∉ t).
inverse(::typeof(log), t::CInterval) = _guard(exp(t.lo), exp(t.hi), t, log)

inverse(::typeof(log), t) = inverse.(log, t)

# CIntervals lifts: distribute over sub-intervals, filter invalids.
# Empty input stays empty (the error state propagates).

function _scale(cis::CIntervals, s::Real)
    cis.n == 0 && return cis
    if cis.n == 1
        r = _scale(cis._solo, s)
        return _is_invalid(r) ? CIntervals() : CIntervals(r)
    end
    CIntervals(filter(!_is_invalid, [_scale(ci, s) for ci in cis._multi]))
end

function _div_into(x::Real, cis::CIntervals)
    cis.n == 0 && return cis
    if iszero(x)
        # x/y = 0 for all y ≠ 0; y = 0 gives ⊥ which must be excluded.
        # Valid y range is ℝ\{0} iff 0 ∈ cis, else empty.
        return 0.0 in cis ? CIntervals([narrow(-Inf, 0.0), narrow(0.0, Inf)]) : CIntervals()
    end
    if cis.n == 1
        t = cis._solo
        _is_invalid(t) && return CIntervals()
        if !iszero(x) && t.lo < 0.0 < t.hi
            # Pre-image is two disjoint rays: (-∞, x/lo] ∪ [x/hi, +∞)
            a, b = x / t.lo, x / t.hi
            return CIntervals([CInterval(-Inf, min(a, b)), CInterval(max(a, b), Inf)])
        end
        r = _div_into(x, t)
        return _is_invalid(r) ? CIntervals() : CIntervals(r)
    end
    results = CInterval[]
    for ci in cis._multi
        _is_invalid(ci) && continue
        if !iszero(x) && ci.lo < 0.0 < ci.hi
            a, b = x / ci.lo, x / ci.hi
            push!(results, CInterval(-Inf, min(a, b)))
            push!(results, CInterval(max(a, b), Inf))
        else
            r = _div_into(x, ci)
            _is_invalid(r) || push!(results, r)
        end
    end
    CIntervals(results)
end

inverse(::typeof(identity), cis::CIntervals) = cis

function _lift_inverse(fn, cis::CIntervals)
    cis.n == 0 && return cis
    if cis.n == 1
        r = inverse(fn, cis._solo)
        return _is_invalid(r) ? CIntervals() : CIntervals(r)
    end
    CIntervals(filter(!_is_invalid, [inverse(fn, ci) for ci in cis._multi]))
end

inverse(::typeof(sqrt), cis::CIntervals) = _lift_inverse(sqrt, cis)
inverse(::typeof(exp),  cis::CIntervals) = _lift_inverse(exp,  cis)
inverse(::typeof(log),  cis::CIntervals) = _lift_inverse(log,  cis)


# cos is periodic with two monotone branches per period in [−π, π]:
#   decreasing [0,  π]: cos(x) ∈ [lo,hi] ↔ x ∈ [acos(hi), acos(lo)]
#   increasing [−π, 0]: cos(x) ∈ [lo,hi] ↔ x ∈ [−acos(lo), −acos(hi)]
# Range of cos is [−1, 1]; targets outside are clamped; fully out-of-range → empty.
function inverse(::typeof(cos), t::CInterval)
    lo_c = max(-1.0, t.lo)
    hi_c = min(1.0, t.hi)
    lo_c > hi_c && return CIntervals()
    CIntervals([narrow(acos(hi_c), acos(lo_c)), narrow(-acos(lo_c), -acos(hi_c))])
end

inverse(::typeof(cos), t) = inverse.(cos, t)

# sin is periodic with two monotone branches per period.
# Three cases based on how the target overlaps [−1, 1]:
#
#   Fully outside [−1,1]: empty preimage.
#
#   Upper-clamped (t.hi > 1): the peak (sin = 1) at π/2 is included, so the two
#   branches merge into one connected arc [asin(lo_c), π − asin(lo_c)].
#
#   Lower-clamped (t.lo < −1): the trough (sin = −1) at 3π/2 is included. On the
#   circle the two branches connect through the trough, forming one arc that wraps
#   across the period boundary. Represented on the real line as the single interval
#   [π − asin(hi_c), asin(hi_c) + 2π].
#
#   Fully inside [−1,1]: two genuinely disjoint arcs, one rising, one falling.
function inverse(::typeof(sin), t::CInterval)
    lo_c = max(-1.0, t.lo)
    hi_c = min(1.0, t.hi)
    lo_c > hi_c && return CIntervals()
    if t.hi > 1.0
        # Peak included: rising and falling arcs merge at π/2.
        return CIntervals(narrow(asin(lo_c), π - asin(lo_c)))
    elseif t.lo < -1.0
        # Trough included: arcs connect through 3π/2, wrapped as a single real interval.
        return CIntervals(narrow(π - asin(hi_c), asin(hi_c) + 2π))
    else
        return CIntervals([narrow(asin(lo_c), asin(hi_c)), narrow(π - asin(hi_c), π - asin(lo_c))])
    end
end

inverse(::typeof(sin), t) = inverse.(sin, t)

# Lift for functions whose CInterval-level inverse returns CIntervals (multi-branch).
# Flattens the sub-interval results from each element into one CIntervals.
function _lift_inverse_multi(fn, cis::CIntervals)
    cis.n == 0 && return cis
    cis.n == 1 && return inverse(fn, cis._solo)  # already a CIntervals (2 branches for cos/sin)
    results = CInterval[]
    for ci in cis._multi
        for sub in inverse(fn, ci)
            _is_invalid(sub) || push!(results, sub)
        end
    end
    CIntervals(results)
end

inverse(::typeof(cos), cis::CIntervals) = _lift_inverse_multi(cos, cis)
inverse(::typeof(sin), cis::CIntervals) = _lift_inverse_multi(sin, cis)

# Circular preimage helpers for periodic functions (sin/cos).
# Returns arcs in [0,2π) as Vector{Tuple{Float64,Float64}}.

is_periodic(fn) = fn === sin || fn === cos

# Reduce arc [start,stop] (real line, length ≤ 2π) to [0,2π), splitting at the seam if needed.
function _emit(start::Float64, stop::Float64)
    s = mod(start, 2π)
    e = s + (stop - start)
    e <= 2π ? Tuple{Float64,Float64}[(s, e)] : Tuple{Float64,Float64}[(s, 2π), (0.0, e - 2π)]
end

# Merge adjacent/overlapping arcs (assumes input already sorted by start).
function _coalesce_arcs(arcs::Vector{Tuple{Float64,Float64}})
    isempty(arcs) && return arcs
    sort!(arcs; by=first)
    out = Tuple{Float64,Float64}[first(arcs)]
    for (s, e) in Iterators.drop(arcs, 1)
        rs, re = last(out)
        s <= re ? (out[end] = (rs, max(re, e))) : push!(out, (s, e))
    end
    out
end

# One-period preimage of target ⊆ [-1,1] under sin, as arcs in [0,2π).
# Implements the algo sketch: rising branch emit(aL,aH) + falling branch emit(π-aH,π-aL).
function _preimage_circular(::typeof(sin), t::CInterval)
    lo_c = max(-1.0, t.lo); hi_c = min(1.0, t.hi)
    lo_c > hi_c && return Tuple{Float64,Float64}[]
    lo_c == -1.0 && hi_c == 1.0 && return Tuple{Float64,Float64}[(0.0, 2π)]
    aL = asin(lo_c); aH = asin(hi_c)   # both in [-π/2, π/2], aL ≤ aH
    _coalesce_arcs(vcat(_emit(aL, aH), _emit(π - aH, π - aL)))
end

# One-period preimage for cos, arcs in [0,2π).
function _preimage_circular(::typeof(cos), t::CInterval)
    lo_c = max(-1.0, t.lo); hi_c = min(1.0, t.hi)
    lo_c > hi_c && return Tuple{Float64,Float64}[]
    lo_c == -1.0 && hi_c == 1.0 && return Tuple{Float64,Float64}[(0.0, 2π)]
    aL = acos(hi_c); aH = acos(lo_c)   # aL ≤ aH in [0, π]
    _coalesce_arcs(vcat(_emit(aL, aH), _emit(2π - aH, 2π - aL)))
end
