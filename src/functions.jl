
export unsafe 
unsafe(_) = false 
unsafe(::typeof(sqrt)) = true 
unsafe(::typeof(log)) = true


# evaluate: real × real → real (broadcasting over vectors)
evaluate(::typeof(+), x, y) = x .+ y
evaluate(::typeof(-), x, y) = x .- y
evaluate(::typeof(*), x, y) = x .* y
evaluate(::typeof(/), x, y) = x ./ y

# leftinverse(op, t, y)  → interval of valid left  args given target t and right arg y
# rightinverse(op, t, x) → interval of valid right args given target t and left  arg x
#
# Scalar arithmetic on CInterval narrows by one ULP (see IntervalSets.jl), so all
# derived interval propagation is conservatively inner-approximated throughout.

leftinverse(::typeof(+), t, y)  = t .- y
rightinverse(::typeof(+), t, x) = t .- x

leftinverse(::typeof(-), t, y)  = t .+ y
rightinverse(::typeof(-), t, x) = x .- t   # scalar - interval flips bounds

# * and / require sign-aware scaling.
# Invalid sentinels are preserved explicitly so they survive downstream filtering.

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

leftinverse(::typeof(*), t, y)  = _scale.(t, inv.(y))
rightinverse(::typeof(*), t, x) = _scale.(t, inv.(x))

leftinverse(::typeof(/), t, y)  = _scale.(t, y)
rightinverse(::typeof(/), t, x) = _div_into.(x, t)

# Unary functions: evaluate(op, x) and inverse(op, t)

evaluate(::typeof(identity), x) = identity.(x)

inverse(::typeof(identity), t::CInterval) = t
inverse(::typeof(identity), t) = inverse.(identity, t)

evaluate(::typeof(sqrt), x) = [v < 0 ? NaN : sqrt(v) for v in x]

# sqrt(x) = z ∈ [l,u]  →  x = z² ∈ [l²,u²] for z≥0.
# Negative targets are unreachable; return the invalid sentinel CInterval(Inf,Inf).
function inverse(::typeof(sqrt), t::CInterval)
    t.hi < 0.0 && return invalid_interval
    l = max(0.0, t.lo)
    narrow(l^2, t.hi^2)
end

inverse(::typeof(sqrt), t) = inverse.(sqrt, t)

evaluate(::typeof(exp), x) = exp.(x)

# exp(x) = z ∈ [l,u]  →  x = log(z) ∈ [log(l), log(u)] for z>0.
# exp is always positive, so u≤0 is unreachable. l≤0 imposes no lower bound (→ -Inf).
function inverse(::typeof(exp), t::CInterval)
    t.hi <= 0.0 && return invalid_interval
    l_log = t.lo <= 0.0 ? -Inf : log(t.lo)
    narrow(l_log, log(t.hi))
end

inverse(::typeof(exp), t) = inverse.(exp, t)

evaluate(::typeof(log), x) = [v <= 0 ? NaN : log(v) for v in x]

# log(x) = z ∈ [l,u]  →  x = exp(z) ∈ [exp(l), exp(u)], always valid.
inverse(::typeof(log), t::CInterval) = narrow(exp(t.lo), exp(t.hi))

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

evaluate(::typeof(cos), x) = [isfinite(v) ? cos(v) : NaN for v in x]
evaluate(::typeof(sin), x) = [isfinite(v) ? sin(v) : NaN for v in x]

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
