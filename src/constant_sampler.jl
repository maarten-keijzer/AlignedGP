using IntervalArithmetic, Random, Statistics

# ------------------------------------------------------------------ config

"""
    ConstantSampler(; γ=1.0, α=0.5, far=1e8, thin_weight=0.0)

Draws a constant from a union of feasible intervals under a Cauchy(0, γ) prior
— i.e. the posterior, given the likelihood is flat on the feasible set.

* `γ`            prior scale: the magnitude below which you are indifferent.
                 Set per node; see [`scale_from_targets`](@ref).
* `α`            tempering exponent, applied to *interval selection* only.
                 1.0 = coherent posterior, 0.0 = uniform over intervals.
* `far`          |x|/γ past which `tan` is swapped for the exact Pareto(1) tail.
* `thin_weight`  selection weight for degenerate (point) intervals; 0 excludes.
"""
Base.@kwdef struct ConstantSampler
    γ::Float64           = 1.0
    α::Float64           = 0.5
    far::Float64         = 1e8
    thin_weight::Float64 = 0.0
    #ConstantSampler() = new()
end

#ConstantSampler(targets::AbstractVector{<:Real}; kw...) =
#    ConstantSampler(; γ = scale_from_targets(targets), kw...)


"""
    scale_from_targets(targets) -> Float64

Half-IQR of the finite targets. For Cauchy(0, γ) the quartiles sit exactly at
±γ, so half-IQR is the *matched* robust scale estimator — not a heuristic.
Falls back to median|t|, then 1.0.
"""
function scale_from_targets(targets)
    v = Float64[t for t in targets if isfinite(t)]
    isempty(v) && return 1.0
    γ = (quantile(v, 0.75) - quantile(v, 0.25)) / 2
    (isfinite(γ) && γ > 0) && return γ
    m = median(abs.(v))
    (isfinite(m) && m > 0) && return m
    return 1.0
end

# ---------------------------------------------------------------- internals

struct _Seg
    lo::Float64
    hi::Float64
    θlo::Float64
    θw::Float64
    mass::Float64   # untempered θ-measure
    tail::Int8      # 0 = trig, +1 = deep +tail, -1 = deep -tail
end

# inward rounding: guarantees the drawn value is inside the *true* interval
# even when bounds are BigFloat/Rational. No-ops for Float64.
_up(x::Float64)   = x
_down(x::Float64) = x
_up(x)   = Float64(x, RoundUp)
_down(x) = Float64(x, RoundDown)

function _segment(lo::Float64, hi::Float64, γ::Float64, far::Float64)
    if lo ≥ far * γ                        # Cauchy tail ≡ Pareto(1) to O(γ²/lo²)
        _Seg(lo, hi, 0.0, 0.0, γ * (inv(lo) - inv(hi)), Int8(1))
    elseif hi ≤ -far * γ
        _Seg(lo, hi, 0.0, 0.0, γ * (inv(-hi) - inv(-lo)), Int8(-1))
    else
        a, b = atan(lo, γ), atan(hi, γ)    # atan(±Inf, γ) == ±π/2 exactly
        _Seg(lo, hi, a, b - a, b - a, Int8(0))
    end
end

function _draw(rng::AbstractRNG, s::_Seg, γ::Float64)
    u = rand(rng)
    x = if s.tail == 0
        γ * tan(s.θlo + u * s.θw)
    elseif s.tail > 0
        inv(inv(s.lo) - u * (inv(s.lo) - inv(s.hi)))     # inv(Inf) === 0.0
    else
        -inv(inv(-s.hi) - u * (inv(-s.hi) - inv(-s.lo)))
    end
    isfinite(x) || (x = isfinite(s.lo) ? s.lo : s.hi)
    return clamp(x, s.lo, s.hi)
end

_each(X::Interval) = (X,)
_each(Xs)          = Xs

# -------------------------------------------------------------------- draw

"""
    draw_constant(rng, s::ConstantSampler, Xs) -> Float64

`Xs` is an `Interval` or an iterable of disjoint `Interval`s. Unbounded,
half-bounded, bounded and zero-crossing segments all take the same path.
Single pass, no allocation (weighted reservoir selection).
"""
function draw_constant(rng::AbstractRNG, s::ConstantSampler, Xs)
    γ = s.γ
    @assert isfinite(γ) && γ > 0 "prior scale γ must be finite and positive"

    W        = 0.0
    chosen   = _Seg(NaN, NaN, 0.0, 0.0, 0.0, Int8(0))
    fallback = NaN

    for X in _each(Xs)
        isempty_interval(X) && continue
        lo, hi = _up(inf(X)), _down(sup(X))
        (isnan(lo) || isnan(hi) || hi < lo) && continue
        isnan(fallback) && (fallback = isfinite(lo) ? lo : (isfinite(hi) ? hi : 0.0))

        seg = _segment(lo, hi, γ, s.far)
        w   = seg.mass > 0 ? seg.mass^s.α : s.thin_weight
        (isfinite(w) && w > 0) || continue

        W += w
        rand(rng) * W < w && (chosen = seg)
    end

    W > 0 || return fallback          # empty, or all-thin with thin_weight=0
    return _draw(rng, chosen, γ)
end

draw_constant(s::ConstantSampler, Xs) = draw_constant(Random.default_rng(), s, Xs)

intv = [interval(0.1, 0.3), interval(5, Inf)]
s = ConstantSampler()
draw_constant(s, intv)
