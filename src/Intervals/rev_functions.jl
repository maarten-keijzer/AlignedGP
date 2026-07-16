# --- Reverse (preimage) functions ------------------------------------------

isperiodic(::Any) = false
#isperiodic(::typeof(sin)) = true

function maybeinterval(lo, hi) 
    (lo != Inf && hi != -Inf && lo <= hi) ? intervaltype(lo, hi) : nothing
end

function sqrt_rev(y::IntervalType) :: Union{Nothing, IntervalType}
    yl, yh = bounds(y)
    yh < 0.0 && return nothing
    yh == 0.0 && return intervaltype(0)

    xl = yl <= 0 ? 0.0 : nextfloat(yl * yl)
    xh = prevfloat(yh * yh)

    return maybeinterval(xl, xh)
end

function add_rev(z::IntervalType, y::Real) :: Union{Nothing, IntervalType}
    !isfinite(y) && return nothing

    zl, zh = bounds(z)
    xl, xh = zl - y, zh - y    

    # Only narrow a side when the round-to-nearest endpoint would push its
    # image outside z; exact shifts keep their (tight) endpoints.
    isfinite(xl) && inf(interval(xl) + y) < zl && (xl = nextfloat(xl))
    isfinite(xh) && sup(interval(xh) + y) > zh && (xh = prevfloat(xh))
    return maybeinterval(xl, xh)
end

umin_rev(z::IntervalType) :: IntervalType = -z 

function mul_rev(z::IntervalType, y::Real) :: Union{Nothing, IntervalType}
    !isfinite(y) && return nothing

    zl, zh = bounds(z)
    if y == 0.0
        (zl > 0 || zh < 0) && return nothing
        return intervaltype(-Inf, Inf)
    end

    xl, xh = minmax(zl / y, zh / y)

    xl = nextfloat(xl)
    xh = prevfloat(xh)
    return maybeinterval(xl, xh)
end

function inv_rev(y::IntervalType) :: Union{Nothing, IntervalType, Tuple{IntervalType, IntervalType}}
    yl, yh = bounds(y)
    if yl < 0 < yh
        intv1 = inv_rev(intervaltype(yl, 0))
        intv2 = inv_rev(intervaltype(0, yh))
        isnothing(intv1) && return intv2
        isnothing(intv2) && return intv1
        return tuple(intv1, intv2)
    elseif yl == yh == 0
        return nothing
    else
        inv_yh = yh == 0 ? -Inf : inv(yh)
        xl, xh = minmax(inv(yl), inv_yh)
        xl = nextfloat(xl)
        xh = prevfloat(xh)
        return maybeinterval(xl, xh)
    end
end

function exp_rev(y::IntervalType) :: Union{Nothing, IntervalType}
    yl, yh = bounds(y)
    yh <= 0.0 && return nothing

    xl = yl <= 0 ? nextfloat(-Inf) : nextfloat(log(yl))
    xh = prevfloat(log(yh))

     return maybeinterval(xl, xh)
end

function log_rev(y::IntervalType) :: Union{Nothing, IntervalType}
    yl, yh = bounds(y)

    xl = nextfloat(exp(yl))
    xh = prevfloat(exp(yh))

    xl == xh && isinf(xl) && return nothing
    return maybeinterval(xl, xh)
end

function sin_rev(y::IntervalType) :: Union{Nothing, IntervalType, Tuple{IntervalType, IntervalType}}
    yl, yh = bounds(y)

    cl = max(yl, -1.0)          # clamp target to sin's range [-1, 1]
    ch = min(yh,  1.0)
    cl > ch && return nothing   # target unreachable by sin

    PI  = intervaltype(π)       # must be a rigorous π enclosure (see note)
    Acl = asin(intervaltype(cl))

    if yh >= 1.0
        # peak covered: rising+falling merge into one arc spanning x = π/2,
        # so no spurious hole at the peak.
        lo = sup(Acl)            # asin(cl) inward (up)
        hi = inf(PI - Acl)       # π−asin(cl) inward (down)
        return lo <= hi ? intervaltype(lo, hi) : nothing
    end

    Ach = asin(intervaltype(ch))

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

const TWOPI = intervaltype(2) * intervaltype(π)
const FULL_CIRCLE = intervaltype(0.0, sup(TWOPI))

function sin_rev_circular(y::IntervalType) :: Union{Nothing, IntervalType, Tuple{IntervalType,IntervalType}}
    yl, yh = bounds(y)
    cl, ch = max(yl, -1.0), min(yh, 1.0)
    cl > ch && return nothing

    peak, trough = yh >= 1.0, yl <= -1.0
    peak && trough && return FULL_CIRCLE

    PI = intervaltype(π)

    if peak                                          # merged around π/2
        A = asin(intervaltype(cl))
        lo, hi = sup(A), inf(PI - A)
        return lo <= hi ? intervaltype(lo, hi) : nothing

    elseif trough                                    # merged around 3π/2, wraps
        B = asin(intervaltype(ch))
        lo, hi = sup(PI - B), inf(TWOPI + B)
        return lo <= hi ? intervaltype(lo, hi) : nothing

    else
        A, B = asin(intervaltype(cl)), asin(intervaltype(ch))
        rlo, rhi = sup(A), inf(B)                    # rising  ⊆ [−π/2, π/2]
        flo, fhi = sup(PI - B), inf(PI - A)          # falling ⊆ [π/2, 3π/2]
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
end
