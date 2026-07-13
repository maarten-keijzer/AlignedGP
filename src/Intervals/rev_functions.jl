# --- Reverse (preimage) functions ------------------------------------------

function sqrt_rev(y::IntervalType) :: Union{Nothing, IntervalType}
    yl, yh = bounds(y)
    yh < 0.0 && return nothing
    yh == 0.0 && return intervaltype(0)

    xl = yl < 0 ? 0.0 : nextfloat(yl * yl)
    xh = prevfloat(yh * yh)

    return xl <= xh ? intervaltype(xl, xh) : nothing
end

function add_rev(z::IntervalType, y::Real) :: Union{Nothing, IntervalType}
    !isfinite(y) && return nothing

    zl, zh = bounds(z)
    xl, xh = zl - y, zh - y    

    # Only narrow a side when the round-to-nearest endpoint would push its
    # image outside z; exact shifts keep their (tight) endpoints.
    isfinite(xl) && inf(interval(xl) + y) < zl && (xl = nextfloat(xl))
    isfinite(xh) && sup(interval(xh) + y) > zh && (xh = prevfloat(xh))
    return xl <= xh ? intervaltype(xl, xh) : nothing
end

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
    return xl <= xh ? intervaltype(xl, xh) : nothing
end

function inv_rev(y::IntervalType) :: Union{Nothing, IntervalType, Tuple{IntervalType, IntervalType}}
    yl, yh = bounds(y)
    if yl < 0 < yh
        intv1 = inv_rev(interval(yl, 0))
        intv2 = inv_rev(interval(0, yh))
        isempty_interval(intv1) && return intv2
        isempty_interval(intv2) && return intv1
        return tuple(intv1, intv2)
    elseif yl == yh == 0
        return nothing
    else
        inv_yh = yh == 0 ? -Inf : inv(yh)
        xl, xh = minmax(inv(yl), inv_yh)
        xl = nextfloat(xl)
        xh = prevfloat(xh)
        return xl <= xh ? intervaltype(xl, xh) : nothing
    end
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

function exp_rev(y::IntervalType) :: Union{Nothing, IntervalType}
    yl, yh = bounds(y)
    yh <= 0.0 && return nothing

    xl = yl <= 0 ? -Inf : nextfloat(log(yl))
    xh = prevfloat(log(yh))

    return xl <= xh ? intervaltype(xl, xh) : nothing
end

function log_rev(y::IntervalType) :: Union{Nothing, IntervalType}
    yl, yh = bounds(y)

    xl = nextfloat(exp(yl))
    xh = prevfloat(exp(yh))

    return xl <= xh ? intervaltype(xl, xh) : nothing
end
