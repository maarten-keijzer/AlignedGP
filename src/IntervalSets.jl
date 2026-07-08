struct CInterval
    lo::Float64
    hi::Float64
    CInterval(l, h) = l <= h ? new(l, h) : new(h, l)
end

const invalid_interval = CInterval(NaN, NaN)
_is_invalid(ci::CInterval) = isnan(ci.lo) || isnan(ci.hi)

# Narrow purposefully disabled for now
#narrow(lo, hi) = CInterval(lo, hi)
function narrow(lo, hi)
    lo, hi = Float64(lo), Float64(hi)
    lo > hi && ((lo, hi) = (hi, lo))
    lo_n = isinf(lo) ? lo : nextfloat(lo)
    hi_n = isinf(hi) ? hi : prevfloat(hi)
    lo_n > hi_n && return invalid_interval  # sub-ULP interval: no float strictly inside
    (hi_n - lo_n) <= 2eps(lo_n) && return invalid_interval # want two eps boundaries
    CInterval(lo_n, hi_n)
end

Base.first(ci::CInterval) = ci.lo
Base.last(ci::CInterval)  = ci.hi

Base.:(==)(a::CInterval, b::CInterval) = a.lo == b.lo && a.hi == b.hi
Base.isapprox(a::CInterval, b::CInterval; kwargs...) =
    isapprox(a.lo, b.lo; kwargs...) && isapprox(a.hi, b.hi; kwargs...)
Base.isapprox(a::AbstractVector{CInterval}, b::AbstractVector{CInterval}; kwargs...) =
    length(a) == length(b) && all(isapprox.(a, b; kwargs...))

# Treat CInterval as a scalar for broadcasting (like a number, not a collection).
Base.broadcastable(ci::CInterval) = Ref(ci)

Base.in(x::Real, ci::CInterval) = ci.lo <= x <= ci.hi

# Interval-interval arithmetic (narrowing)
Base.:+(ci1::CInterval, ci2::CInterval) = narrow(ci1.lo + ci2.lo, ci1.hi + ci2.hi)
Base.:-(ci1::CInterval, ci2::CInterval) = narrow(ci1.lo - ci2.hi, ci1.hi - ci2.lo)

# Scalar arithmetic (narrowing) — for broadcasting in leftinverse/rightinverse.
# Invalid sentinels pass through unchanged so they survive downstream filtering.
Base.:+(ci::CInterval, s::Real) = _is_invalid(ci) ? ci : narrow(ci.lo + s, ci.hi + s)
Base.:+(s::Real, ci::CInterval) = _is_invalid(ci) ? ci : narrow(s + ci.lo, s + ci.hi)
Base.:-(ci::CInterval, s::Real) = _is_invalid(ci) ? ci : narrow(ci.lo - s, ci.hi - s)
Base.:-(s::Real, ci::CInterval) = _is_invalid(ci) ? ci : narrow(s - ci.hi, s - ci.lo)

Base.:*(ci::CInterval, s::Real) = _is_invalid(ci) ? ci : narrow(ci.lo * s, ci.hi * s)
Base.:*(s::Real, ci::CInterval) = ci * s


const MAX_SUBINTERVALS = 16

# CIntervals stores one interval inline (_solo, n==1) to avoid heap allocation in the
# common single-interval case. A heap vector (_multi) is used only when n > 1.
struct CIntervals
    _solo::CInterval                           # valid only when n == 1
    _multi::Union{Nothing, Vector{CInterval}}  # valid only when n > 1
    n::Int
end

CIntervals() = CIntervals(CInterval(0.0, 0.0), nothing, 0)
CIntervals(ci::CInterval) = CIntervals(ci, nothing, 1)

function CIntervals(v::Vector{CInterval})
    n = length(v)
    if n > MAX_SUBINTERVALS
        v = v[1:MAX_SUBINTERVALS]
        n = MAX_SUBINTERVALS
    end
    n == 0 && return CIntervals()
    n == 1 && return CIntervals(v[1])
    CIntervals(CInterval(0.0, 0.0), v, n)
end

CIntervals(lo::Float64, hi::Float64) = CIntervals(CInterval(lo, hi))
CIntervals(t::NTuple{N, CInterval}) where {N} = CIntervals(collect(t))

Base.convert(::Type{CIntervals}, ci::CInterval) = CIntervals(ci)

function Base.iterate(c::CIntervals, state=1)
    state > c.n && return nothing
    ci = c.n == 1 ? c._solo : c._multi[state]
    return (ci, state + 1)
end

Base.length(c::CIntervals) = c.n
Base.isempty(c::CIntervals) = c.n == 0

function Base.in(x::Real, cis::CIntervals)
    cis.n == 0 && return false
    cis.n == 1 && return x in cis._solo
    any(ci -> x in ci, cis._multi)
end

function Base.:(==)(a::CIntervals, b::CIntervals)
    a.n != b.n && return false
    a.n == 0 && return true
    a.n == 1 && return a._solo == b._solo
    a._multi == b._multi
end

function Base.isapprox(a::CIntervals, b::CIntervals; kwargs...)
    a.n != b.n && return false
    a.n == 0 && return true
    a.n == 1 && return isapprox(a._solo, b._solo; kwargs...)
    all(isapprox(ai, bi; kwargs...) for (ai, bi) in zip(a._multi, b._multi))
end

# .items returns the internal Vector for backward compatibility.
# For n==1 this allocates a fresh 1-element vector; only used in tests / select_constant.
function Base.getproperty(cis::CIntervals, sym::Symbol)
    if sym === :items
        cis.n == 0 && return CInterval[]
        cis.n == 1 && return [cis._solo]
        return getfield(cis, :_multi)
    end
    return getfield(cis, sym)
end

function flatten(v::AbstractVector{CIntervals})
    result = CInterval[]
    for cis in v
        if cis.n == 1
            push!(result, cis._solo)
        elseif cis.n > 1
            append!(result, cis._multi)
        end
    end
    result
end

# Treat CIntervals as a scalar for broadcasting (like CInterval, not a collection).
Base.broadcastable(cis::CIntervals) = Ref(cis)

# Scalar arithmetic: n==1 fast path avoids heap allocation entirely.
function Base.:+(cis::CIntervals, s::Real)
    cis.n == 0 && return cis
    cis.n == 1 && return CIntervals(cis._solo + s)
    CIntervals([ci + s for ci in cis._multi])
end

function Base.:+(s::Real, cis::CIntervals)
    cis.n == 0 && return cis
    cis.n == 1 && return CIntervals(s + cis._solo)
    CIntervals([s + ci for ci in cis._multi])
end

function Base.:-(cis::CIntervals, s::Real)
    cis.n == 0 && return cis
    cis.n == 1 && return CIntervals(cis._solo - s)
    CIntervals([ci - s for ci in cis._multi])
end

function Base.:-(s::Real, cis::CIntervals)
    cis.n == 0 && return cis
    cis.n == 1 && return CIntervals(s - cis._solo)
    CIntervals([s - ci for ci in cis._multi])
end

function Base.:*(cis::CIntervals, s::Real)
    cis.n == 0 && return cis
    cis.n == 1 && return CIntervals(cis._solo * s)
    CIntervals([ci * s for ci in cis._multi])
end

Base.show(io::IO, ci::CInterval) = print(io, "<$(ci.lo),$(ci.hi)>")
function Base.show(io::IO, cis::CIntervals)
    if cis.n == 0
        print(io, "()")
    elseif cis.n == 1
        print(io, cis._solo)
    else
        print(io, "(")
        for i in eachindex(cis._multi)
            print(io, cis._multi[i])
            if i != cis.n 
                print(io, ", ")
            end
        end
        print(io, ")")
    end
end
