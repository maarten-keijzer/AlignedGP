using Intervals

# evaluate: real × real → real (broadcasting over vectors)
evaluate(::typeof(+), x, y) = x .+ y
evaluate(::typeof(-), x, y) = x .- y
evaluate(::typeof(*), x, y) = x .* y
evaluate(::typeof(/), x, y) = x ./ y

# Intervals.jl defines Interval ± scalar and scalar - Interval (with bound flip),
# so + and - inverses reduce to those operations and broadcast naturally.
#
# leftinverse(op, t, y)  → interval of valid left  args given target t and right arg y
# rightinverse(op, t, x) → interval of valid right args given target t and left  arg x

leftinverse(::typeof(+), t, y)  = t .- y          # [l-y, u-y]
rightinverse(::typeof(+), t, x) = t .- x          # [l-x, u-x]

leftinverse(::typeof(-), t, y)  = t .+ y          # [l+y, u+y]
rightinverse(::typeof(-), t, x) = x .- t          # [x-u, x-l]  (scalar-Interval flips bounds)

# * and / require sign-aware scaling; not defined in Intervals.jl.
# Using min/max handles sign flip without a branch, and lets float arithmetic
# propagate Inf/NaN for zero divisors (deliberate — caller filters degenerate events).

_scale(t::Interval{T,Closed,Closed}, s::Real) where {T} =
    (a = first(t) * s; b = last(t) * s; Interval(min(a, b), max(a, b)))

_div_into(x::Real, t::Interval{T,Closed,Closed}) where {T} =
    (a = x / first(t); b = x / last(t); Interval(min(a, b), max(a, b)))

leftinverse(::typeof(*), t, y)  = _scale.(t, inv.(y))
rightinverse(::typeof(*), t, x) = _scale.(t, inv.(x))

leftinverse(::typeof(/), t, y)  = _scale.(t, y)
rightinverse(::typeof(/), t, x) = _div_into.(x, t)

# Unary functions: evaluate(op, x) and inverse(op, t)

evaluate(::typeof(sqrt), x) = sqrt.(x)

# sqrt(x) = z ∈ [l,u]  →  x = z² ∈ [l²,u²] for z≥0.
# Negative targets are unreachable; return the invalid sentinel Interval(typemax,typemax).
function inverse(::typeof(sqrt), t::Interval{T,Closed,Closed}) where {T}
    last(t) < zero(T) && return Interval(typemax(T), typemax(T))
    l = max(zero(T), first(t))
    Interval(l^2, last(t)^2)
end

inverse(::typeof(sqrt), t) = inverse.(sqrt, t)

evaluate(::typeof(exp), x) = exp.(x)

# exp(x) = z ∈ [l,u]  →  x = log(z) ∈ [log(l), log(u)] for z>0.
# exp is always positive, so u≤0 is unreachable. l≤0 imposes no lower bound (→ -Inf).
function inverse(::typeof(exp), t::Interval{T,Closed,Closed}) where {T}
    last(t) <= zero(T) && return Interval(typemax(T), typemax(T))
    l_log = first(t) <= zero(T) ? oftype(last(t), -Inf) : log(first(t))
    Interval(l_log, log(last(t)))
end

inverse(::typeof(exp), t) = inverse.(exp, t)

evaluate(::typeof(log), x) = log.(x)   # DomainError for x≤0 caught at Tree level

# log(x) = z ∈ [l,u]  →  x = exp(z) ∈ [exp(l), exp(u)], always valid.
inverse(::typeof(log), t::Interval{T,Closed,Closed}) where {T} = Interval(exp(first(t)), exp(last(t)))

inverse(::typeof(log), t) = inverse.(log, t)
