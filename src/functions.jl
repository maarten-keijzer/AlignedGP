


# evaluate: real × real → real (broadcasting over vectors)
evaluate(::typeof(+), x, y) = x .+ y
evaluate(::typeof(-), x, y) = x .- y
evaluate(::typeof(*), x, y) = x .* y
evaluate(::typeof(/), x, y) = x ./ y

# evaluate: real → real , prevent DomainErrors from being thrown
evaluate(::typeof(-), x) = -x
evaluate(::typeof(inv), x) = inv.(x)

evaluate(::typeof(sqrt), x) = [v < 0 ? NaN : sqrt(v) for v in x]
evaluate(::typeof(log), x) = [v <= 0 ? NaN : log(v) for v in x]
evaluate(::typeof(exp), x) = exp.(x)
evaluate(::typeof(sin), x) = [isfinite(v) ? sin(v) : NaN for v in x]
evaluate(::typeof(cos), x) = [isfinite(v) ? cos(v) : NaN for v in x]


inverse(::typeof(+), t::IntervalVector, y::Vector{<:Real}) = invert(t, y, add_rev)
inverse(::typeof(+), x::Vector{<:Real}, t::IntervalVector) = invert(t, x, add_rev)
inverse(::typeof(*), t::IntervalVector, y::Vector{<:Real}) = invert(t, y, mul_rev)
inverse(::typeof(*), x::Vector{<:Real}, t::IntervalVector) = invert(t, x, mul_rev)

inverse(::typeof(-), t::IntervalVector, y::Vector{<:Real}) = invert(t, -y, add_rev)
inverse(::typeof(-), x::Vector{<:Real}, t::IntervalVector) = invert(invert(t, x, add_rev), umin_rev) # TODO: extra ulp
 
inverse(::typeof(/), t::IntervalVector, y::Vector{<:Real}) = invert(t, inv.(y), mul_rev) # TODO: extra ulp
inverse(::typeof(/), x::Vector{<:Real}, t::IntervalVector) = invert(invert(t, x, mul_rev), inv_rev)# TODO: extra ulp

inverse(::typeof(-), t::IntervalVector) = invert(t, umin_rev)
inverse(::typeof(inv), t::IntervalVector) = invert(t, inv_rev)
inverse(::typeof(sqrt), t::IntervalVector) = invert(t, sqrt_rev)
inverse(::typeof(log), t::IntervalVector) = invert(t, log_rev)
inverse(::typeof(exp), t::IntervalVector) = invert(t, exp_rev)
inverse(::typeof(sin), t::IntervalVector) = invert(t, sin_rev)
