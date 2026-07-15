module ReverseIntervals

using IntervalArithmetic

# Interface to either Interval or BareInterval
# const IntervalType = Interval{Float64}
# intervaltype(x::Real) = interval(x)
# intervaltype(lo::Real, hi::Real) = interval(lo, hi)

const IntervalType = BareInterval{Float64}
intervaltype(x::Real) = bareinterval(x)
intervaltype(lo::Real, hi::Real) = bareinterval(lo, hi)
Base.:+(iv::IntervalType, x::Real) = iv + bareinterval(x)
Base.:*(iv::IntervalType, x::Real) = iv * bareinterval(x) 

export IntervalType, intervaltype
export sqrt_rev, add_rev, mul_rev, inv_rev, sin_rev, exp_rev, log_rev, umin_rev
export IntervalVector, invert

include("rev_functions.jl")

struct IntervalVector 
    intervals::Vector{IntervalType}
    offsets::Vector{Int}
end
IntervalVector() = IntervalVector(IntervalType[], [1])
IntervalVector(iv::Vector{IntervalType}) = IntervalVector(iv, collect(1:length(iv)+1))

function IntervalVector(tuples::Tuple{T, T}...) where T <: Real 
    intvec = Vector{IntervalType}(undef, length(tuples))
    for (i, tup) in enumerate(tuples)
        intvec[i] = intervaltype(first(tup), last(tup))
    end
    IntervalVector(intvec)
end

caserange(off, i) = off[i]:(off[i+1]-1)      # count = off[i+1]-off[i], 0 when empty
caseview(iv, off, i) = @view iv[caserange(off, i)]   # no allocation
ncases(off) = length(off) - 1

Base.getindex(iv::IntervalVector, i::Int) = caseview(iv.intervals, iv.offsets, i)
Base.eachindex(iv::IntervalVector) = 1:ncases(iv.offsets)
Base.length(iv::IntervalVector) = ncases(iv.offsets)

function invert(iv::IntervalVector, rev_fun)
    intervals = IntervalType[]
    offsets = [1]
    for i in eachindex(iv)
        for out in rev_fun.(iv[i])
            if !isnothing(out)
                if out isa IntervalType
                    push!(intervals, out)
                else
                    for o in out
                        push!(intervals, o)
                    end
                end
            end
        end
        push!(offsets, length(intervals)+1)
    end
    return IntervalVector(intervals, offsets)
end

function invert(iv::IntervalVector, arg::Vector{<:Real}, rev_fun)
    intervals = IntervalType[]
    offsets = [1]
    for i in eachindex(iv)
        for out in rev_fun.(iv[i], arg[i])
            if !isnothing(out)
                push!(intervals, out)
            end
        end
        push!(offsets, length(intervals)+1)
    end
    return IntervalVector(intervals, offsets)
end

# Re-use add-rev as the inverse of addition is actually the sought for minus function
Base.:-(iv::IntervalVector, y::Vector{<:Real}) = invert(iv, y, add_rev)
function Base.:-(intervalvec::IntervalVector, y::Real) 
    yv = intervaltype(y)
    intv = [i - yv for i in intervalvec.intervals]
    IntervalVector(intv, intervalvec.offsets)
end

end # module Intervals