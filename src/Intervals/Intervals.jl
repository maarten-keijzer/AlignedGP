module RevIntervals

using IntervalArithmetic
export issubset_interval, in_interval, isequal_interval, sup, inf # for testing

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
export sqrt_rev, add_rev, mul_rev, inv_rev, sin_rev, exp_rev, log_rev
export Intervals, invert

include("rev_functions.jl")

struct Intervals 
    intervals::Vector{IntervalType}
    offsets::Vector{Int}
end

Intervals(iv::Vector{IntervalType}) = Intervals(iv, collect(1:length(iv)+1))

caserange(off, i) = off[i]:(off[i+1]-1)      # count = off[i+1]-off[i], 0 when empty
caseview(iv, off, i) = @view iv[caserange(off, i)]   # no allocation
ncases(off) = length(off) - 1

Base.getindex(iv::Intervals, i::Int) = caseview(iv.intervals, iv.offsets, i)
Base.eachindex(iv::Intervals) = 1:ncases(iv.offsets)
Base.length(iv::Intervals) = ncases(iv.offsets)

function invert(iv::Intervals, rev_fun)
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
    return Intervals(intervals, offsets)
end

function invert(iv::Intervals, arg::Vector{<:Real}, rev_fun)
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
    return Intervals(intervals, offsets)
end


end # module Intervals