using Intervals 

const CI = Interval{Float64,Closed,Closed}
const CISet = IntervalSet{CI}
const CInterval = Union{CI, CISet}

function f(targets::Vector{CInterval}, evals::Vector{Float64})
    return targets .- evals
end

t = [ CI(1.0, 2.0), CI(3.0, 4.0) ]
ev = [1.5, 3.5]

f(t, ev)

a = CISet(CI(1.0, 2.0))

t2 = [a]

