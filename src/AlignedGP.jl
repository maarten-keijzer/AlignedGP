module AlignedGP

using Random
using Statistics: median, mean
using PythonCall

include("Intervals/Intervals.jl")
include("IntervalSets.jl")
include("problem_setup.jl")
include("efforts.jl")
include("interval_alignment.jl")
include("functions.jl")
include("nodes.jl")
include("tree.jl")
include("insert_with_alignment.jl")
include("linear_scaling.jl")
include("init.jl")
include("variation.jl")
include("lexicase_algo.jl")

export CInterval, CIntervals, narrow, _is_invalid, invalid_interval, flatten
export max_overlap_region, select_constant
export evaluate, leftinverse, rightinverse, inverse, _scale, _div_into
export AddedValue, Node, Var, Constant, BinaryNode, UnaryNode, Tree, insert, insert_with_alignment, complexity
export scaled_evaluate, linear_scale, linear_scale_l1
export SymbolTable,init, initstrata, iteratestrata!, print_report
export valid_init, initial_fit, optimize
export ProblemSetup, GPParams

export simple_regression, coordinate_descent

end # module AlignedGP
