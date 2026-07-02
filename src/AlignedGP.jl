module AlignedGP

using Intervals

include("interval_alignment.jl")
include("functions.jl")
include("trees.jl")
include("insert_with_alignment.jl")
include("linear_scaling.jl")
include("init.jl")
include("problem_setup.jl")
include("variation.jl")
include("fitness.jl")

export max_overlap_region, select_constant
export evaluate, leftinverse, rightinverse, inverse
export AddedValue, Var, Constant, BinaryNode, UnaryNode, Tree, insert, insert_with_alignment, complexity
export scaled_evaluate, linear_scale
export SymbolTable,init, initstrata, iteratestrata!, print_report

export simple_regression, coordinate_descent

end # module AlignedGP
