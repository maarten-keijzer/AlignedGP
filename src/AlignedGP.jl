module AlignedGP

using Random
using Statistics: median, mean
using PythonCall

include("Intervals/Intervals.jl")

using .ReverseIntervals

include("problem_setup.jl")
include("efforts.jl")
include("interval_alignment.jl")
include("functions.jl")
include("nodes.jl")
include("tree.jl")
include("simplify.jl")
include("two_band.jl")
include("insert_with_alignment.jl")
include("linear_scaling.jl")
include("init.jl")
include("variation.jl")
include("lexicase_algo.jl")
include("front.jl")

export max_overlap_region, select_constant
export evaluate, inverse
export compute_hits, hitvector
export two_band_hits, two_band_loss, two_band_score, two_band_lexicase, two_band_lexicase_pool
export residual_eps_lexicase, select_parent
export retarget
export AddedValue, Node, Var, Constant, BinaryNode, UnaryNode, Tree, insert, insert_with_alignment, complexity
export model_evaluations, linear_scale, linear_scale_l1
export SymbolTable,init, initstrata, iteratestrata!
export valid_init, initial_fit, optimize
export ProblemSetup, GPParams, simplify
export min_tolerance, pathlen_complexity

export sin_rev_circular
export coordinate_descent

export Front, add_to_front!, merge_with_front!, complexities_front, errors_front


end # module AlignedGP
