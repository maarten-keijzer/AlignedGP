module AlignedGP

using Intervals


include("interval_alignment.jl")
include("functions.jl")
include("trees.jl")

export max_overlap_region, select_constant
export evaluate, leftinverse, rightinverse, inverse
export AddedValue, Var, Constant, BinaryNode, UnaryNode, Tree, insert

end # module AlignedGP
