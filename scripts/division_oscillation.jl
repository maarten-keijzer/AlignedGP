
using AlignedGP, Random, Printf

const CI = CInterval

# Access unexported internals via module prefix
const _sl       = AlignedGP._surrogate_left
const _sr       = AlignedGP._surrogate_right
const _cav      = AlignedGP.compute_added_value

inputs = [[-0.1, 0.0, 0.2, 0.8]]
targets = [
    CIntervals(0.1, 0.2), 
    CIntervals(-0.1, 0.1),
    CIntervals(0.1, 0.3),
    CIntervals(-0.4, 0.2)
]

c1 = ones(size(inputs[1])) * 0.1 



