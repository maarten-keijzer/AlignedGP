using AlignedGP
using Test

setup = keijzer1([exp, log, sqrt, sin], [+, *, -, /])

for i in 1:10_000
    global tree
    tree, _ = valid_init(setup.symboltable, 15, setup.inputs)
    tree, _ = optimize(tree, setup.inputs, setup.interval_targets)
end


