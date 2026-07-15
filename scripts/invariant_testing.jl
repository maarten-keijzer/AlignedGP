using AlignedGP
using Test

setup = keijzer1([exp], [+,*]);

for i in 1:1000
    global tree
    tree, _ = valid_init(setup.symboltable, 15, setup.inputs)
    tree, _ = optimize(tree, setup.inputs, setup.interval_targets)
end


