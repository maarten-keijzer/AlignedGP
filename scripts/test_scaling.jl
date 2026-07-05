using AlignedGP
using Test
setup = keijzer1()

sub, _ = valid_init(setup.symboltable, 51, setup.inputs)

tree = sub
ev = evaluate(tree, setup.inputs)
hits = sum(ev[i] ∈ setup.interval_targets[i] for i in eachindex(ev))

# Optimize full tree
for i in 1:length(tree)
    global tree
    tree, _ = insert_with_alignment(tree, tree[i], i, 1, setup.inputs, setup.interval_targets, true)
end

ev = evaluate(tree, setup.inputs)
hits = sum(ev[i] ∈ setup.interval_targets[i] for i in eachindex(ev))

for i in 1:length(tree)
    add = tree[i].addition
    @show add.hits, tree[i]
    @test add.hits == hits
end
