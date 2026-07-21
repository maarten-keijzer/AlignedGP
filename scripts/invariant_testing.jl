using AlignedGP
using Test

setup = keijzer4(unaries=[exp, log, sqrt], binaries=[+,*,-,/]);

# Re-fit exactly ONE additive constant along a random root->i path: descend to the
# subtree at index i and re-align its top added-value against the surrogate target,
# leaving every other constant on the path untouched (recursive_stabbing = false).
optimize_random_path(root, i, inputs, targets) =
    insert_with_alignment(root, root[i], i, 1, inputs, targets, false)[1:2]

for i in 1:1000
    global tree
    tree, evals = valid_init(setup.symboltable, 25, setup.inputs)
    tree, _ = optimize(tree, setup.inputs, setup.interval_targets)
    hits = compute_hits(evaluate(tree, setup.inputs), setup.interval_targets)

    # for _ in 1:100
    #     tree, evaluations = optimize_random_path(tree, rand(1:length(tree)), setup.inputs, setup.interval_targets)
    #     newhits = compute_hits(evaluations, setup.interval_targets)
    #     #@show hits, newhits
    #     hits = newhits
    # end
end


