using Random

function draw_subtree(node::Node, rng=Random.GLOBAL_RNG)
    if length(node) == 1
        return node, 1
    end

    want_internal = rand(rng) < 0.9

    while true 
        index = rand(rng, 1:length(node))
        subtree = node[index]
        if want_internal && length(subtree) > 1 
            return subtree, index 
        elseif !want_internal && length(subtree) == 1
            return subtree, index 
        end
    end
end

function standard_crossover(root1::Node, root2::Node, rng::AbstractRNG=Random.GLOBAL_RNG)
    _, insertion_point = draw_subtree(root1, rng)
    donation, _ = draw_subtree(root2, rng)
    insert(root1, donation, insertion_point)
end

function aligned_crossover(root1::Node, root2::Node, inputs, targets, rng::AbstractRNG=Random.GLOBAL_RNG)
    _, insertion_point = draw_subtree(root1, rng)
    donation, _ = draw_subtree(root2, rng)
    insert_with_alignment(root1, donation, insertion_point, inputs, targets) |> first
    #insert_with_alignment_nonrecursive(root1, donation, insertion_point, inputs, targets) 
end

function sizefair_mutation(root::Node, symboltable::SymbolTable, rng=Random.GLOBAL_RNG)
    subtree, insertion_point = draw_subtree(root, rng)
    donation = init(symboltable, length(subtree), rng)
    insert(root, donation, insertion_point)
end

function aligned_mutation(root::Node, symboltable::SymbolTable, inputs, targets, rng=Random.GLOBAL_RNG)
    subtree, insertion_point = draw_subtree(root, rng)
    donation = init(symboltable, length(subtree), rng)
    insert_with_alignment(root, donation, insertion_point, inputs, targets) |> first
    #insert_with_alignment_nonrecursive(root, donation, insertion_point, inputs, targets) 
end
