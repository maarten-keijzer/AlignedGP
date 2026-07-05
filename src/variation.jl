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

function standard_crossover(root1::Node, root2::Node, effort::EffortStats, rng::AbstractRNG=Random.GLOBAL_RNG)
    _, insertion_point = draw_subtree(root1, rng)
    donation, _ = draw_subtree(root2, rng)
    tree, d = insert(root1, donation, insertion_point, 1)
    
    effort.sum_evals += d - 1 
    
    tree
end

function sizefair_mutation(root::Node, symboltable::SymbolTable, effort::EffortStats, rng=Random.GLOBAL_RNG)
    subtree, insertion_point = draw_subtree(root, rng)
    donation = init(symboltable, length(subtree), rng)
    tree, d = insert(root, donation, insertion_point, 1)
    
    effort.sum_evals += length(donation)
    effort.sum_evals += d - 1 
    
    tree
end

function aligned_crossover(root1::Node, root2::Node, inputs, targets, effort::EffortStats, rng::AbstractRNG=Random.GLOBAL_RNG)
    _, insertion_point = draw_subtree(root1, rng)
    donation, _ = draw_subtree(root2, rng)
    tree, _, d = insert_with_alignment(root1, donation, insertion_point, 1, inputs, targets, false)
    
    effort.sum_stabs += 1
    effort.sum_evals += 3d 
    
    tree 
end

function aligned_mutation(root::Node, symboltable::SymbolTable, inputs, targets, effort::EffortStats, rng=Random.GLOBAL_RNG)
    subtree, insertion_point = draw_subtree(root, rng)
    donation = init(symboltable, length(subtree), rng)
    tree, _, d = insert_with_alignment(root, donation, insertion_point, 1, inputs, targets, false)
    
    effort.sum_evals += length(donation)
    effort.sum_stabs += 1
    effort.sum_evals += 3d 
    
    tree 
end

function recursive_aligned_crossover(root1::Node, root2::Node, inputs, targets, effort::EffortStats, rng::AbstractRNG=Random.GLOBAL_RNG)
    _, insertion_point = draw_subtree(root1, rng)
    donation, _ = draw_subtree(root2, rng)
    tree, _, d = insert_with_alignment(root1, donation, insertion_point, 1, inputs, targets, true)
    
    effort.sum_stabs += d
    
    tree 
end

function recursive_aligned_mutation(root::Node, symboltable::SymbolTable, inputs, targets, effort::EffortStats, rng=Random.GLOBAL_RNG)
    subtree, insertion_point = draw_subtree(root, rng)
    donation = init(symboltable, length(subtree), rng)
    tree, _, d = insert_with_alignment(root, donation, insertion_point, 1, inputs, targets, true)
    
    effort.sum_evals += length(donation)
    effort.sum_stabs += d
    
    tree 
end
