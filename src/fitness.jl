
include("find_dominated.jl")

function evaluate_to_tree(node::Node, setup::ProblemSetup)
    inputs = setup.inputs
    targets = setup.interval_targets

    output = evaluate(node, inputs)
    if all(isfinite, output)
        hits = BitVector([output[i] ∈ targets[i] for i in 1:length(targets)])
        slope, intercept, mse = linear_scale(output, setup.noisy_targets)
        scaled_out = slope .* output .+ intercept
        scaled_hitcount = sum(scaled_out[i] ∈ targets[i] for i in eachindex(targets))
        return Tree(node, hits, slope, intercept, mse, scaled_hitcount)
    else
        return Tree(node, BitVector([false for _ in targets]))
    end
end

function initstrata(setup::ProblemSetup)
    strata = [Tree[] for _ in 1:setup.max_complexity]

    total = 0
    while total < setup.population_size 
        targetsize = mod1(total, setup.max_complexity)
        node = init(setup.symboltable, targetsize, setup.rng)
        tree = evaluate_to_tree(node, setup)
        if isnothing(tree)
            continue
        end
        cmp = complexity(tree)
        if cmp <= length(strata)
            push!(strata[cmp], tree)
            total += 1
        end
    end
    return strata
end

function lexicase(pop::Vector{Tree}, max_lexicase::Int, rng)
    remaining = collect(1:length(pop))
    cases = randperm(rng, length(first(pop).hits))
    for i in 1:max_lexicase 
        case = cases[i] 
        new_remaining = []
        for j in remaining
            if pop[j].hits[case] 
                push!(new_remaining, j)
            end
        end
        if length(new_remaining) == 1 
            return pop[new_remaining[1]]
        elseif length(new_remaining) == 0 
            break
        end
        remaining = new_remaining
    end
    return pop[rand(rng,remaining)]
end

function find_replacement(pop::Vector{Tree})
   indy1 = rand(1:length(pop))
    indy2 = rand(1:length(pop))
    while indy1 == indy2 
        indy2 = rand(1:length(pop))
    end

    if sum(pop[indy1].hits) < sum(pop[indy2].hits)
        indy2 = indy1
    end
    return indy2
end

function find_least_contributor(pop::Vector{Tree})
    distr = zeros(Int, length(first(pop).hits))
    for indy in pop
        distr .+= indy.hits
    end
    highest_index=0
    highest_common=0
    for i in eachindex(pop)
        hits = pop[i].hits
        mn = minimum(distr[hits], init=typemax(Int))
        if mn > highest_common
            highest_common = mn
            highest_index = i 
        end
    end
    return highest_index
end

function iteratestrata!(strata::Vector{Vector{Tree}}, setup::ProblemSetup; standard_gp=false)
    # pick a random population
    pop1 = rand(strata)
    while isempty(pop1)
        pop1 = rand(strata)
    end
    # select next pop uniformly from a window
    range = first(pop1).complexity .+ (-10:10)
    selection = rand(range)
    while selection > length(strata) || selection < 1 || isempty(strata[selection])
        selection = rand(range)
    end
    pop2 = strata[selection]

    indy1 = lexicase(pop1, setup.max_lexicase_comparisons, setup.rng)
    indy2 = lexicase(pop2, setup.max_lexicase_comparisons, setup.rng)
    child = try
        if rand(setup.rng) < setup.cross_mut_prob 
            if standard_gp 
                standard_crossover(indy1.root, indy2.root, setup.rng)
            else
                aligned_crossover(indy1.root, indy2.root, setup.inputs, setup.interval_targets, setup.rng) 
            end
        else
            if standard_gp
                sizefair_mutation(indy1.root, setup.symboltable, setup.rng)
            else 
                aligned_mutation(indy1.root, setup.symboltable, setup.inputs, setup.interval_targets, setup.rng) 
            end
        end
    catch e 
        if e isa DomainError 
            return
        else
            rethrow(e)
        end
    end
    tree = evaluate_to_tree(child, setup)
    if isnothing(tree)
        return
    end
    
    if complexity(tree) <= length(strata)
        push!(strata[complexity(tree)], tree)

        # remove an individual from random pop
        heavy_pop = argmax([length(stratum) for stratum in strata])
        pop3 = strata[heavy_pop]
        while length(pop3) < 2 pop3 = rand(strata) end
        
        # result = find_dominated_trees(pop3)
        # if isnothing(result)
        #     to_delete = find_replacement(pop3)
        # else
        #     to_delete = result[1]
        # end

        to_delete = find_replacement(pop3)
        #to_delete = find_least_contributor(pop3)
        
        deleteat!(pop3, to_delete)
    end
end

function print_report(strata::Vector{Vector{Tree}})
    nIndies = 0
    nHits = 0
    besthits = 0
    best = nothing
    for stratum in strata
        for indy in stratum
            nIndies += 1
            nHits += sum(indy.hits)
            if sum(indy.hits) > besthits 
                besthits = sum(indy.hits)
                best = indy
            end
        end
    end
    println("best $besthits, avg $(round(nHits/nIndies, digits=2)), mse $(round(best.mse, digits=7)) complexity $(best.complexity) avg-pathlen $(round(best.pathlen_complexity / best.complexity, digits=3))")
    return besthits
end

function coordinate_descent(tree::Tree, setup::ProblemSetup, ntries)
    x = setup.inputs 
    t = setup.interval_targets

    node = tree.root
    besthits = sum(tree.hits)
    println("Original: $besthits")
    for _ in 1:ntries 
        sub = rand(1:length(node))
        new_node, evals = insert_with_alignment(node, node[sub], sub, x,t)
        if all(isfinite, evals)
            hits = sum(evals[j] ∈ t[j] for j in eachindex(t))
            if sum(hits) > besthits
                node = new_node 
                besthits = sum(hits)
                println("Improved: $besthits")
            end
        end
    end
    evaluate_to_tree(node, setup)
end


