
include("find_dominated.jl")

function initstrata(setup::ProblemSetup)
    strata = [Tree[] for _ in 1:setup.params.max_complexity]

    total = 0
    while total < setup.params.population_size
        targetsize = mod1(total, setup.params.max_complexity)
        node,_ = valid_init(setup.symboltable, targetsize, setup.inputs, setup.rng)
        tree = evaluate_to_tree(node, setup)
        if isnothing(tree)
            continue
        end
        cmp = min(length(strata), complexity(tree))
        push!(strata[cmp], tree)
        total += 1
    end
    return strata, EffortStats(0, 0)
end

function initial_fit(strata::Vector{Vector{Tree}}, setup::ProblemSetup)
    for stratum in strata 
        for i in eachindex(stratum)
            node = stratum[i].root;
            node,_ = insert_with_alignment(node, node, 1, 1, setup.inputs, setup.interval_targets, false)
            stratum[i] = evaluate_to_tree(node, setup)      
        end
    end
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


function inverse_lexicase(pop::Vector{Tree}, max_lexicase::Int, rng)
    remaining = collect(1:length(pop))
    cases = randperm(rng, length(first(pop).hits))
    for i in 1:max_lexicase 
        case = cases[i] 
        new_remaining = []
        for j in remaining
            if !pop[j].hits[case] 
                push!(new_remaining, j)
            end
        end
        if length(new_remaining) == 1 
            return new_remaining[1]
        elseif length(new_remaining) == 0 
            break
        end
        remaining = new_remaining
    end
    return rand(rng,remaining)
end

function tourney(pop::Vector{Tree}, rng::AbstractRNG, t=5)
    champion = rand(rng, pop)
    if t > 1
        challenger = tourney(pop, rng, t-1)
        if sum(challenger.hits) > sum(champion.hits)
            champion = challenger
        end
    end
    return champion
end


function find_replacement(pop::Vector{Tree}, rng::AbstractRNG)
    indy1 = rand(rng, 1:length(pop))
    indy2 = rand(rng, 1:length(pop))
    while indy1 == indy2
        indy2 = rand(rng, 1:length(pop))
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

function iteratestrata!(strata::Vector{Vector{Tree}}, setup::ProblemSetup, effort::EffortStats; method::OptMethod=RecursiveStab)
    # pick a random population
    pop1 = rand(setup.rng, strata)
    while isempty(pop1)
        pop1 = rand(setup.rng, strata)
    end

    indy1 = lexicase(pop1, setup.params.max_lexicase_comparisons, setup.rng)
    # indy1 = tourney(pop1, setup.rng)
    # indy2 = tourney(pop2, setup.rng)
        if rand(setup.rng) < setup.params.cross_mut_prob
            # select next pop uniformly from a window
            range = first(pop1).complexity .+ (-10:10)
            selection = rand(setup.rng, range)
            while selection > length(strata) || selection < 1 || isempty(strata[selection])
                selection = rand(setup.rng, range)
            end
            pop2 = strata[selection]
            indy2 = lexicase(pop2, setup.params.max_lexicase_comparisons, setup.rng)
            if method == Standard 
                child = standard_crossover(indy1.root, indy2.root, effort, setup.rng)
            elseif method == Stab
                child = aligned_crossover(indy1.root, indy2.root, setup.inputs, setup.interval_targets, effort, setup.rng)
            else
                child = recursive_aligned_crossover(indy1.root, indy2.root, setup.inputs, setup.interval_targets, effort, setup.rng)
            end
        else
            if method == Standard
                child = sizefair_mutation(indy1.root, setup.symboltable, effort, setup.rng)
            elseif method == Stab
                child = aligned_mutation(indy1.root, setup.symboltable, setup.inputs, setup.interval_targets, effort, setup.rng) 
            else 
                child = recursive_aligned_mutation(indy1.root, setup.symboltable, setup.inputs, setup.interval_targets, effort, setup.rng) 
            end
        end
    tree = evaluate_to_tree(child, setup)
    if isnothing(tree) || sum(tree.hits) == 0
        return 0
    end
    
    if complexity(tree) <= length(strata)
        push!(strata[complexity(tree)], tree)

        # remove an individual from random pop
        heavy_pop = argmax([length(stratum) for stratum in strata])
        pop3 = strata[heavy_pop]
        while length(pop3) < 2 pop3 = rand(setup.rng, strata) end
        
        # result = find_dominated_trees(pop3)
        # if isnothing(result)
        #     to_delete = find_replacement(pop3)
        # else
        #     to_delete = result[1]
        # end

        to_delete = find_replacement(pop3, setup.rng)
        #to_delete = inverse_lexicase(pop3, setup.max_lexicase_comparisons, setup.rng)
        #to_delete = find_least_contributor(pop3)
        
        deleteat!(pop3, to_delete)
    end
    sum(tree.hits)
end

function print_report(strata::Vector{Vector{Tree}}, effort)
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
    println("best $besthits, avg $(round(nHits/nIndies, digits=2)), mse $(round(best.mse, digits=7)) complexity $(best.complexity) avg-pathlen $(round(best.pathlen_complexity / best.complexity, digits=3)), effort $(round(log10(effort), digits=4))")
    return besthits
end

function coordinate_descent(tree::Tree, setup::ProblemSetup, ntries)
    x = setup.inputs 
    t = setup.interval_targets

    node = tree.root
    besthits = sum(tree.hits)
    println("Original: $besthits")
    for _ in 1:ntries
        sub = rand(setup.rng, 1:length(node))
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


