
include("find_dominated.jl")

function initstrata(setup::ProblemSetup)
    strata = [Tree[] for _ in 1:setup.params.max_complexity]

    total = 0
    nevals = 0
    while total < setup.params.population_size
        targetsize = mod1(total, setup.params.max_complexity)
        node, _ = valid_init(setup.symboltable, targetsize, setup.inputs, setup.rng)
        tree = evaluate_to_tree(node, setup)
        if isnothing(tree)
            continue
        end
        if complexity(tree) > length(strata) 
            continue 
        end
        push!(strata[complexity(tree)], tree)
        total += 1
        nevals += tree.complexity
    end
    return strata, EffortStats(nevals, 0)
end

function initial_fit(strata::Vector{Vector{Tree}}, setup::ProblemSetup)
    for stratum in strata
        for i in eachindex(stratum)
            node = stratum[i].root
            node, _ = insert_with_alignment(node, node, 1, 1, setup.inputs, setup.interval_targets, false)
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
    return pop[rand(rng, remaining)]
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
    return rand(rng, remaining)
end

# Dormant (unused) selector: still on the raw maximize framing (primary-hit count),
# not two-band loss. Update to `.loss` if it is ever wired back into iteratestrata!.
function tourney(pop::Vector{Tree}, t=5)
    champion = rand(pop)
    if t > 1
        challenger = tourney(pop, t - 1)
        if sum(challenger.hits) > sum(champion.hits)
            champion = challenger
        end
    end
    return champion
end


function find_replacement(pop::Vector{Tree})
    @assert length(pop) > 1
    indy1 = rand(1:length(pop))
    indy2 = rand(1:length(pop))
    while indy1 == indy2
        indy2 = rand(1:length(pop))
    end

    # Route replacement through two-band loss (lower is better): delete the worse.
    if pop[indy1].loss > pop[indy2].loss
        indy2 = indy1
    end
    return indy2
end

function find_least_contributor(pop::Vector{Tree})
    distr = zeros(Int, length(first(pop).hits))
    for indy in pop
        distr .+= indy.hits
    end
    highest_index = 0
    highest_common = 0
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

# Parent selector dispatch: residual ε-lexicase when enabled, else two-band lexicase.
# Only PARENT selection is affected; constant fitting, replacement, and hit/loss
# bookkeeping are untouched.
function select_parent(pop::Vector{Tree}, setup::ProblemSetup)
    if setup.params.use_residual_lexicase
        return residual_eps_lexicase(pop, setup.noisy_targets,
                                     setup.params.max_lexicase_comparisons, setup.rng)
    else
        return two_band_lexicase(pop, setup.params.max_lexicase_comparisons, setup.rng)
    end
end

function iteratestrata!(strata::Vector{Vector{Tree}}, setup::ProblemSetup, effort::EffortStats) :: Int
    method = setup.params.method

    if setup.params.use_tournament_stratum
        # Dormant path (off by default): stratum ranking still uses the raw maximize
        # framing (max primary-hit count), not two-band loss.
        maxhits = 0
        maxhitstratum = 1
        for i in eachindex(strata)
            if isempty(strata[i])
                continue
            end
            hits = maximum(sum(indy.hits) for indy in strata[i])
            if hits > maxhits
                maxhits = hits 
                maxhitstratum = i 
            end
        end
        i1 = rand(eachindex(strata))
        while isempty(strata[i1])
            i1 = rand(eachindex(strata))
        end
        i2 = rand(eachindex(strata))
        while isempty(strata[i2])
            i2 = rand(eachindex(strata))
        end
        # pick lowest
        if abs(i1 - maxhitstratum) > abs(i2 - maxhitstratum) 
            i1 = i2
        end
        pop1 = strata[i1]
    else
        # pick a random population
        pop1 = rand(strata)
        while isempty(pop1)
            pop1 = rand(strata)
        end
    end

    indy1 = select_parent(pop1, setup)
    # indy1 = tourney(pop1)
    # indy2 = tourney(pop2)
    if rand(setup.rng) < setup.params.cross_mut_prob
        # select next pop uniformly from a window
        range = first(pop1).complexity .+ (-10:10)
        selection = rand(range)
        while selection > length(strata) || selection < 1 || isempty(strata[selection])
            selection = rand(range)
        end
        pop2 = strata[selection]
        indy2 = select_parent(pop2, setup)
        if method == Standard || method == ConstantStab
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
        elseif method == ConstantStab
            if rand(setup.rng) < setup.params.constant_stab_probability
                child = constant_stab_mutation(indy1, setup.inputs, setup.interval_targets, effort, setup.rng)
            else
                child = sizefair_mutation(indy1.root, setup.symboltable, effort, setup.rng)
            end
        else
            child = recursive_aligned_mutation(indy1.root, setup.symboltable, setup.inputs, setup.interval_targets, effort, setup.rng)
        end
    end
    tree = evaluate_to_tree(child, setup)
    
    #invalid?
    if isnothing(tree) || sum(tree.hits) == 0 || complexity(tree) > length(strata)
        return 0
    end

    push!(strata[complexity(tree)], tree)

    # remove an individual from random pop
    heavy_pop = argmax([length(stratum) for stratum in strata])
    pop3 = strata[heavy_pop]
    while length(pop3) < 2
        pop3 = rand(strata)
    end

    # result = find_dominated_trees(pop3)
    # if isnothing(result)
    #     to_delete = find_replacement(pop3)
    # else
    #     to_delete = result[1]
    # end

    to_delete = find_replacement(pop3)
    #to_delete = inverse_lexicase(pop3, setup.max_lexicase_comparisons, setup.rng)
    #to_delete = find_least_contributor(pop3)

    deleteat!(pop3, to_delete)
    return sum(tree.hits)
end

function coordinate_descent(tree::Tree, setup::ProblemSetup, ntries)
    x = setup.inputs
    t = setup.interval_targets

    nmse = 0
    for i in 1:ntries
        sub = rand(1:length(tree))
        new_node, _ = insert_with_alignment(tree.root, tree.root[sub], sub, 1, x, t, false)
        new_tree = evaluate_to_tree(new_node, setup)
        if complexity(new_tree) < complexity(tree)
            tree = new_tree
            println("1) Improved $i: hits=$(sum(tree.hits)) complexity=$(tree.complexity), mse=$(tree.mse)")
        end

        if sum(new_tree.hits) < sum(tree.hits)
            continue 
        end

        if sum(new_tree.hits) > sum(tree.hits)
            tree = new_tree
            println("2) Improved $i: hits=$(sum(tree.hits)) complexity=$(tree.complexity), mse=$(tree.mse)")
        end

        if new_tree.mse <= tree.mse 
            tree = new_tree
            nmse += 1
            if nmse > 100
                println("3) Improved $i: hits=$(sum(tree.hits)) complexity=$(tree.complexity), mse=$(tree.mse)")
                nmse = 0
            end
        end

    end
    tree
end



