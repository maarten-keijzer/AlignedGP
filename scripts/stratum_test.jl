using AlignedGP
using Random

mutable struct Stratum
    pop::Vector{Tree}
    maxhits::Int 
end

function Base.push!(stratum::Stratum, tree::Tree) 
    stratum.maxhits = max(stratum.maxhits, sum(tree.hits))
    push!(stratum.pop, tree)
end

mutable struct SizeStrata
    strata::Vector{Stratum}
    nindividuals::Int
    maxcomplexity::Int
    beststratumlocation::Int
end

Base.length(strata::SizeStrata) = strata.nindividuals

function Base.iterate(s::SizeStrata, state=(1, 1))
    si, ti = state
    while si <= length(s.strata)
        stratum = s.strata[si]
        if ti <= length(stratum.pop)
            return (stratum.pop[ti], (si, ti + 1))
        end
        si += 1
        ti = 1
    end
    return nothing
end

function Base.push!(strata::SizeStrata, tree::Tree)
    cmpl = complexity(tree)
    if cmpl <= strata.maxcomplexity
        strata.nindividuals += 1
        while length(strata.strata) < cmpl
            push!(strata.strata, Stratum([], 0))
        end
        push!(strata.strata[cmpl], tree)
        best = strata.beststratumlocation
        new_hits = strata.strata[cmpl].maxhits
        if best == -1 || new_hits > strata.strata[best].maxhits ||
                (new_hits == strata.strata[best].maxhits && cmpl < best)
            strata.beststratumlocation = cmpl
        end
        return
    end
end


# create a distribution tailing of from the first bin
function cull!(strata, rng = Random.GLOBAL_RNG)
    best = strata.beststratumlocation
    distr = zeros(length(strata.strata))
    for i in 1:best 
        distr[i] = 1
    end
    for i in best+1:length(strata.strata)
        distr[i] = distr[i-1] * 0.99
    end
    target = distr .* strata.nindividuals ./ sum(distr)

    todel = argmax(
        length(strata.strata[i].pop) <= 1 ? -Inf : length(strata.strata[i].pop) - target[i]
        for i in eachindex(strata.strata)
    )
    
    stratum = strata.strata[todel]
    pop = stratum.pop 

    @assert length(pop) > 1

    # find a deletion point
    t1 = rand(rng, eachindex(pop))
    t2 = t1
    while t2 == t1 
        t2 = rand(rng, eachindex(pop))
    end
    if pop[t2].hits > pop[t1].hits
        t2 = t1
    end
    deleteat!(pop, t2)
    strata.nindividuals -= 1
end

setup = keijzer1()
strata = SizeStrata([], 0, 150, -1)
for i = 1:10_000 

    node, _ = valid_init(setup.symboltable, rand(1:150), setup.inputs)
    tree = AlignedGP.evaluate_to_tree(node, setup)
    push!(strata, tree)
    if length(strata) > 5000 
        cull!(strata, setup.rng)
    end
end 
length(strata)
using Plots
bar(length.(stratum.pop for stratum in strata.strata))