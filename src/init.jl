using Random


function init(symboltable::SymbolTable, targetsize::Int, rng::AbstractRNG=Random.GLOBAL_RNG)
    # Initialize a random tree with the given target size
    if targetsize == 1
        if rand(rng, 1:symboltable.nvars+1) == 1
            return Constant(randn(rng) / randn(rng)) # Cauchy distributrion for ERC
        else
            return Var(rand(rng, 1:symboltable.nvars))
        end
    else
        if targetsize == 2 && length(symboltable.unaries) > 0
            # Only unary node can have size 2
            fun = rand(rng, symboltable.unaries)
            child = init(symboltable, 1, rng)
            return UnaryNode(fun, child)
        else
            # Binary node
            selection = rand(rng, 1:length(symboltable.binaries) + length(symboltable.unaries))
            if selection <= length(symboltable.binaries)
                # stopgap for when we have no unaries and size 2 is the target
                targetsize = max(3, targetsize)
                # Select binary node
                fun = symboltable.binaries[selection]
                leftsize = rand(rng, 1:targetsize-2)
                rightsize = targetsize - 1 - leftsize
                leftchild = init(symboltable, leftsize, rng)
                rightchild = init(symboltable, rightsize, rng)
                return BinaryNode(fun, leftchild, rightchild)
            else # Select unary node
                fun = symboltable.unaries[selection - length(symboltable.binaries)]
                child = init(symboltable, targetsize - 1, rng)
                return UnaryNode(fun, child)
            end
        end
    end
end

# Initialize to a tree that is at most targetsize, but for which the outputs are finite
function valid_init(symboltable::SymbolTable, targetsize::Int, inputs, rng::AbstractRNG=Random.GLOBAL_RNG) :: Tuple{Node, Vector{Float64}}
    # Initialize a random tree with the given target size
    if targetsize == 1
        if rand(rng, 1:symboltable.nvars+1) == 1
            node = Constant(randn(rng) / randn(rng)) # Cauchy distributrion for ERC
            return node, evaluate(node, inputs)
        else
            node = Var(rand(rng, 1:symboltable.nvars))
            return node, evaluate(node, inputs)
        end
    else
        if targetsize == 2 && length(symboltable.unaries) > 0
            # Only unary node can have size 2
            fun = rand(rng, symboltable.unaries)
            child, childevals = valid_init(symboltable, targetsize - 1, inputs, rng)
            evals = evaluate(fun, childevals)
            if !all(isfinite, evals)
                return child, childevals
            end

            return UnaryNode(fun, child), evals
        else
            # Binary node
            selection = rand(rng, 1:length(symboltable.binaries) + length(symboltable.unaries))
            if selection <= length(symboltable.binaries)
                # stopgap for when we have no unaries and size 2 is the target
                targetsize = max(3, targetsize)
                # Select binary node
                fun = symboltable.binaries[selection]
                leftsize = rand(rng, 1:targetsize-2)
                rightsize = targetsize - 1 - leftsize
                leftchild, leftevals = valid_init(symboltable, leftsize, inputs, rng)
                rightchild, rightevals = valid_init(symboltable, rightsize, inputs, rng)
                evals = evaluate(fun, leftevals, rightevals)
                if !all(isfinite, evals)
                    return rand(rng, Bool) ? (leftchild, leftevals) : (rightchild, rightevals)
                end
                return BinaryNode(fun, leftchild, rightchild), evals
            else # Select unary node
                fun = symboltable.unaries[selection - length(symboltable.binaries)]
                child, childevals = valid_init(symboltable, targetsize - 1, inputs, rng)
                evals = evaluate(fun, childevals)
                if !all(isfinite, evals)
                    return child, childevals
                end
                return UnaryNode(fun, child), evals
            end
        end
    end
end

