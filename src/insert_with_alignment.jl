"""
insert with alignment: insert a node into a tree at index i, and update the alignment of the recipient node to include the allowed intervals of the donation node.
It will also update the added values of all nodes upward
"""
function insert_with_alignment(recipient::BinaryNode, donation::Node, i::Int, inputs, targets::Vector{Interval{Float64,Closed,Closed}})

    if i == 1
        return aligned_node(donation, targets, inputs) # find best additive constant
    elseif i <= 1 + length(recipient.left)
        # Evaluate right child -- keep
        right_evals = evaluate(recipient.right, inputs)
        # compute surrogate targets for left child
        new_targets = _surrogate_left(recipient.fun, targets .- recipient.addition.value, right_evals)
        # insert donation into left child, receive new child and left child evaluations
        new_left, left_evals = insert_with_alignment(recipient.left, donation, i - 1, inputs, new_targets)
        # Re-evaluate the new node without its added value
        new_evals = evaluate(recipient.fun, left_evals, right_evals)

        # compute an added value for the new node that aligns with the target (new_evals gets updated with that constant value)
        added_value, updated_evals = compute_added_value(new_evals, targets)
        # Return the re-optimized node with the new left child as well as the evaluations
        return BinaryNode(recipient.fun, new_left, recipient.right, added_value), updated_evals
    else
        # Evaluate left child -- keep
        left_evals = evaluate(recipient.left, inputs)
        # compute surrogate targets for right child
        new_targets = _surrogate_right(recipient.fun, targets .- recipient.addition.value, left_evals)
        # insert donation into right child, receive new child and right child evaluations
        new_right, right_evals = insert_with_alignment(recipient.right, donation, i - 1 - length(recipient.left), inputs, new_targets)
        # Re-evaluate the new node without its added value
        new_evals = evaluate(recipient.fun, left_evals, right_evals)

        # compute an added value for the new node that aligns with the target (new_evals gets updated with that constant value)
        added_value, updated_evals = compute_added_value(new_evals, targets)
        # Return the re-optimized node with the new right child as well as the evaluations
        return BinaryNode(recipient.fun, recipient.left, new_right, added_value), updated_evals
    end
end

function insert_with_alignment(recipient::UnaryNode, donation::Node, i::Int, inputs, targets::Vector{Interval{Float64,Closed,Closed}})
    
    if i == 1
        return aligned_node(donation, targets, inputs)
    else
        new_targets = inverse(recipient.fun, targets .- recipient.addition.value)
        new_child, child_evals = insert_with_alignment(recipient.child, donation, i - 1, inputs, new_targets)
        new_evals = evaluate(recipient.fun, child_evals)

        added_value, updated_evals = compute_added_value(new_evals, targets)
        return UnaryNode(recipient.fun, new_child, added_value), updated_evals
    end
end

function insert_with_alignment(::Union{Var, Constant}, donation::Node, i::Int, inputs, targets::Vector{Interval{Float64,Closed,Closed}})
    if i == 1
        return aligned_node(donation, targets, inputs)
    else
        error("Index out of bounds for node")
    end
end

function compute_added_value(evals, targets::Vector{Interval{Float64,Closed,Closed}}, rng::AbstractRNG=Random.GLOBAL_RNG)
    valid = findall(isfinite, evals)
    if isempty(valid)
        return AddedValue(IntervalSet(Interval{Float64,Closed,Closed}[]), 0.0), evals
    end
    res, _ = max_overlap_region(targets[valid] .- evals[valid])
    value = select_constant(res, rng)
    return AddedValue(res, value), evals .+ value
end

# Map non-finite sibling evals to the invalid sentinel so they don't corrupt the
# interval inverse computation (Interval(NaN,NaN) is not a valid Intervals.jl type).
_surrogate_left(fun, targets, se)  = [isfinite(se[j]) ? leftinverse(fun, targets[j], se[j])  : Interval(typemax(Float64), typemax(Float64)) for j in eachindex(targets)]
_surrogate_right(fun, targets, se) = [isfinite(se[j]) ? rightinverse(fun, targets[j], se[j]) : Interval(typemax(Float64), typemax(Float64)) for j in eachindex(targets)]

function aligned_node(node::Node, targets::Vector{Interval{Float64,Closed,Closed}}, inputs)

    # Evaluate node without its added value
    out = evaluate(node, inputs) .- node.addition.value
    # Find an added value that optimizes hit count
    added_value, new_outputs = compute_added_value(out, targets)

    if node isa Var
        node = Var(node.index, added_value)
    elseif node isa Constant
        node = Constant(added_value)
    elseif node isa BinaryNode
        node = BinaryNode(node.fun, node.left, node.right, added_value)
    elseif node isa UnaryNode
        node = UnaryNode(node.fun, node.child, added_value)
    else
        error("Unknown node type")
    end
    # Return the re-optimized node with the new added value as well as the evaluations
    return node, new_outputs
end

function insert_with_alignment_nonrecursive(recipient::BinaryNode, donation::Node, i::Int, inputs, targets::Vector{Interval{Float64,Closed,Closed}})

    if i == 1
        return aligned_node(donation, targets, inputs) |> first # find best additive constant
    elseif i <= 1 + length(recipient.left)
        # Evaluate right child -- keep
        right_evals = evaluate(recipient.right, inputs)
        # compute surrogate targets for left child
        new_targets = _surrogate_left(recipient.fun, targets .- recipient.addition.value, right_evals)
        # insert donation into left child, receive new child and left child evaluations
        new_left = insert_with_alignment_nonrecursive(recipient.left, donation, i - 1, inputs, new_targets)

        # Return the re-optimized node with the new left child as well as the evaluations
        return BinaryNode(recipient.fun, new_left, recipient.right, recipient.addition)
    else
        # Evaluate left child -- keep
        left_evals = evaluate(recipient.left, inputs)
        # compute surrogate targets for right child
        new_targets = _surrogate_right(recipient.fun, targets .- recipient.addition.value, left_evals)
        # insert donation into right child, receive new child and right child evaluations
        new_right = insert_with_alignment_nonrecursive(recipient.right, donation, i - 1 - length(recipient.left), inputs, new_targets)
        return BinaryNode(recipient.fun, recipient.left, new_right, recipient.addition)
    end
end

function insert_with_alignment_nonrecursive(recipient::UnaryNode, donation::Node, i::Int, inputs, targets::Vector{Interval{Float64,Closed,Closed}})
    
    if i == 1
        return aligned_node(donation, targets, inputs) |> first
    else
        new_targets = inverse(recipient.fun, targets .- recipient.addition.value)
        new_child = insert_with_alignment_nonrecursive(recipient.child, donation, i - 1, inputs, new_targets)
        return UnaryNode(recipient.fun, new_child, recipient.addition)
    end
end

function insert_with_alignment_nonrecursive(::Union{Var, Constant}, donation::Node, i::Int, inputs, targets::Vector{Interval{Float64,Closed,Closed}})
    if i == 1
        return aligned_node(donation, targets, inputs) |> first
    else
        error("Index out of bounds for node")
    end
end
