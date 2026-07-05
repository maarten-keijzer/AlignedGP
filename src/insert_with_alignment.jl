"""
insert with alignment: insert a node into a tree at index i, and update the alignment of the recipient node to include the allowed intervals of the donation node.
It will also update the added values of all nodes upward
"""
function insert_with_alignment(recipient::BinaryNode, donation::Node, i::Int, d::Int, inputs, targets::Vector{CIntervals}, recursive_stabbing)

    if i == 1
        return align_node(donation, targets, inputs)..., d # find best additive constant
    elseif i <= 1 + length(recipient.left)
        # Evaluate right child -- keep
        right_evals = evaluate(recipient.right, inputs)
        # compute surrogate targets for left child
        new_targets = _surrogate_left(recipient.fun, targets .- recipient.addition.value, right_evals)
        # insert donation into left child, receive new child and left child evaluations
        new_left, left_evals, d = insert_with_alignment(recipient.left, donation, i - 1, d + 1, inputs, new_targets, recursive_stabbing)
        # Re-evaluate the new node without its added value
        new_evals = evaluate(recipient.fun, left_evals, right_evals)

        # compute an added value for the new node that aligns with the target (new_evals gets updated with that constant value)
        if recursive_stabbing
            added_value, updated_evals = compute_added_value(new_evals, targets)
        else
            added_value = recipient.addition
            updated_evals = new_evals .+ recipient.addition.value
        end

        # Return the re-optimized node with the new left child as well as the evaluations
        return BinaryNode(recipient.fun, new_left, recipient.right, added_value), updated_evals, d
    else
        # Evaluate left child -- keep
        left_evals = evaluate(recipient.left, inputs)
        # compute surrogate targets for right child
        new_targets = _surrogate_right(recipient.fun, targets .- recipient.addition.value, left_evals)
        # insert donation into right child, receive new child and right child evaluations
        new_right, right_evals, d = insert_with_alignment(recipient.right, donation, i - 1 - length(recipient.left), d + 1, inputs, new_targets, recursive_stabbing)
        # Re-evaluate the new node without its added value
        new_evals = evaluate(recipient.fun, left_evals, right_evals)

        # compute an added value for the new node that aligns with the target (new_evals gets updated with that constant value)
        if recursive_stabbing
            added_value, updated_evals = compute_added_value(new_evals, targets)
        else
            added_value = recipient.addition
            updated_evals = new_evals .+ recipient.addition.value
        end
        # Return the re-optimized node with the new right child as well as the evaluations
        return BinaryNode(recipient.fun, recipient.left, new_right, added_value), updated_evals, d
    end
end

function insert_with_alignment(recipient::UnaryNode, donation::Node, i::Int, d::Int, inputs, targets::Vector{CIntervals}, recursive_stabbing)

    if i == 1
        return align_node(donation, targets, inputs)..., d
    else
        new_targets = inverse(recipient.fun, targets .- recipient.addition.value)
        new_child, child_evals, d = insert_with_alignment(recipient.child, donation, i - 1, d+1, inputs, new_targets, recursive_stabbing)
        new_evals = evaluate(recipient.fun, child_evals)

        if recursive_stabbing
            added_value, updated_evals = compute_added_value(new_evals, targets)
        else
            added_value = recipient.addition
            updated_evals = new_evals .+ recipient.addition.value
        end
        return UnaryNode(recipient.fun, new_child, added_value), updated_evals, d
    end
end

function insert_with_alignment(::Union{Var, Constant}, donation::Node, i::Int, d::Int, inputs, targets::Vector{CIntervals}, _)
    return align_node(donation, targets, inputs)..., d
end

_hits(x, t) = sum(x[i] ∈ t[i] for i in eachindex(x))

function compute_added_value(evals, targets::Vector{CIntervals}, rng::AbstractRNG=Random.GLOBAL_RNG)
    # Use exact (non-narrowing) shift so that point targets [x,x] with eval==x
    # produce [0,0] (not invalid_interval), preserving exact hits in the depth count.
    res, depth = max_overlap_region(targets .- evals)
    value = select_constant(res, rng)
    added_value = AddedValue(res, value, depth)
    evals = evals .+ value

    @assert _hits(evals, targets) >= depth "$(_hits(evals, targets)), $depth $(maximum(evals)) $(minimum(evals)) $value"

    return added_value, evals
end

# Map non-finite sibling evals to the invalid sentinel so they don't corrupt the
# interval inverse computation (NaN is not a valid CInterval).
_surrogate_left(fun, targets, se)  = CIntervals[isfinite(se[j]) ? leftinverse(fun, targets[j], se[j])  : CIntervals() for j in eachindex(targets)]
_surrogate_right(fun, targets, se) = CIntervals[isfinite(se[j]) ? rightinverse(fun, targets[j], se[j]) : CIntervals() for j in eachindex(targets)]

function _optimize(node::BinaryNode, targets::Vector{CIntervals}, inputs, rng::AbstractRNG)
    right_evals_pre = evaluate(node.right, inputs)
    left_targets = _surrogate_left(node.fun, targets .- node.addition.value, right_evals_pre)

    new_left, left_evals = _optimize(node.left, left_targets, inputs, rng)

    right_targets = _surrogate_right(node.fun, targets .- node.addition.value, left_evals)
    new_right, right_evals_post = _optimize(node.right, right_targets, inputs, rng)
    new_evals = evaluate(node.fun, left_evals, right_evals_post)
    added_value, updated_evals = compute_added_value(new_evals, targets, rng)

    left_hits   = new_left.addition.hits
    right_hits  = new_right.addition.hits

    parent_hits = added_value.hits
    if parent_hits < left_hits || parent_hits < right_hits
        @show parent_hits, left_hits, right_hits
    end
    @assert parent_hits >= right_hits  "[DEBUG-a4f2] T3 < T2: parent=$parent_hits right=$right_hits left=$left_hits, fun=$(node.fun)"
    @assert parent_hits == added_value.hits "$parent_hits != $(added_value.hits)"

    return BinaryNode(node.fun, new_left, new_right, added_value), updated_evals
end

function _optimize(node::UnaryNode, targets::Vector{CIntervals}, inputs, rng::AbstractRNG)
    child_targets = inverse(node.fun, targets .- node.addition.value)
    new_child, child_evals = _optimize(node.child, child_targets, inputs, rng)
    new_evals = evaluate(node.fun, child_evals)
    added_value, updated_evals = compute_added_value(new_evals, targets, rng)

    parent_hits = _hits(updated_evals, targets)

    @assert parent_hits >= new_child.addition.hits

    return UnaryNode(node.fun, new_child, added_value), updated_evals
end

function _optimize(node::Var, targets::Vector{CIntervals}, inputs, rng::AbstractRNG)
    out = evaluate(node, inputs) .- node.addition.value
    added_value, evals = compute_added_value(out, targets, rng)
    return Var(node.index, added_value), evals
end

function _optimize(node::Constant, targets::Vector{CIntervals}, inputs, rng::AbstractRNG)
    out = evaluate(node, inputs) .- node.addition.value
    added_value, evals = compute_added_value(out, targets, rng)
    return Constant(added_value), evals
end

"""
    optimize(root, inputs, targets, rng) -> (Node, evals)

Recursively set the best additive constant on every node in the tree using
interval inversion. Targets are pushed down via left/right-inverse and
`inverse`, evaluations are returned back up. Assumes all subtree evaluations
are finite.
"""
function optimize(root::Node, inputs, targets::Vector{CIntervals}, rng::AbstractRNG=Random.GLOBAL_RNG)
    _optimize(root, targets, inputs, rng)
end

function align_node(node::Node, targets::Vector{CIntervals}, inputs)

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

