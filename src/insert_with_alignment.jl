"""
insert with alignment: insert a node into a tree at index i, and update the alignment of the recipient node to include the allowed intervals of the donation node.
It will also update the added values of all nodes upward
"""
function insert_with_alignment(recipient::BinaryNode, donation::Node, i::Int, d::Int, inputs, targets::IntervalVector, recursive_stabbing; circular::Bool = false)

    if nintervals(targets) > 10length(targets)
        #@warn "Many intervals generated: $(nintervals(targets))"
    end

    if i == 1 
        return align_node(donation, targets, inputs; circular=circular)..., d # find best additive constant
    elseif i <= 1 + length(recipient.left)
        # Evaluate right child -- keep
        right_evals = evaluate(recipient.right, inputs)

        # compute surrogate targets for left child 
        new_targets = inverse(recipient.fun, targets - recipient.addition.value, right_evals)

        # insert donation into left child, receive new child and left child evaluations
        new_left, left_evals, d = insert_with_alignment(recipient.left, donation, i - 1, d + 1, inputs, new_targets, recursive_stabbing)
        # Re-evaluate the new node without its added value
        new_evals = evaluate(recipient.fun, left_evals, right_evals)

        # compute an added value for the new node that aligns with the target (new_evals gets updated with that constant value)
        if recursive_stabbing
            added_value, updated_evals = compute_added_value(new_evals, targets; circular=circular)
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
        new_targets = inverse(recipient.fun, left_evals, targets - recipient.addition.value)

        # insert donation into right child, receive new child and right child evaluations
        new_right, right_evals, d = insert_with_alignment(recipient.right, donation, i - 1 - length(recipient.left), d + 1, inputs, new_targets, recursive_stabbing)
        # Re-evaluate the new node without its added value
        new_evals = evaluate(recipient.fun, left_evals, right_evals)

        # compute an added value for the new node that aligns with the target (new_evals gets updated with that constant value)
        if recursive_stabbing
            added_value, updated_evals = compute_added_value(new_evals, targets; circular=circular)
        else
            added_value = recipient.addition
            updated_evals = new_evals .+ recipient.addition.value
        end
        # Return the re-optimized node with the new right child as well as the evaluations
        return BinaryNode(recipient.fun, recipient.left, new_right, added_value), updated_evals, d
    end
end

function insert_with_alignment(recipient::UnaryNode, donation::Node, i::Int, d::Int, inputs, targets::IntervalVector, recursive_stabbing; circular::Bool=false)

    if i == 1
        return align_node(donation, targets, inputs; circular=circular)..., d
    else
        new_targets = inverse(recipient.fun, targets - recipient.addition.value)

        # The child is the direct descendant of this node; if this node is `sin`
        # its constant lives on a circle (mod 2π), so stab it circularly (Situation 1).
        child_circular = isperiodic(recipient.fun)
        new_child, child_evals, d = insert_with_alignment(recipient.child, donation, i - 1, d+1, inputs, new_targets, recursive_stabbing, circular=child_circular)

        new_evals = evaluate(recipient.fun, child_evals)

        if recursive_stabbing
            added_value, updated_evals = compute_added_value(new_evals, targets; circular=circular)
        else
            added_value = recipient.addition
            updated_evals = new_evals .+ recipient.addition.value
        end
        return UnaryNode(recipient.fun, new_child, added_value), updated_evals, d
    end
end

function insert_with_alignment(::Union{Var, Constant}, donation::Node, i::Int, d::Int, inputs, targets::IntervalVector, recursive_stabbing; circular::Bool=false,)
    return align_node(donation, targets, inputs; circular=circular)..., d
end

canonicalize_circular(c::Real, C::Real = 2π) = mod(c + C / 2, C) - C / 2  # into (-C/2, C/2]

function compute_added_value(evals, targets::IntervalVector, rng::AbstractRNG=Random.GLOBAL_RNG; circular::Bool = false)
    if circular
        carcs = targets - evals                    # c-space arcs, per case
        res, depth = fold_stab(carcs)
        value = select_constant(res, rng)
        hits  = circular_hits(carcs, value)         # independent recount (validates per-case disjointness)
        value = value == 0.0 ? value : canonicalize_circular(value)  # store canonical; hits mod-C invariant
        evals = evals .+ value
    else
        res, depth = max_overlap_region((targets - evals).intervals)
        value = select_constant(res, rng)
        evals = evals .+ value
        hits  = compute_hits(evals, targets)
    end

    correct = hits >= depth 
    correct |= abs(value) > 1e+18
    @assert correct "$hits, $depth $(minimum(evals)) $(maximum(evals)) val=$value intv=$res"
    
    added_value = AddedValue(res, value, hits)

    return added_value, evals
end

function _optimize(node::BinaryNode, targets::IntervalVector, inputs, rng::AbstractRNG)
    
    right_evals_pre = evaluate(node.right, inputs)
    left_targets = inverse(node.fun, targets, right_evals_pre)

    new_left, left_evals = _optimize(node.left, left_targets, inputs, rng)

    right_targets = inverse(node.fun, left_evals, targets)
    
    new_right, right_evals_post = _optimize(node.right, right_targets, inputs, rng)
    new_evals = evaluate(node.fun, left_evals, right_evals_post)
    added_value, updated_evals = compute_added_value(new_evals, targets, rng)

    left_hits   = new_left.addition.hits
    right_hits  = new_right.addition.hits

    parent_hits = added_value.hits
    if parent_hits < left_hits || parent_hits < right_hits
        phits = compute_hits(updated_evals, targets)
        #@show node
        @show phits, parent_hits, left_hits, right_hits, node.fun, added_value.value
    end
    @assert parent_hits >= right_hits  "[DEBUG-a4f2] T3 < T2: parent=$parent_hits right=$right_hits left=$left_hits, fun=$(node.fun)"
    @assert parent_hits == added_value.hits "$parent_hits != $(added_value.hits)"

    return BinaryNode(node.fun, new_left, new_right, added_value), updated_evals
end

function _optimize(node::UnaryNode, targets::IntervalVector, inputs, rng::AbstractRNG)
    #child_targets = inverse(node.fun, targets .- node.addition.value)
    child_targets = inverse(node.fun, targets)

    new_child, child_evals = _optimize(node.child, child_targets, inputs, rng)

    new_evals = evaluate(node.fun, child_evals)
    added_value, updated_evals = compute_added_value(new_evals, targets, rng)

    parent_hits = compute_hits(updated_evals, targets)

    @assert parent_hits >= new_child.addition.hits

    return UnaryNode(node.fun, new_child, added_value), updated_evals
end

function _optimize(node::Var, targets::IntervalVector, inputs, rng::AbstractRNG)
    out = evaluate(node, inputs) .- node.addition.value
    added_value, evals = compute_added_value(out, targets, rng)
    return Var(node.index, added_value), evals
end

function _optimize(node::Constant, targets::IntervalVector, inputs, rng::AbstractRNG)
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
function optimize(root::Node, inputs, targets::IntervalVector, rng::AbstractRNG=Random.GLOBAL_RNG)
    _optimize(root, targets, inputs, rng)
end

function align_node(node::Node, targets::IntervalVector, inputs; circular::Bool=false)

    # Evaluate node without its added value
    out = evaluate(node, inputs) .- node.addition.value
    # Find an added value that optimizes hit count
    added_value, new_outputs = compute_added_value(out, targets; circular=circular)

    # An empty allowed region means the inverse chain collapsed and there is no valid
    # constant for this node: e.g. the path launders an Inf back to a finite value
    # (`x/0 = Inf`, then `exp(-Inf) = 0`), so inverting through it multiplies by
    # `inv(0) = Inf` and yields no preimage; or a sub-ULP surrogate narrowed to nothing.
    # `compute_added_value` falls back to 0.0 there, which would discard the node's
    # current (possibly hit-scoring) constant and can turn `x/c` into `x/0`. Keep the
    # node unchanged instead.
    if isempty(added_value.allowed_intervals)
        return node, out .+ node.addition.value
    end

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
