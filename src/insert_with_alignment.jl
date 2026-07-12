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
        #new_targets = _surrogate_left(recipient.fun, targets, right_evals)        

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
        #new_targets = _surrogate_right(recipient.fun, targets, left_evals)

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
        #new_targets = inverse(recipient.fun, targets)

        new_child, child_evals, d = insert_with_alignment(recipient.child, donation, i - 1, d+1, inputs, new_targets, recursive_stabbing)

        if recursive_stabbing && is_periodic(recipient.fun)
            raw_child = child_evals .- new_child.addition.value
            circ_av, child_evals = circular_compute_added_value(raw_child, targets .- recipient.addition.value, recipient.fun)
            new_child = set_addition(new_child, circ_av)
        end

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
    # `targets .- evals` narrows (CInterval - Real → narrow), so a point target [x,x]
    # with eval==x collapses to invalid_interval and drops out of the depth count; the
    # exact hit is still recovered below via `unnarrowed_hits` against the raw targets.
    res, depth = max_overlap_region(targets .- evals)
    value = select_constant(res, rng)
    evals = evals .+ value

    unnarrowed_hits = _hits(evals, targets)
    @assert unnarrowed_hits >= depth "$(_hits(evals, targets)), $depth $(minimum(evals)) $(maximum(evals)) val=$value intv=$res"

    added_value = AddedValue(res, value, unnarrowed_hits)
    
    return added_value, evals
end

# Reconstruct any node type with a replacement AddedValue (nodes are immutable).
set_addition(n::Var,        av::AddedValue) = Var(n.index, av)
set_addition(n::Constant,   av::AddedValue) = Constant(av)
set_addition(n::UnaryNode,  av::AddedValue) = UnaryNode(n.fun, n.child, av)
set_addition(n::BinaryNode, av::AddedValue) = BinaryNode(n.fun, n.left, n.right, av)

# Circular stab: find the additive constant c ∈ [0,2π) that maximises hits when
# fn(raw_evals[i] + c) ∈ targets[i].  Uses the mod-2π alignment + unrolling trick
# described in specs/sine-inversion.md.
function circular_compute_added_value(raw_evals, targets::Vector{CIntervals}, fn, rng::AbstractRNG=Random.GLOBAL_RNG)
    TWO_PI = 2π
    arcs = Tuple{Float64,Float64}[]

    for i in eachindex(raw_evals)
        isfinite(raw_evals[i]) || continue
        s_mod = mod(raw_evals[i], TWO_PI)
        for ci in targets[i].items
            for (a, b) in _preimage_circular(fn, ci)
                c_s = mod(a - s_mod, TWO_PI)
                c_e = c_s + (b - a)           # preserve arc length
                if c_e <= TWO_PI
                    push!(arcs, (c_s, c_e))
                    # Duplicate in [2π,4π). If c_s==0 the duplicate would start at exactly
                    # 2π and collide (enter before leave) with any arc ending at 2π — use
                    # nextfloat to break the tie.
                    dup_start = c_s > 0.0 ? c_s + TWO_PI : nextfloat(TWO_PI)
                    push!(arcs, (dup_start, c_e + TWO_PI))
                else
                    # wrapping arc: split, then duplicate each piece.
                    # Piece 4 must start at nextfloat(2π) so its enter event fires AFTER
                    # piece 1's leave event at 2π (same-coordinate: leave < enter sort order
                    # would otherwise inflate depth by 2× at the seam).
                    push!(arcs, (c_s,             TWO_PI));          push!(arcs, (c_s + TWO_PI, 2*TWO_PI))
                    push!(arcs, (0.0,              c_e - TWO_PI));   push!(arcs, (nextfloat(TWO_PI), c_e))
                end
            end
        end
    end

    if isempty(arcs)
        return AddedValue(CIntervals(), 0.0, 0), copy(raw_evals)
    end

    ci_list = CIntervals[CIntervals(CInterval(a, b)) for (a, b) in arcs]
    region_unrolled, depth = max_overlap_region(ci_list)

    # Project optimal region back to [0, 2π): keep [0,2π) slice and map [2π,4π)→[0,2π).
    parts = CInterval[]
    for ci in region_unrolled.items
        lo1 = max(0.0, ci.lo);    hi1 = min(TWO_PI, ci.hi);    lo1 < hi1 && push!(parts, CInterval(lo1, hi1))
        lo2 = max(TWO_PI, ci.lo) - TWO_PI; hi2 = min(2*TWO_PI, ci.hi) - TWO_PI
        lo2 < hi2 && push!(parts, CInterval(lo2, hi2))
    end
    sort!(parts; by=x -> x.lo)
    merged = CInterval[]
    for ci in parts
        if isempty(merged) || ci.lo > last(merged).hi
            push!(merged, ci)
        else
            merged[end] = CInterval(last(merged).lo, max(last(merged).hi, ci.hi))
        end
    end

    region = CIntervals(merged)
    value  = select_constant(region, rng)
    evals  = raw_evals .+ value
    return AddedValue(region, value, depth), evals
end

# Map non-finite sibling evals to the invalid sentinel so they don't corrupt the
# interval inverse computation (NaN is not a valid CInterval).
_surrogate_left(fun, targets, se)  = CIntervals[isfinite(se[j]) ? leftinverse(fun, targets[j], se[j])  : CIntervals() for j in eachindex(targets)]
_surrogate_right(fun, targets, se) = CIntervals[isfinite(se[j]) ? rightinverse(fun, targets[j], se[j]) : CIntervals() for j in eachindex(targets)]

function _optimize(node::BinaryNode, targets::Vector{CIntervals}, inputs, rng::AbstractRNG)
    right_evals_pre = evaluate(node.right, inputs)
    #left_targets = _surrogate_left(node.fun, targets .- node.addition.value, right_evals_pre)
    left_targets = _surrogate_left(node.fun, targets, right_evals_pre)

    new_left, left_evals = _optimize(node.left, left_targets, inputs, rng)

    #right_targets = _surrogate_right(node.fun, targets .- node.addition.value, left_evals)
    right_targets = _surrogate_right(node.fun, targets, left_evals)
    
    new_right, right_evals_post = _optimize(node.right, right_targets, inputs, rng)
    new_evals = evaluate(node.fun, left_evals, right_evals_post)
    added_value, updated_evals = compute_added_value(new_evals, targets, rng)

    left_hits   = new_left.addition.hits
    right_hits  = new_right.addition.hits

    parent_hits = added_value.hits
    if parent_hits < left_hits || parent_hits < right_hits
        phits = _hits(updated_evals, targets)
        @show phits, parent_hits, left_hits, right_hits, node.fun
    end
    @assert parent_hits >= right_hits  "[DEBUG-a4f2] T3 < T2: parent=$parent_hits right=$right_hits left=$left_hits, fun=$(node.fun)"
    @assert parent_hits == added_value.hits "$parent_hits != $(added_value.hits)"

    return BinaryNode(node.fun, new_left, new_right, added_value), updated_evals
end

function _optimize(node::UnaryNode, targets::Vector{CIntervals}, inputs, rng::AbstractRNG)
    #child_targets = inverse(node.fun, targets .- node.addition.value)
    child_targets = inverse(node.fun, targets)

    new_child, child_evals = _optimize(node.child, child_targets, inputs, rng)

    if is_periodic(node.fun)
        raw_child = child_evals .- new_child.addition.value
        circ_av, child_evals = circular_compute_added_value(raw_child, targets, node.fun, rng)
        new_child = set_addition(new_child, circ_av)
    end

    new_evals = evaluate(node.fun, child_evals)

    # Sync circular-stab depth with actual fn evaluation hits.  FP boundary effects can
    # make depth > _hits(fn(child_evals), targets) by 1; parent_hits (which adds its own
    # constant on top) must be >= new_child.addition.hits, and parent_hits >= actual_sin_hits.
    if is_periodic(node.fun)
        actual_sin_hits = _hits(new_evals, targets)
        if actual_sin_hits < new_child.addition.hits
            new_child = set_addition(new_child, AddedValue(new_child.addition.allowed_intervals, new_child.addition.value, actual_sin_hits))
        end
    end

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

