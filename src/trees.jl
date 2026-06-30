

struct AddedValue
    allowed_intervals::IntervalSet
    value::Real 

    AddedValue(v::Real) = new(IntervalSet(Interval(v, v)), v)
    AddedValue(allowed_intervals::IntervalSet{Interval{Float64,Closed,Closed}}, value::Real) = new(allowed_intervals, value)
end
zeroval = AddedValue(IntervalSet(Interval(0.0, 0.0)), 0.0)

abstract type Node end 

struct BinaryNode <: Node
    left::Node
    right::Node
    fun::Function
    addition::AddedValue
    size::Int
    function BinaryNode(fun::Function, left::Node, right::Node, addition::AddedValue)
        size = length(left) + length(right) + 1 + length(addition)
        new(left, right, fun, addition, size)
    end
    function BinaryNode(fun::Function, left::Node, right::Node) 
        size = length(left) + length(right) + 1
        new(left, right, fun, zeroval, size)
    end
end

struct Var <: Node 
    index::Int
    addition::AddedValue
    Var(index::Int, addition::AddedValue) = new(index, addition)
    Var(index::Int) = new(index, zeroval)
end

struct Constant <: Node
    addition::AddedValue
    Constant(addition::AddedValue) = new(addition)
    Constant(v::Real) = new(AddedValue(v))
    Constant() = new(zeroval)
end

struct UnaryNode <: Node
    child::Node
    fun::Function
    addition::AddedValue
    size::Int
    function UnaryNode(fun::Function, child::Node, addition::AddedValue)
        new(child, fun, addition, 1 + length(child) + length(addition))
    end
    function UnaryNode(fun::Function, child::Node)
        new(child, fun, zeroval, 1 + length(child))
    end
end

struct Tree 
    root::Node
    hits::BitVector
    Tree(root::Node, hits::BitVector) = new(root, hits)
end

Base.length(c::AddedValue) = c.value == 0 ? 0 : 1
Base.length(b::BinaryNode) = b.size
Base.length(u::UnaryNode) = u.size
Base.length(v::Var) = 1 + length(v.addition)
Base.length(z::Constant) = 1
Base.length(t::Tree) = length(t.root)

evaluate(node::BinaryNode, x) = evaluate(node.fun, evaluate(node.left, x), evaluate(node.right, x)) .+ node.addition.value
evaluate(node::UnaryNode, x) = evaluate(node.fun, evaluate(node.child, x)) .+ node.addition.value
evaluate(node::Var, x) = x[node.index] .+ node.addition.value
evaluate(node::Constant, x) = zero(first(x)) .+ node.addition.value

function evaluate(tree::Tree, x)
    try
        evaluate(tree.root, x)
    catch e
        e isa DomainError ? fill(Inf, length(x[1])) : rethrow(e)
    end
end

"""
getindex 
"""
Base.getindex(tree::Tree, i::Int) = getindex(tree.root, i)

function Base.getindex(node::Var, i::Int) 
    @assert i == 1
    return node
end
function Base.getindex(node::Constant, i::Int) 
    @assert i == 1
    return node
end

function Base.getindex(node::BinaryNode, i::Int)
    if i == 1
        return node
    elseif i <= 1 + length(node.left)
        return getindex(node.left, i - 1)
    else
        return getindex(node.right, i - 1 - length(node.left))
    end
end

function Base.getindex(node::UnaryNode, i::Int)
    i == 1 ? node : getindex(node.child, i - 1)
end

"""
insert
"""

function insert(recipient::BinaryNode, donation::Node, i::Int)
    if i == 1
        return BinaryNode(recipient.fun, donation, recipient.right, recipient.addition)
    elseif i <= 1 + length(recipient.left)
        new_left = insert(recipient.left, donation, i - 1)
        return BinaryNode(recipient.fun, new_left, recipient.right, recipient.addition)
    else
        new_right = insert(recipient.right, donation, i - 1 - length(recipient.left))
        return BinaryNode(recipient.fun, recipient.left, new_right, recipient.addition)
    end
end

insert(recipient::Var, donation::Node, i::Int) = i == 1 ? BinaryNode(+, donation, recipient) : error("Index out of bounds for Var node")
insert(recipient::Constant, donation::Node, i::Int) = i == 1 ? BinaryNode(+, donation, recipient) : error("Index out of bounds for Zero node")

function insert(recipient::UnaryNode, donation::Node, i::Int)
    if i == 1
        return UnaryNode(recipient.fun, donation, recipient.addition)
    else
        return UnaryNode(recipient.fun, insert(recipient.child, donation, i - 1), recipient.addition)
    end
end

"""
insert with alignment: insert a node into a tree at index i, and update the alignment of the recipient node to include the allowed intervals of the donation node.
It will also update the added values of all nodes upward  
"""
function insert_with_alignment(recipient::BinaryNode, donation::Node, i::Int, targets::Vector{Interval{Float64,Closed,Closed}}, inputs)

    if i == 1 
        return aligned_node(donation, targets, inputs) # find best additive constant
    elseif i <= 1 + length(recipient.left)
        # Evaluate right child -- keep
        right_evals = evaluate(recipient.right, inputs)
        # compute surrogate targets for left child
        new_targets = leftinverse(recipient.fun, targets, right_evals)
        # insert donation into left child, receive new child and left child evaluations
        new_left, left_evals = insert_with_alignment(recipient.left, donation, i - 1, new_targets, inputs)
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
        new_targets = rightinverse(recipient.fun, targets, left_evals)
        # insert donation into right child, receive new child and right child evaluations
        new_right, right_evals = insert_with_alignment(recipient.right, donation, i - 1 - length(recipient.left), new_targets, inputs)
        # Re-evaluate the new node without its added value
        new_evals = evaluate(recipient.fun, left_evals, right_evals)
        # compute an added value for the new node that aligns with the target (new_evals gets updated with that constant value)
        added_value, updated_evals = compute_added_value(new_evals, targets)
        # Return the re-optimized node with the new right child as well as the evaluations
        return BinaryNode(recipient.fun, recipient.left, new_right, added_value), updated_evals
    end
end

function insert_with_alignment(recipient::UnaryNode, donation::Node, i::Int, targets::Vector{Interval{Float64,Closed,Closed}}, inputs)
    if i == 1
        return aligned_node(donation, targets, inputs)
    else
        new_targets = inverse(recipient.fun, targets)
        new_child, child_evals = insert_with_alignment(recipient.child, donation, i - 1, new_targets, inputs)
        new_evals = evaluate(recipient.fun, child_evals)
        added_value, updated_evals = compute_added_value(new_evals, targets)
        return UnaryNode(recipient.fun, new_child, added_value), updated_evals
    end
end

function insert_with_algnment(::Union{Var, Constant}, donation::Node, i::Int, targets::Vector{Interval{Float64,Closed,Closed}}, inputs)
    if i == 1
        return aligned_node(donation, targets, inputs)
    else
        error("Index out of bounds for node")
    end
end

function compute_added_value(evals, targets::Vector{Interval{Float64,Closed,Closed}})
    res, _ = max_overlap_region(targets .- evals)
    value = select_constant(res)
    return AddedValue(res, value), evals .+ value
end

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

"""
Show
"""

function Base.show(io::IO, c::AddedValue)
    if c.value != 0.0 
        print(io, " + $(c.value)")
    end
end

Base.show(io::IO, node::Var) = print(io, "x[$(node.index)]", node.addition)
Base.show(io::IO, node::Constant) = print(io, node.addition.value)
Base.show(io::IO, node::BinaryNode) = print(io, "(", node.left, " ", node.fun, " ", node.right, ")", node.addition)
Base.show(io::IO, node::UnaryNode) = print(io, node.fun, "(", node.child, ")", node.addition)