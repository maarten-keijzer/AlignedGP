

using Random

struct AddedValue
    allowed_intervals::IntervalSet
    value::Real

    function AddedValue(v::Real) 
        @assert isfinite(v)
        new(IntervalSet(Interval(v, v)), v)
    end
    AddedValue(allowed_intervals::IntervalSet{Interval{Float64,Closed,Closed}}, value::Real) = new(allowed_intervals, value)
end
zeroval = AddedValue(IntervalSet(Interval(0.0, 0.0)), 0.0)

complexity(c::AddedValue) = c.value == 0 ? 0 : 1

abstract type Node end

struct BinaryNode <: Node
    left::Node
    right::Node
    fun::Function
    addition::AddedValue
    size::Int
    complexity::Int
    function BinaryNode(fun::Function, left::Node, right::Node, addition::AddedValue)
        size = length(left) + length(right) + 1
        cmplx = complexity(left) + complexity(right) + 1 + complexity(addition)
        new(left, right, fun, addition, size, cmplx)
    end
    function BinaryNode(fun::Function, left::Node, right::Node)
        size = length(left) + length(right) + 1
        cmplx = complexity(left) + complexity(right) + 1
        new(left, right, fun, zeroval, size, cmplx)
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
    complexity::Int
    function UnaryNode(fun::Function, child::Node, addition::AddedValue)
        new(child, fun, addition, 1 + length(child), 1 + complexity(child) + complexity(addition))
    end
    function UnaryNode(fun::Function, child::Node)
        new(child, fun, zeroval, 1 + length(child), 1 + complexity(child))
    end
end

struct Tree
    root::Node
    hits::BitVector
    complexity::Int
    pathlen_complexity::Int
    slope::Float64
    intercept::Float64
    mse::Float64
    scaled_hitcount::Int
    Tree(root::Node, hits::BitVector) = new(root, hits, complexity(root), pathlen_complexity(root), 1.0, 0.0, Inf, 0)
    Tree(root::Node, hits::BitVector, slope::Float64, intercept::Float64, mse::Float64, scaled_hitcount::Int) =
        new(root, hits, complexity(root), pathlen_complexity(root), slope, intercept, mse, scaled_hitcount)
end

Base.length(b::BinaryNode) = b.size
Base.length(u::UnaryNode) = u.size
Base.length(v::Var) = 1
Base.length(z::Constant) = 1
Base.length(t::Tree) = length(t.root)

complexity(b::BinaryNode) = b.complexity
complexity(u::UnaryNode) = u.complexity
complexity(v::Var) = 1 + complexity(v.addition)
complexity(::Constant) = 1
complexity(t::Tree) = t.complexity

pathlen_complexity(b::BinaryNode) = complexity(b) + pathlen_complexity(b.left) + pathlen_complexity(b.right)
pathlen_complexity(u::UnaryNode) = complexity(u) + pathlen_complexity(u.child)
pathlen_complexity(v::Var) = complexity(v)
pathlen_complexity(::Constant) = 1
pathlen_complexity(t::Tree) = t.pathlen_complexity

evaluate(node::BinaryNode, x) = evaluate(node.fun, evaluate(node.left, x), evaluate(node.right, x)) .+ node.addition.value
evaluate(node::UnaryNode, x) = evaluate(node.fun, evaluate(node.child, x)) .+ node.addition.value
evaluate(node::Var, x) = x[node.index] .+ node.addition.value
evaluate(node::Constant, x) = zero(first(x)) .+ node.addition.value

evaluate(tree::Tree, x) = evaluate(tree.root, x)
scaled_evaluate(tree::Tree, x) = tree.slope .* evaluate(tree.root, x) .+ tree.intercept

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
Show
"""

function Base.show(io::IO, c::AddedValue)
    if c.value != 0.0
        # itv = c.allowed_intervals.items[1]
        # lo = first(itv)
        # hi = last(itv)
        #print(io, " + <$lo, $hi> $(c.value)")
        print(io, " + $(c.value)")
    end
end

Base.show(io::IO, node::Var) = print(io, "x[$(node.index)]", node.addition)
Base.show(io::IO, node::Constant) = print(io, node.addition.value)
Base.show(io::IO, node::BinaryNode) = print(io, "(", node.left, " ", node.fun, " ", node.right, ")", node.addition)
Base.show(io::IO, node::UnaryNode) = print(io, node.fun, "(", node.child, ")", node.addition)

