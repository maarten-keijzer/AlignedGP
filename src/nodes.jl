

using Random

struct AddedValue
    allowed_intervals::Vector{IntervalType}
    value::Real
    hits::Int

    function AddedValue(v::Real)
        @assert isfinite(v)
        new([], v, -1)
    end
    AddedValue(allowed_intervals::Vector{IntervalType}, value::Real) = new(allowed_intervals, value, -1)
    AddedValue(allowed_intervals::Vector{IntervalType}, value::Real, hits::Int) = new(allowed_intervals, value, hits)
end
zeroval = AddedValue(0.0)

complexity(c::AddedValue) = c.value == 0 ? 0 : 1

abstract type Node end

struct BinaryNode <: Node
    left::Node
    right::Node
    fun::Function
    addition::AddedValue
    size::Int
    complexity::Int
    pathlen_complexity::Int
    function BinaryNode(fun::Function, left::Node, right::Node, addition::AddedValue)
        size = length(left) + length(right) + 1
        cmplx = complexity(left) + complexity(right) + 1 + complexity(addition)
        plc = cmplx + pathlen_complexity(left) + pathlen_complexity(right)
        new(left, right, fun, addition, size, cmplx, plc)
    end
    function BinaryNode(fun::Function, left::Node, right::Node)
        size = length(left) + length(right) + 1
        cmplx = complexity(left) + complexity(right) + 1
        plc = cmplx + pathlen_complexity(left) + pathlen_complexity(right)
        new(left, right, fun, zeroval, size, cmplx, plc)
    end
end

struct Var <: Node
    index::Int
    addition::AddedValue
    Var(index::Int, addition::AddedValue) = new(index, addition)
    Var(index::Int, c::Float64) = new(index, AddedValue(c))
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
    pathlen_complexity::Int
    function UnaryNode(fun::Function, child::Node, addition::AddedValue)
        cmplx = 1 + complexity(child) + complexity(addition)
        new(child, fun, addition, 1 + length(child), cmplx, cmplx + pathlen_complexity(child))
    end
    function UnaryNode(fun::Function, child::Node)
        cmplx = 1 + complexity(child)
        new(child, fun, zeroval, 1 + length(child), cmplx, cmplx + pathlen_complexity(child))
    end
end

Base.length(b::BinaryNode) = b.size
Base.length(u::UnaryNode) = u.size
Base.length(v::Var) = 1
Base.length(z::Constant) = 1

complexity(b::BinaryNode) = b.complexity
complexity(u::UnaryNode) = u.complexity
complexity(v::Var) = 1 + complexity(v.addition)
complexity(::Constant) = 1

pathlen_complexity(b::BinaryNode) = b.pathlen_complexity
pathlen_complexity(u::UnaryNode) = u.pathlen_complexity
pathlen_complexity(v::Var) = complexity(v)
pathlen_complexity(::Constant) = 1

evaluate(node::BinaryNode, x) = evaluate(node.fun, evaluate(node.left, x), evaluate(node.right, x)) .+ node.addition.value
evaluate(node::UnaryNode, x) = evaluate(node.fun, evaluate(node.child, x)) .+ node.addition.value
evaluate(node::Var, x) = x[node.index] .+ node.addition.value
evaluate(node::Constant, x) = zero(first(x)) .+ node.addition.value

"""
getindex
"""

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

insert(node::Node, donation::Node, i::Int) = insert(node, donation, i, 1) |> first

function insert(recipient::BinaryNode, donation::Node, i::Int, d::Int)
    if i == 1
        return BinaryNode(recipient.fun, donation, recipient.right, recipient.addition), d
    elseif i <= 1 + length(recipient.left)
        new_left, d = insert(recipient.left, donation, i - 1, d+1)
        return BinaryNode(recipient.fun, new_left, recipient.right, recipient.addition), d
    else
        new_right, d = insert(recipient.right, donation, i - 1 - length(recipient.left), d+1)
        return BinaryNode(recipient.fun, recipient.left, new_right, recipient.addition), d
    end
end

insert(recipient::Var, donation::Node, i::Int, d::Int) = BinaryNode(+, donation, recipient), d
insert(recipient::Constant, donation::Node, i::Int, d::Int) = BinaryNode(+, donation, recipient), d

function insert(recipient::UnaryNode, donation::Node, i::Int, d::Int)
    if i == 1
        return UnaryNode(recipient.fun, donation, recipient.addition), d
    else
        new_node, d = insert(recipient.child, donation, i-1, d+1)
        return UnaryNode(recipient.fun, new_node, recipient.addition), d
    end
end

"""
Collect nodes
"""
function collect_nodes!(node::Node, nodes::Vector{Node}) 
    push!(nodes, node)
    if node isa BinaryNode 
        collect_nodes!(node.left, nodes)
        collect_nodes!(node.right, nodes)
    elseif node isa UnaryNode 
        collect_nodes!(node.child, nodes)
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

export print_constructive
print_constructive(v::Var) = "Var($(v.index), $(v.addition.value))"
print_constructive(c::Constant) = "Constant($(c.addition.value))"
print_constructive(u::UnaryNode) = "UnaryNode($(u.fun), $(print_constructive(u.child)), AddedValue($(u.addition.value)))"
print_constructive(b::BinaryNode) = "BinaryNode($(b.fun), $(print_constructive(b.left)), $(print_constructive(b.right)), AddedValue($(b.addition.value)))"
