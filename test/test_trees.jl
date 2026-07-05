using AlignedGP
using Test

import AlignedGP: evaluate

const CI = CInterval

@testset "AddedValue" begin
    av = AddedValue(3.0)
    @test av.value == 3.0
    @test complexity(av) == 1

    av_zero = AddedValue(0.0)
    @test complexity(av_zero) == 0

    region = CIntervals(CI(1.0, 2.0))
    av2 = AddedValue(region, 1.5)
    @test av2.value == 1.5
    @test av2.allowed_intervals == region
end

@testset "Node construction, length and complexity" begin
    v = Var(1)
    @test v.index == 1
    @test length(v) == 1
    @test complexity(v) == 1

    v_added = Var(2, AddedValue(1.0))
    @test length(v_added) == 1
    @test complexity(v_added) == 2

    c = Constant(5.0)
    @test length(c) == 1
    @test complexity(c) == 1

    c_zero = Constant()
    @test length(c_zero) == 1
    @test complexity(c_zero) == 1  # zero constant still counts — it IS the value

    bn = BinaryNode(+, Var(1), Var(2))
    @test length(bn) == 3
    @test complexity(bn) == 3

    bn_nested = BinaryNode(*, bn, Constant(2.0))
    @test length(bn_nested) == 5
    @test complexity(bn_nested) == 5

    # nonzero addition adds 1 to complexity but not to structural length
    bn_added = BinaryNode(+, Var(1), Var(2), AddedValue(1.0))
    @test length(bn_added) == 3
    @test complexity(bn_added) == 4
end

@testset "evaluate nodes" begin
    x = [collect(1.0:3.0), collect(4.0:6.0)]

    @test evaluate(Var(1), x) == [1.0, 2.0, 3.0]
    @test evaluate(Var(2), x) == [4.0, 5.0, 6.0]

    @test evaluate(Var(1, AddedValue(10.0)), x) == [11.0, 12.0, 13.0]

    @test evaluate(Constant(), x) == [0.0, 0.0, 0.0]
    @test evaluate(Constant(7.0), x) == [7.0, 7.0, 7.0]

    add_node = BinaryNode(+, Var(1), Var(2))
    @test evaluate(add_node, x) == [5.0, 7.0, 9.0]

    mul_node = BinaryNode(*, Var(1), Var(2))
    @test evaluate(mul_node, x) == [4.0, 10.0, 18.0]

    sub_node = BinaryNode(-, Var(2), Var(1))
    @test evaluate(sub_node, x) == [3.0, 3.0, 3.0]

    div_node = BinaryNode(/, Var(2), Var(1))
    @test evaluate(div_node, x) == [4.0, 2.5, 2.0]

    add_with_addition = BinaryNode(+, Var(1), Var(2), AddedValue(1.0))
    @test evaluate(add_with_addition, x) == [6.0, 8.0, 10.0]
end

@testset "Tree evaluate" begin
    x = [collect(1.0:3.0), collect(4.0:6.0)]

    tree = Tree(BinaryNode(+, Var(1), Var(2)), BitVector())
    @test evaluate(tree, x) == [5.0, 7.0, 9.0]
    @test length(tree) == 3
end

@testset "Tree DomainError catch" begin
    x_neg = [[-1.0, -2.0, -3.0]]
    sqrt_tree = Tree(UnaryNode(sqrt, Var(1)), BitVector())
    result = evaluate(sqrt_tree, x_neg)
    @test all(isnan, result)
    @test length(result) == 3
end

@testset "getindex" begin
    # tree: (x[1] + x[2])
    #  index 1: BinaryNode
    #  index 2: Var(1)  (left child)
    #  index 3: Var(2)  (right child)
    bn = BinaryNode(+, Var(1), Var(2))
    tree = Tree(bn, BitVector())

    @test tree[1] isa BinaryNode
    @test tree[2] isa Var && tree[2].index == 1
    @test tree[3] isa Var && tree[3].index == 2

    # Var and Constant only support index 1
    v = Var(1)
    @test v[1] === v

    c = Constant(3.0)
    @test c[1] === c

    # Nested: (x[1] + x[2]) * x[3]
    #  1: BinaryNode(*)
    #  2: BinaryNode(+)  (left subtree root)
    #  3: Var(1)
    #  4: Var(2)
    #  5: Var(3)
    nested = BinaryNode(*, BinaryNode(+, Var(1), Var(2)), Var(3))
    @test nested[1] isa BinaryNode && nested[1].fun == (*)
    @test nested[2] isa BinaryNode && nested[2].fun == (+)
    @test nested[3] isa Var && nested[3].index == 1
    @test nested[4] isa Var && nested[4].index == 2
    @test nested[5] isa Var && nested[5].index == 3
end

@testset "insert" begin
    v = Var(1)
    donated = Var(2)

    # Inserting into a leaf node at index 1 wraps it: BinaryNode(+, donation, recipient)
    result = insert(v, donated, 1)
    @test result isa BinaryNode
    @test result.fun == (+)
    @test result.left === donated
    @test result.right === v

    bn = BinaryNode(+, Var(1), Var(2))

    # Replace root (index 1): donation becomes left child, original right is kept
    result2 = insert(bn, Var(3), 1)
    @test result2 isa BinaryNode
    @test result2.left === Var(3)
    @test result2.right === bn.right

    # Replace left child (index 2): Var(1) becomes BinaryNode(+, Var(3), Var(1))
    result3 = insert(bn, Var(3), 2)
    @test result3.left isa BinaryNode
    @test result3.left.left === Var(3)
    @test result3.left.right === bn.left   # original left (Var(1)) kept as right
    @test result3.right === bn.right

    # Replace right child (index 3): Var(2) becomes BinaryNode(+, Var(3), Var(2))
    result4 = insert(bn, Var(3), 3)
    @test result4.left === bn.left
    @test result4.right isa BinaryNode
    @test result4.right.left === Var(3)

    # Nested: (x[1] + x[2]) * x[3], insert at index 3 (Var(1), the first leaf of left subtree)
    nested = BinaryNode(*, BinaryNode(+, Var(1), Var(2)), Var(3))
    result5 = insert(nested, Constant(99.0), 3)
    @test result5.left isa BinaryNode        # left subtree replaced
    @test result5.left.left isa BinaryNode   # Var(1) was wrapped into a new BinaryNode
    @test result5.left.left.left isa Constant
    @test result5.left.left.left.addition.value == 99.0
    @test result5.right === nested.right     # right subtree unchanged
end

@testset "show" begin
    @test sprint(show, Var(1)) == "x[1]"
    @test sprint(show, Var(1, AddedValue(2.0))) == "x[1] + 2.0"
    @test sprint(show, Constant(3.0)) == "3.0"
    @test sprint(show, Constant()) == "0.0"
    @test sprint(show, BinaryNode(+, Var(1), Var(2))) == "(x[1] + x[2])"
    @test sprint(show, BinaryNode(+, Var(1), Var(2), AddedValue(1.0))) == "(x[1] + x[2]) + 1.0"
end

@testset "UnaryNode construction, length and complexity" begin
    u = UnaryNode(sqrt, Var(1))
    @test u.fun === sqrt
    @test u.child === Var(1)
    @test length(u) == 2   # 1 (self) + 1 (Var)
    @test complexity(u) == 2

    u_added = UnaryNode(exp, Var(2), AddedValue(3.0))
    @test u_added.fun === exp
    @test length(u_added) == 2   # nonzero addition does not affect structural length
    @test complexity(u_added) == 3   # but does add 1 to complexity

    # UnaryNode wrapping a BinaryNode child
    u_deep = UnaryNode(log, BinaryNode(+, Var(1), Var(2)))
    @test length(u_deep) == 4   # 1 + 3
    @test complexity(u_deep) == 4
end

@testset "UnaryNode evaluate" begin
    x = [collect(1.0:3.0)]

    @test evaluate(UnaryNode(sqrt, Var(1)), x) ≈ sqrt.([1.0, 2.0, 3.0])
    @test evaluate(UnaryNode(exp, Var(1)), x) ≈ exp.([1.0, 2.0, 3.0])
    @test evaluate(UnaryNode(log, Var(1)), x) ≈ log.([1.0, 2.0, 3.0])

    # addition is applied after the unary function
    u_shift = UnaryNode(sqrt, Var(1), AddedValue(10.0))
    @test evaluate(u_shift, x) ≈ sqrt.([1.0, 2.0, 3.0]) .+ 10.0

    # nested: sqrt(x[1] + x[2])
    x2 = [collect(1.0:3.0), collect(3.0:-1.0:1.0)]
    u_nested = UnaryNode(sqrt, BinaryNode(+, Var(1), Var(2)))
    @test evaluate(u_nested, x2) ≈ sqrt.([4.0, 4.0, 4.0])
end

@testset "UnaryNode DomainError via Tree" begin
    x_neg = [[-1.0, -2.0, -3.0]]
    @test all(isnan, evaluate(Tree(UnaryNode(sqrt, Var(1)), BitVector()), x_neg))
    @test all(isnan, evaluate(Tree(UnaryNode(log, Var(1)), BitVector()), x_neg))
end

@testset "UnaryNode getindex" begin
    u = UnaryNode(sqrt, Var(1))
    # index 1: the UnaryNode itself
    @test u[1] === u
    # index 2: the child
    @test u[2] isa Var && u[2].index == 1

    # Deeper: sqrt((x[1] + x[2]))
    #  1: UnaryNode(sqrt)
    #  2: BinaryNode(+)
    #  3: Var(1)
    #  4: Var(2)
    u_deep = UnaryNode(sqrt, BinaryNode(+, Var(1), Var(2)))
    @test u_deep[1] === u_deep
    @test u_deep[2] isa BinaryNode
    @test u_deep[3] isa Var && u_deep[3].index == 1
    @test u_deep[4] isa Var && u_deep[4].index == 2
end

@testset "UnaryNode insert" begin
    u = UnaryNode(sqrt, Var(1))

    # index 1: donation replaces the child, fun and addition are preserved
    result = insert(u, Var(2), 1)
    @test result isa UnaryNode
    @test result.fun === sqrt
    @test result.child === Var(2)

    # index 2: descend into child (Var(1) gets wrapped)
    result2 = insert(u, Var(3), 2)
    @test result2 isa UnaryNode
    @test result2.fun === sqrt
    @test result2.child isa BinaryNode
    @test result2.child.left === Var(3)
    @test result2.child.right === u.child

    # addition is preserved across insert
    u_added = UnaryNode(exp, Var(1), AddedValue(5.0))
    result3 = insert(u_added, Var(2), 1)
    @test result3.addition.value == 5.0
end

@testset "UnaryNode show" begin
    @test sprint(show, UnaryNode(sqrt, Var(1))) == "sqrt(x[1])"
    @test sprint(show, UnaryNode(exp, Var(2))) == "exp(x[2])"
    @test sprint(show, UnaryNode(log, Var(1), AddedValue(3.0))) == "log(x[1]) + 3.0"
    @test sprint(show, UnaryNode(sqrt, BinaryNode(+, Var(1), Var(2)))) == "sqrt((x[1] + x[2]))"
end
