using AlignedGP
using Test

import AlignedGP.draw_subtree

@testset "Draw node" begin
    tree = BinaryNode(+, Constant(), Constant())
    tree = BinaryNode(+, tree, tree)
    tree = BinaryNode(+, tree, tree)
    tree = BinaryNode(+, tree, tree)
        
    nterminals = 0
    for _ in 1:1000
        sub, point = draw_subtree(tree)
        if length(sub) == 1
            nterminals += 1
        end
    end
    @test nterminals < 130
end


