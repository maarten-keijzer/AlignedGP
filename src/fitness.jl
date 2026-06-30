
struct SizeStratum
    pop :: Vector{Tree}
end

function Tree(node::Node, targets::Vector{Interval{Float64,Closed,Closed}}, inputs)
   output = evaluate(node, inputs)
   hits::BitVector = [output[i] ∈ targets[i] for i in 1:length(targets)]
   return Tree(node, hits)
end



