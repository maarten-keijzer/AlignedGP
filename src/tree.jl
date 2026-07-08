
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


Base.length(t::Tree) = length(t.root)
complexity(t::Tree) = t.complexity
pathlen_complexity(t::Tree) = t.pathlen_complexity
evaluate(tree::Tree, x) = evaluate(tree.root, x)
scaled_evaluate(tree::Tree, x) = tree.slope .* evaluate(tree.root, x) .+ tree.intercept
Base.getindex(tree::Tree, i::Int) = getindex(tree.root, i)

function evaluate_to_tree(node::Node, setup::ProblemSetup)
    inputs = setup.inputs
    targets = setup.interval_targets

    output = evaluate(node, inputs)
    hits = BitVector(output[i] ∈ targets[i] for i in 1:length(targets))
    
    if setup.params.use_l2_scaling
        slope, intercept, mse = linear_scale(output, setup.noisy_targets)
    else 
        slope, intercept, mse, _ = linear_scale_l1(output, setup.noisy_targets)
    end
    scaled_out = slope .* output .+ intercept
    
    scaled_hits = BitVector(scaled_out[i] ∈ targets[i] for i in eachindex(targets))
    scaled_hitcount = sum(scaled_hits)

    return Tree(node, scaled_hitcount > sum(hits) ? scaled_hits : hits, slope, intercept, mse, scaled_hitcount)
end




