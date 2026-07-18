
struct Tree
    root::Node
    
    hits::BitVector             # primary (inner-band) hits
    secondary_hits::BitVector   # outer-band hits; ⊇ primary under the containment invariant
    
    loss::Float64               # two_band_loss(hits, secondary_hits); lower is better
    isscaledloss::Bool          # Is the loss obtained through scaling?    
    
    complexity::Int
    pathlen_complexity::Int
    
    slope::Float64
    intercept::Float64
    mse::Float64
    
    evals::Vector{Float64}      # cached raw node output; tol-invariant, drives retarget
    
    #Lightweight wrapper (e.g. for `evaluate` and testing) — no fitness computed yet.
    Tree(root::Node, hits::BitVector) =
        new(root, hits, BitVector(), Inf, false, complexity(root), pathlen_complexity(root),
            1.0, 0.0, Inf, Float64[])
    Tree(root::Node, hits::BitVector, secondary_hits::BitVector, loss::Float64, isscaledloss::Bool,
         slope::Float64, intercept::Float64, mse::Float64,
         evals::Vector{Float64}) =
        new(root, hits, secondary_hits, loss, isscaledloss, complexity(root), pathlen_complexity(root),
            slope, intercept, mse, evals)
end


Base.length(t::Tree) = length(t.root)
complexity(t::Tree) = t.complexity
pathlen_complexity(t::Tree) = t.pathlen_complexity
evaluate(tree::Tree, x) = evaluate(tree.root, x)
model_evaluations(tree::Tree) = tree.isscaledloss ? tree.slope .* tree.evals .+ tree.intercept : tree.evals
Base.getindex(tree::Tree, i::Int) = getindex(tree.root, i)

function hitvector(ev::Vector{<:Real}, t::IntervalVector)
    b = Vector{Bool}(undef, length(t))
    for i in eachindex(t)
        slice = t[i]
        b[i] = any(intv -> in_interval(ev[i], intv), slice)
    end
    BitVector(b)
end

function compute_hits(ev::Vector{<:Real}, t::IntervalVector)
    hits = 0
    for i in eachindex(t)
        slice = t[i]
        for elem in slice 
            if in_interval(ev[i], elem)
                hits+=1
            end
        end
    end
    return hits
end

function evaluate_to_tree(node::Node, setup::ProblemSetup)
    output = evaluate(node, setup.inputs)

    if setup.params.use_l2_scaling
        slope, intercept, mse = linear_scale(output, setup.noisy_targets)
    else
        slope, intercept, mse, _ = linear_scale_l1(output, setup.noisy_targets)
    end

    # Cache the raw output on the tree; `retarget` rebuilds hits from it (and the
    # tol-invariant slope/intercept) without re-evaluating the node.
    return _tree_from_evals(node, output, slope, intercept, mse, setup)
end

# Rebuild a tree's two-band fitness against the current bands from its cached raw
# evals, reusing slope/intercept/mse (all tol-invariant). This is what `set_tol!`
# calls after a ratchet, replacing a full re-evaluation.
function retarget(tree::Tree, setup::ProblemSetup)
    _tree_from_evals(tree.root, tree.evals, tree.slope, tree.intercept, tree.mse, setup)
end

# Shared core: pick the better of the raw vs linearly-scaled eval series by primary
# hit count, keeping that series' secondary hits (so the containment invariant is never
# broken by mixing series), then compute the two-band loss.
function _tree_from_evals(node::Node, output::Vector{Float64},
                          slope::Float64, intercept::Float64, mse::Float64,
                          setup::ProblemSetup)
    targets = setup.interval_targets
    noisy = setup.noisy_targets
    tau_outer = setup.params.tau_outer

    p_raw, s_raw = two_band_hits(output, targets, noisy, tau_outer)
    raw_loss = two_band_loss(p_raw, s_raw)

    scaled_out = slope .* output .+ intercept
    p_sc, s_sc = two_band_hits(scaled_out, targets, noisy, tau_outer)
    scaled_loss = two_band_loss(p_sc, s_sc)

    isloss_scaled = scaled_loss < raw_loss
    loss = isloss_scaled ? scaled_loss : raw_loss
    hits = isloss_scaled ? p_sc : p_raw 
    secondary = isloss_scaled ? s_sc : s_raw

    return Tree(node, hits, secondary, loss, isloss_scaled, slope, intercept, mse, output)
end




