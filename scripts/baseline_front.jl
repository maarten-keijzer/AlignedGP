using AlignedGP

struct BaseIndy 
    root::Node 
    slope::Float64 
    intercept::Float64 
    mse::Float64
    previous::Union{Nothing, BaseIndy}

    length::Int 

    function BaseIndy(root, slope, intercept, mse, previous) 
        # (b1 X + c1) + (b2 Y + c2) = b1 X + b2 Y + c3 
        len = isnothing(previous) ? length(root) : length(root) + 2 + length(previous)
        new(root, slope, intercept, mse, previous, len)
    end
end

Base.length(indy::BaseIndy) = indy.length

function create_indy(parent, node, evals, targets)
    slope, intercept, mse = linear_scale(evals, targets)
    BaseIndy(node, slope, intercept, mse, parent)
end

function baseline_fit(parent, symboltable, inputs, targets)
    node, evals = valid_init(symboltable, rand(3:13), inputs)
    all(isfinite, evals) || return nothing, nothing
    indy = create_indy(parent, node, evals, targets)
    return indy, indy.slope .* evals .+ indy.intercept
end


function create_baseline_front(setup, minerror)
    targets = copy(setup.noisy_targets)
    inputs = setup.inputs 

    indy = argmin( indy -> indy.mse, create_indy(nothing, Var(i), evaluate(Var(i), inputs), targets) for i in 1:length(inputs))
    evals = evaluate(indy.root, inputs)

    getmse(indy::BaseIndy) = indy.mse
    front = Front{BaseIndy}(errfn=getmse)

    while indy.mse > minerror 
        targets .-= evals
        best, scaled_evals = baseline_fit(indy, setup.symboltable, inputs, targets)
        for _ in 2:50 
            contender, scaled_evals2 = baseline_fit(indy, setup.symboltable, inputs, targets)
            if isnothing(best) || (!isnothing(contender) && contender.mse < best.mse)
                best = contender
                scaled_evals = scaled_evals2 
            end     
        end
        indy = best 
        evals = scaled_evals
        add_to_front!(front, indy)
        if length(front.front) % 20 == 0 
            @show indy.mse, length(indy)
        end
    end
    return front
end

