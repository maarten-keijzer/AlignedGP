
function linear_scale_l1(output::Vector{Float64}, targets::Vector{Float64})
    n = length(output)
    p_mean = sum(output) / n

    vp = sum((output[i] - p_mean)^2 for i in 1:n) / n
    if vp < eps()
        b = median(targets)
        mae = sum(abs(targets[i] - b) for i in 1:n) / n
        return 0.0, b, mae, 0
    end

    # initialise with OLS solution
    slope, intercept, _ = linear_scale(output, targets)

    local iter
    for outer iter in 1:50
        weights = [1.0 / max(abs(slope * output[i] + intercept - targets[i]), eps()) for i in 1:n]
        sw  = sum(weights)
        swx = sum(weights[i] * output[i] for i in 1:n)
        swy = sum(weights[i] * targets[i] for i in 1:n)
        swxx = sum(weights[i] * output[i]^2 for i in 1:n)
        swxy = sum(weights[i] * output[i] * targets[i] for i in 1:n)

        denom = sw * swxx - swx^2
        new_slope     = abs(denom) < eps() ? 0.0 : (sw * swxy - swx * swy) / denom
        new_intercept = (swy - new_slope * swx) / sw

        converged = abs(new_slope - slope) < 1e-8 && abs(new_intercept - intercept) < 1e-8
        slope, intercept = new_slope, new_intercept
        converged && break
    end

    mae = sum(abs(slope * output[i] + intercept - targets[i]) for i in 1:n) / n
    return slope, intercept, mae, iter
end

function linear_scale(output::Vector{Float64}, targets::Vector{Float64})
    n = length(output)
    p_mean = sum(output) / n
    t_mean = sum(targets) / n

    vp = sum((output[i] - p_mean)^2 for i in 1:n) / n

    slope, intercept = if vp < eps()
        0.0, t_mean
    else
        cov_pt = sum((output[i] - p_mean) * (targets[i] - t_mean) for i in 1:n) / n
        a = cov_pt / vp
        a, t_mean - a * p_mean
    end

    mse = sum((slope * output[i] + intercept - targets[i])^2 for i in 1:n) / n
    return slope, intercept, mse
end
