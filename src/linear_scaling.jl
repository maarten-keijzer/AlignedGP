
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
