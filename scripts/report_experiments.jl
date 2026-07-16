using AlignedGP
using CairoMakie

const METHOD_MAP = Dict(
    "Standard"      => Standard,
    "Stab"          => Stab,
    "RecursiveStab" => RecursiveStab,
    "ConstantStab"  => ConstantStab,
)

function parse_params(comments)
    get(k) = strip(last(split(comments[k], "=")))
    GPParams(
        method                    = METHOD_MAP[get("method")],
        population_size           = parse(Int,     get("population_size")),
        max_complexity            = parse(Int,     get("max_complexity")),
        cross_mut_prob            = parse(Float64, get("cross_mut_prob")),
        max_lexicase_comparisons  = parse(Int,     get("max_lexicase_comparisons")),
        use_l2_scaling            = parse(Bool,    get("use_l2_scaling")),
        constant_stab_probability = parse(Float64, get("constant_stab_probability")),
        use_tournament_stratum    = parse(Bool,    get("use_tournament_stratum")),
    )
end

function read_run(path)
    comments = Dict{String,String}()
    rows = Vector{NTuple{5,Float64}}()
    for line in eachline(path)
        if startswith(line, "#")
            k, v = split(lstrip(line, '#'), "="; limit=2)
            comments[strip(k)] = strip(v)
        elseif startswith(line, "time_running")
            continue  # header
        else
            parts = split(line, ",")
            length(parts) < 3 && continue
            t    = parse(Float64, parts[1])
            e    = parse(Float64, parts[2]) / 1_000_000 #Mnods
            h    = parse(Float64, parts[3])
            iter = length(parts) >= 4 ? parse(Float64, parts[4]) : NaN
            mse  = length(parts) >= 5 ? parse(Float64, parts[5]) : NaN
            push!(rows, (t, e, h, iter, mse))
        end
    end
    comments, rows
end

function read_setup(dir)
    experiments = Dict{String, @NamedTuple{params::GPParams, n_targets::Int, runs::Vector{Matrix{Float64}}}}()

    for fname in readdir(dir; join=true)
        endswith(fname, ".csv") || continue
        base = basename(fname)
        m = match(r"^([0-9a-f]{8})_\d+\.csv$", base)
        isnothing(m) && continue
        hash_key = m[1]

        comments, rows = read_run(fname)
        isempty(rows) && continue

        data = Matrix{Float64}(undef, length(rows), 5)
        for (i, (t, e, h, iter, mse)) in enumerate(rows)
            data[i, :] = [t, e, h, iter, mse]
        end

        if !haskey(experiments, hash_key)
            params = parse_params(comments)
            n_targets = parse(Int, strip(comments["n_targets"]))
            experiments[hash_key] = (; params, n_targets, runs=Matrix{Float64}[])
        end
        push!(experiments[hash_key].runs, data)
    end

    experiments
end

# xcol: column index for x-axis values
# ycol: column index for y-axis values (3=hits, 5=best_mse)
# y_agg: how to aggregate multiple y values in a bin (maximum for hits, minimum for mse)
# y_transform: applied after aggregation (e.g. success rate threshold)
function avg_by_effort(runs, bin_edges, xcol, ycol, y_agg, y_transform)
    nbins = length(bin_edges) - 1
    avg = zeros(Float64, nbins)
    counts = zeros(Int, nbins)

    for data in runs
        effs = data[:, xcol]
        yvals = data[:, ycol]
        last_y = NaN
        for b in 1:nbins
            lo, hi = bin_edges[b], bin_edges[b+1]
            mask = (effs .>= lo) .& (effs .< hi)
            if any(mask)
                finite_vals = filter(isfinite, yvals[mask])
                isempty(finite_vals) || (last_y = y_agg(finite_vals))
            end
            isnan(last_y) && continue
            avg[b] += y_transform(last_y)
            counts[b] += 1
        end
    end

    bin_centers = [(bin_edges[b] + bin_edges[b+1]) / 2 for b in 1:nbins]
    valid = counts .> 0
    bin_centers[valid], avg[valid] ./ counts[valid]
end

function plot_experiments(experiments; nbins=50, log_scale=true, x_axis=:effort, y_axis=:hits,
                          output="experiments.pdf")
    xcol = x_axis == :time ? 1 : x_axis == :iterations ? 4 : 2
    xlabel = x_axis == :time ? "Time (s)" : x_axis == :iterations ? "Individuals processed" : "Effort (M nods)"
    ylabel = y_axis == :success_rate ? "Success rate (%)" :
             y_axis == :best_mse     ? "Best MSE"          : "Average max hits"

    all_x = Float64[]
    for exp in values(experiments)
        for data in exp.runs
            append!(all_x, data[:, xcol])
        end
    end
    isempty(all_x) && error("No data found")

    pos = all_x[all_x .> 0]
    x_min = isempty(pos) ? minimum(all_x) : (log_scale ? 10^floor(log10(minimum(pos))) : minimum(pos))
    x_max = log_scale ? 10^ceil(log10(maximum(all_x))) : maximum(all_x)
    bin_edges = log_scale ?
        10 .^ range(log10(x_min), log10(x_max); length=nbins+1) :
        range(x_min, x_max; length=nbins+1)

    fig = Figure(size=(500, 400), fontsize=14)
    ax = Axis(fig[1, 1];
        xlabel=xlabel, ylabel=ylabel,
        xscale=log_scale ? log10 : identity,
        title="Experiment comparison",
    )

    for (hash_key, exp) in experiments
        if y_axis == :success_rate
            ycol, y_agg = 3, maximum
            y_transform = h -> h >= exp.n_targets ? 100.0 : 0.0
        elseif y_axis == :best_mse
            ycol, y_agg = 5, minimum
            y_transform = identity
        else
            ycol, y_agg = 3, maximum
            y_transform = identity
        end

        label = exp.params.method == Stab ? "SingleStab" : "$(exp.params.method)"
        centers, avgs = avg_by_effort(exp.runs, bin_edges, xcol, ycol, y_agg, y_transform)
        isempty(centers) && continue
        lines!(ax, centers, avgs; label=label, linewidth=2)
    end

    if y_axis == :success_rate
        ylims!(ax, 0, 100)
    elseif y_axis == :hits
        max_targets = maximum(exp.n_targets for exp in values(experiments))
        ylims!(ax, 0, max_targets)
    end

    axislegend(ax; position=:rb)

    save(output, fig)
    display(fig)
    fig
end

experiments = read_setup("data/keijzer4_0.025_9.5")
delete!(experiments, "fbf651e8")
delete!(experiments, "5005c278")
experiments2 = read_setup("data/keijzer4_noncircular_0.025_9.5")
experiments2["c4751128"].params.method = Stab
experiments["bla"] = experiments2["c4751128"]

plot_experiments(experiments, log_scale=false, x_axis=:effort, y_axis=:success_rate,
                 output="doc/figs/experiments_success.pdf")

experiments = read_setup("data/keijzer4_circular_0.025_9.5")
plot_experiments(experiments, log_scale=false, x_axis=:effort, y_axis=:success_rate,
                 output="doc/figs/experiments_success.pdf")

plot_experiments(experiments, log_scale=false, x_axis=:effort, y_axis=:hits,
                 output="doc/figs/experiments_hits.pdf")
plot_experiments(experiments, log_scale=false, x_axis=:time, y_axis=:success_rate,
                 output="doc/figs/experiments_success_time.pdf")
