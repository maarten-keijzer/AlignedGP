using AlignedGP
using Plots

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
    rows = Vector{NTuple{3,Float64}}()
    for line in eachline(path)
        if startswith(line, "#")
            k, v = split(lstrip(line, '#'), "="; limit=2)
            comments[strip(k)] = strip(v)
        elseif startswith(line, "time_running")
            continue  # header
        else
            parts = split(line, ",")
            length(parts) == 3 || continue
            push!(rows, (parse(Float64, parts[1]),
                         parse(Float64, parts[2]),
                         parse(Float64, parts[3])))
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

        data = Matrix{Float64}(undef, length(rows), 3)
        for (i, (t, e, h)) in enumerate(rows)
            data[i, :] = [t, e, h]
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

function avg_hits_by_effort(runs, bin_edges, xcol, y_transform)
    nbins = length(bin_edges) - 1
    avg = zeros(Float64, nbins)
    counts = zeros(Int, nbins)

    for data in runs
        effs = data[:, xcol]
        hits = data[:, 3]
        last_hits = NaN
        for b in 1:nbins
            lo, hi = bin_edges[b], bin_edges[b+1]
            mask = (effs .>= lo) .& (effs .< hi)
            if any(mask)
                last_hits = maximum(hits[mask])
            end
            isnan(last_hits) && continue  # run hasn't started yet
            avg[b] += y_transform(last_hits)
            counts[b] += 1
        end
    end

    bin_centers = [(bin_edges[b] + bin_edges[b+1]) / 2 for b in 1:nbins]
    valid = counts .> 0
    bin_centers[valid], avg[valid] ./ counts[valid]
end

function plot_experiments(experiments; nbins=50, log_scale=true, x_axis=:effort, y_axis=:hits)
    xcol = x_axis == :time ? 1 : 2
    xlabel = x_axis == :time ? "Time (s)" : "Effort"
    ylabel = y_axis == :success_rate ? "Success rate" : "Average max hits"

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

    plt = plot(
        xlabel=xlabel, ylabel=ylabel,
        xscale=log_scale ? :log10 : :identity,
        legend=:topleft, title="Experiment comparison"
    )

    for (hash_key, exp) in experiments
        y_transform = y_axis == :success_rate ? (h -> h >= exp.n_targets ? 1.0 : 0.0) : identity
        label = "$(exp.params.method)"
        centers, avgs = avg_hits_by_effort(exp.runs, bin_edges, xcol, y_transform)
        isempty(centers) && continue
        plot!(plt, centers, avgs; label=label, lw=2)
    end

    plt
end

experiments = read_setup("data/keijzer4_dup_9")
display(plot_experiments(experiments, log_scale=false, x_axis=:effort, y_axis=:success_rate))
