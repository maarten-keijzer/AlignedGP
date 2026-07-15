using AlignedGP
using Random

with_rng(setup::ProblemSetup, rng) = ProblemSetup(
    setup.inputs, setup.ideal_targets, setup.noisy_targets, setup.interval_targets,
    setup.symboltable, setup.params, rng
)

with_params(setup::ProblemSetup, params) = ProblemSetup(
    setup.inputs, setup.ideal_targets, setup.noisy_targets, setup.interval_targets,
    setup.symboltable, params, setup.rng
)

function params_hash(setup)
    p = setup.params
    h = hash((p.method, p.population_size, p.max_complexity, p.cross_mut_prob,
              p.max_lexicase_comparisons, p.use_l2_scaling, p.constant_stab_probability,
              p.use_tournament_stratum, length(setup.ideal_targets)))
    string(h, base=16)[1:8]
end

function next_run_index(dir, experiment_id)
    existing = filter(readdir(dir)) do f
        !isnothing(match(Regex("^$(experiment_id)_(\\d+)\\.csv\$"), f))
    end
    isempty(existing) && return 1
    nums = [parse(Int, match(r"_(\d+)\.csv$", f)[1]) for f in existing]
    maximum(nums) + 1
end

function add_experiments!(futures, setup, dir; nruns, maxeffort, master_seed=42)
    mkpath(dir)
    # Snapshot params: GPParams is mutable and shared by reference across all
    # spawned tasks. The caller mutates setup.params.method between batches, and
    # @spawn only queues the tasks — they read p.method later, after the mutation.
    # Copy once here (on the calling thread, before any task runs) so this batch
    # sees a stable, private params.
    setup = with_params(setup, deepcopy(setup.params))
    experiment_id = params_hash(setup)
    start_index = next_run_index(dir, experiment_id)
    for i in start_index:(start_index + nruns - 1)
        future = Threads.@spawn begin
            filename = joinpath(dir, "$(experiment_id)_$(lpad(i, 3, '0')).csv")
            run_setup = with_rng(setup, Xoshiro(master_seed + i))
            println("Starting run $filename")
            open(filename, "w") do io
                p = run_setup.params
                println(io, "# method=$(p.method)")
                println(io, "# population_size=$(p.population_size)")
                println(io, "# max_complexity=$(p.max_complexity)")
                println(io, "# cross_mut_prob=$(p.cross_mut_prob)")
                println(io, "# max_lexicase_comparisons=$(p.max_lexicase_comparisons)")
                println(io, "# use_l2_scaling=$(p.use_l2_scaling)")
                println(io, "# constant_stab_probability=$(p.constant_stab_probability)")
                println(io, "# use_tournament_stratum=$(p.use_tournament_stratum)")
                println(io, "# n_targets=$(length(run_setup.ideal_targets))")
                println(io, "time_running,eff,maxhits,iterations,best_mse")

                strata, effort = initstrata(run_setup)
                start_time = time()
                eff = AlignedGP.compute_effort(effort, length(run_setup.interval_targets))
                last_effort = eff
                iterations = 0
                while log10(eff) < maxeffort
                    iteratestrata!(strata, run_setup, effort)
                    iterations += 1
                    eff = AlignedGP.compute_effort(effort, length(run_setup.interval_targets))

                    if eff - last_effort > 1e7
                        maxhits = maximum(sum(indy.hits) for stratum in strata for indy in stratum)
                        time_running = time() - start_time
                        all_mses = (indy.mse for stratum in strata for indy in stratum if isfinite(indy.mse))
                        best_mse = isempty(all_mses) ? NaN : minimum(all_mses)

                        println(io, "$time_running,$eff,$maxhits,$iterations,$best_mse")
                        flush(io)

                        if maxhits == length(run_setup.ideal_targets)
                            break
                        end
                        last_effort = eff
                    end
                end
            end
        end
        push!(futures, future)
    end
end

tol=0.025
effort = 9.5
nruns = 20
setup = keijzer4(tol=tol)
dir = "data/keijzer4_$(tol)_$(effort)"

futures = []
setup.params.method = RecursiveStab
add_experiments!(futures, setup, dir, nruns=nruns, maxeffort=effort)
setup.params.method = Standard
add_experiments!(futures, setup, dir, nruns=nruns, maxeffort=effort)
setup.params.method = Stab
add_experiments!(futures, setup, dir, nruns=nruns, maxeffort=effort)

# wait until everything is done
[wait(future) for future in futures]

