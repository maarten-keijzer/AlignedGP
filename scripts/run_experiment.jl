using AlignedGP
using Random

with_rng(setup::ProblemSetup, rng) = ProblemSetup(
    setup.inputs, setup.ideal_targets, setup.noisy_targets, setup.interval_targets,
    setup.symboltable, setup.params, rng
)

function params_hash(setup)
    p = setup.params
    h = hash((p.method, p.population_size, p.max_complexity, p.cross_mut_prob,
              p.max_lexicase_comparisons, p.use_l2_scaling, p.constant_stab_probability,
              p.use_tournament_stratum, length(setup.ideal_targets)))
    string(h, base=16)[1:8]
end

function run_experiments(setup, dir; nruns, maxeffort, master_seed=42)
    mkpath(dir)
    experiment_id = params_hash(setup)
    Threads.@threads for i in 1:nruns
        filename = joinpath(dir, "$(experiment_id)_$(lpad(i, 3, '0')).csv")
        run_setup = with_rng(setup, Xoshiro(master_seed + i))
        println("Starting run $filename")
        open(filename, "w") do io
            p = run_setup.params
            println(io, "# method=$(p.method)")
            println(io, "# population_size=$(p.population_size)")
            println(io, "# max_complexity=$(p.max_complexity)")
            println(io, "# cross_mut_prob=$(p.cross_mut_prob)")
            println(io, "# n_targets=$(length(run_setup.ideal_targets))")
            println(io, "time_running,eff,maxhits")

            strata, effort = initstrata(run_setup)
            start_time = time()
            eff = AlignedGP.compute_effort(effort, length(run_setup.interval_targets))
            last_effort = eff
            while log10(eff) < maxeffort
                iteratestrata!(strata, run_setup, effort)
                eff = AlignedGP.compute_effort(effort, length(run_setup.interval_targets))

                if eff - last_effort > 1e7 
                    maxhits = maximum(sum(indy.hits) for stratum in strata for indy in stratum)
                    time_running = time() - start_time

                    println(io, "$time_running,$eff,$maxhits")
                    flush(io)

                    if maxhits == length(run_setup.ideal_targets)
                        break
                    end
                    last_effort = eff
                end
            end
        end
    end
end

setup = keijzer4(tol=0.025)
dir = "data/keijzer4_0.025"

run_experiments(setup, dir, nruns=10, maxeffort=8.1)
setup.params.method = Standard
run_experiments(setup, dir, nruns=10, maxeffort=8.1)

