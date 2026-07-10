using AlignedGP

function dorun(setup, logeffort)
    strata, effort = initstrata(setup);

    starttime=time()
    while(true)
        hits = iteratestrata!(strata, setup, effort)
        # if time() - starttime > seconds
        #     break;
        # end

        eff = compute_effort(effort, length(setup.interval_targets))
        if log10(eff) > logeffort 
            break
        end

        if hits == length(setup.interval_targets)
            break
        end
    end

    maximum(sum(indy.hits) for stratum in strata for indy in stratum)
end

function compare_two(setup1, setup2, effort, ntries)

    wins1 = 0
    wins2 = 0 
    ties = 0

    sumhits1 = 0
    sumhits2 = 0 

    n = 0

    for i = 1:ntries
        task1 = Threads.@spawn dorun(setup1, effort)
        task2 = Threads.@spawn dorun(setup2, effort)
        hits1, hits2 = fetch.( (task1, task2) )
        if hits1 > hits2
            wins1 += 1
        elseif hits2 > hits1
            wins2 += 1
        else
            ties += 1
        end

        sumhits1 += hits1
        sumhits2 += hits2
        n += 1

        println("$wins1 $wins2 $ties $(wins1/(wins1 + wins2))  $(sumhits1/n)  $(sumhits2/n) $hits1 $hits2")
    end

end

setup1 = keijzer4(tol=0.025)#([sqrt, exp, log], [+,*,/,-])
setup2 = deepcopy(setup1)

setup1.params.method = RecursiveStab
setup1.params.use_tournament_stratum = true
setup2.params.method = Standard
setup2.params.use_tournament_stratum = true

compare_two(setup1, setup2, 8.6, 200)


