using AlignedGP

function dorun(method, setup, seconds)
    strata, effort = initstrata(setup);

    starttime=time()
    while(true)
        hits = iteratestrata!(strata, setup, effort, method=method)
        if time() - starttime > seconds
            break;
        end

        if hits == length(setup.interval_targets)
            break
        end
    end

    maximum(sum(indy.hits) for stratum in strata for indy in stratum)
end

function compare_two(method1, method2, setup, seconds, ntries)

    wins1 = 0
    wins2 = 0 
    ties = 0

    sumhits1 = 0
    sumhits2 = 0 

    n = 0

    for i = 1:ntries
        task1 = Threads.@spawn dorun(method1, setup, seconds)
        task2 = Threads.@spawn dorun(method2, setup, seconds)
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

setup = keijzer4()#([sqrt, exp, log], [+,*,/,-])
compare_two(RecursiveStab, Standard, setup, 25, 200)
