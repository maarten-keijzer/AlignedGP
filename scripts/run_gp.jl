using AlignedGP
using Plots

function plot_best(strata)
    pop = [indy for stratum in strata for indy in stratum]

    besthits = 0 
    bestindy = first(pop)
    totalhits = copy(bestindy.hits)
    distr = zeros(Int, length(totalhits))
    for indy in pop 
        if sum(indy.hits) > besthits || (sum(indy.hits) == besthits && indy.complexity < bestindy.complexity)
            besthits = sum(indy.hits)
            bestindy = indy        
        end
        distr .+= indy.hits
    end

    #bar(distr)

    ev = evaluate(bestindy, setup.inputs)

    lo = first.(setup.interval_targets |> flatten)
    hi = last.(setup.interval_targets |> flatten)
    plot(setup.inputs[1], setup.ideal_targets, label="ideal", color=:black)
    plot!(setup.inputs[1], lo, fillrange=hi, fillalpha=0.2, alpha=0, label="interval")
    plt = plot!(setup.inputs[1], ev)
    display(plt)
    pop, bestindy, distr
end

setup = AlignedGP.keijzer1();
method = RecursiveStab
strata, effort = initstrata(setup);

begin
    lasttime=time()
    lasteffort = 0; last_processed = 0
    @time for nindies in 1:typemax(Int)
        iteratestrata!(strata, setup, effort, method=method)
        eff = AlignedGP.compute_effort(effort, length(setup.interval_targets))
        
        if log10(eff) > 8.6
            break 
        end

        if time() - lasttime > 2
            lasttime=time() 
            print("Δeffort $(round(log10(eff - lasteffort), digits=4)), Δprocessed $(nindies-last_processed) ")
            besthits = print_report(strata, eff)
            # plot_best(strata)
            if besthits==length(setup.interval_targets) break end
            lasteffort = eff
            last_processed = nindies
        end
    end
end 

pop, bestindy, distr = plot_best(strata)
bar(distr)

using Statistics
mean(abs2, setup.ideal_targets .- ev)

#bar([length(stratum) for stratum in strata])

complexities = complexity.(pop)
mse = [sum(indy.mse) for indy in pop]
scatter(complexities, mse)
