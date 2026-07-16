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
    #ev = bestindy.slope .* ev .+ bestindy.intercept

    lo = getproperty.(setup.interval_targets.intervals, :lo)
    hi = getproperty.(setup.interval_targets.intervals, :hi)
    plot(setup.inputs[1], setup.ideal_targets, label="ideal", color=:black)
    plot!(setup.inputs[1], lo, fillrange=hi, fillalpha=0.2, alpha=0, label="interval")
    plt = plot!(setup.inputs[1], ev)
    display(plt)
    pop, bestindy, distr, ev
end

#setup = keijzer1([sqrt, sin], [+,*,/,-])
setup = AlignedGP.keijzer4(tol=0.05, unaries=[exp, log, sqrt, sin]);
#setup = load_pmlb("1027_ESL", tol=0.5)
#setup = load_pmlb("706_sleuth_case1202", tol=5.0)
# setup = load_pmlb("1096_FacultySalaries", tol=0.1, 
    # unaries=[sqrt])
#setup = load_pmlb("560_bodyfat", tol=0.05)

setup.params.method = RecursiveStab
strata, effort = initstrata(setup);

# Run inside a function so the loop bookkeeping (lasttime, lasteffort, …) are locals;
# at top level a `for` body is soft scope and reassigning these globals is ambiguous.
function evolve!(strata, setup, effort; effortcap = 10)
    lasttime = time()
    lasteffort = 0.0
    last_processed = 0
    @time for nindies in 1:typemax(Int)
        iteratestrata!(strata, setup, effort)
        eff = AlignedGP.compute_effort(effort, length(setup.interval_targets))

        log10(eff) > effortcap && break

        if time() - lasttime > 2
            lasttime = time()
            print("Δeffort $(round(log10(eff - lasteffort), digits=4)), Δprocessed $(nindies - last_processed) ")
            besthits = print_report(strata, eff)
            # plot_best(strata)
            besthits == length(setup.interval_targets) && break
            lasteffort = eff
            last_processed = nindies
        end
    end
end

evolve!(strata, setup, effort)

pop, bestindy, distr, ev = plot_best(strata);
bar(distr);

using Statistics
mean(abs2, setup.ideal_targets .- ev)

#bar([length(stratum) for stratum in strata])

complexities = complexity.(pop)
mse = [sum(indy.mse) for indy in pop]
scatter(complexities, mse)
