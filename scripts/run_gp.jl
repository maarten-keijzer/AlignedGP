using AlignedGP
using Plots

function print_report(strata::Vector{Vector{Tree}}, effort)
    nIndies = 0
    nHits = 0
    besthits = 0
    best = nothing        # min-loss individual (driver reporting routes through loss)
    bestloss = Inf
    for stratum in strata
        for indy in stratum
            nIndies += 1
            nHits += sum(indy.hits)
            besthits = max(besthits, sum(indy.hits))
            if indy.loss < bestloss
                bestloss = indy.loss
                best = indy
            end
        end
    end
    ntargets = length(first(first(strata)).hits)
    besthits = sum(best.hits)   # primary hits of the reported (min-loss) individual
    println("best $besthits/$ntargets (loss $(round(bestloss, digits=4)), secondary $(sum(best.secondary_hits))/$ntargets), avg $(round(nHits/nIndies, digits=2)), mse $(round(best.mse, digits=7)) complexity $(best.complexity) avg-pathlen $(round(best.pathlen_complexity / best.complexity, digits=3)), effort $(round(log10(effort), digits=4))")
    return besthits
end

function get_best(strata)
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
    ev = evaluate(bestindy, setup.inputs)
    pop, bestindy, distr, ev
end

function plot_best(strata)
   
    pop, bestindy, distr, ev = get_best(strata)
    lo = getproperty.(setup.interval_targets.intervals, :lo)
    hi = getproperty.(setup.interval_targets.intervals, :hi)
    scatter(setup.inputs[1], setup.noisy_targets, label="target", color=:black)
    plot!(setup.inputs[1], setup.ideal_targets, label="ideal", color=:black)
    plot!(setup.inputs[1], lo, fillrange=hi, fillalpha=0.2, alpha=0, label="interval")
    plt = plot!(setup.inputs[1], ev, label="y")
    display(plt)
    pop, bestindy, distr, ev
end

function scatter_best(strata)
   
    pop, bestindy, distr, ev = get_best(strata)
    @show compute_hits(ev, setup.interval_targets), sum(bestindy.hits)
    lo = getproperty.(setup.interval_targets.intervals, :lo)
    hi = getproperty.(setup.interval_targets.intervals, :hi)
    noisy = setup.noisy_targets
    order = sortperm(noisy)
    xs = noisy[order]
    plot(xs, lo[order], fillrange=hi[order], fillalpha=0.2, alpha=0, label="interval")
    plot!(xs, xs, color=:black, label="y=x")
    plt = scatter!(ev, noisy, label="model")
    display(plt)
    pop, bestindy, distr, ev
end


# Run inside a function so the loop bookkeeping (lasttime, lasteffort, …) are locals;
# at top level a `for` body is soft scope and reassigning these globals is ambiguous.
# Returns `true` when an individual solved every case in the inner band (feasible ⇒
# ratchet), `false` when the global effort cap was hit first (true stop).
function evolve!(strata, setup, effort; effortcap = 10, stoponmaxhits=true)
    lasttime = time()
    lasteffort = 0.0
    last_processed = 0
    for nindies in 1:typemax(Int)
        hits = iteratestrata!(strata, setup, effort)
        reportandbreak = stoponmaxhits && hits == length(setup.interval_targets)
        eff = AlignedGP.compute_effort(effort, length(setup.interval_targets))

        log10(eff) > effortcap && return false

        if time() - lasttime > 2 || reportandbreak
            lasttime = time()
            print("Δeffort $(round(log10(eff - lasteffort), digits=4)), Δprocessed $(nindies - last_processed) ")
            besthits = print_report(strata, eff)
            # plot_best(strata)
            lasteffort = eff
            last_processed = nindies
            if reportandbreak
                return true
            end
        end
    end
    return false
end

using AlignedGP.ReverseIntervals
function set_tol!(strata, setup, tol)

    noisy = setup.noisy_targets
    newintervals = IntervalVector(intervaltype.(noisy .- tol, noisy .+ tol))
    setup.interval_targets.intervals .= newintervals.intervals

    # Re-target from each tree's cached raw evals — no node re-evaluation. Recomputes
    # primary/secondary hits + loss against the new inner band and current outer floor,
    # reusing the tol-invariant slope/intercept/mse.
    for stratum in strata
        for i in eachindex(stratum)
            stratum[i] = AlignedGP.retarget(stratum[i], setup)
        end
    end
end


tol = 10.0
#tol=0.2
#setup = keijzer1(tol=tol, noise=0.0)
setup = keijzer4(tol=tol, noise = 0.01);
#setup = load_pmlb("1027_ESL", tol=tol)
#setup = load_pmlb("706_sleuth_case1202", tol=tol)
#setup = load_pmlb("1096_FacultySalaries", tol=tol)
#setup = load_pmlb("560_bodyfat", tol=tol)
#setup = load_pmlb("601_fri_c1_250_5", tol=tol)

setup.params.method = RecursiveStab
#setup.params.max_complexity=300
# Two-band ratchet starts with the outer floor equal to the inner band (single-band
# search); each ratchet freezes the floor at the tol just solved. The problem
# constructor already set tau_outer == tol, so no two-band pressure until the first shrink.
strata, effort = initstrata(setup);

# Feasibility-gated ratchet: only tighten the inner band once an individual has solved
# every case in it. The effort cap inside evolve! is the true global stop.
while tol > 1e-8
    solved = evolve!(strata, setup, effort; effortcap=1000)
    if length(setup.inputs) == 1
        plot_best(strata);
    else
        scatter_best(strata)
    end
    solved || break                          # effort exhausted before feasible → stop

    setup.params.tau_outer = tol             # freeze outer floor at the solved tol
    tol *= 0.61803398875
    @show tol
    set_tol!(strata, setup, tol)             # tighten inner band + retarget all
end

pop, bestindy, distr, ev = plot_best(strata);

