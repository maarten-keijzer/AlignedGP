using AlignedGP
using Plots
using Statistics: median, quantile
using IntervalArithmetic: sup, inf

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

    bestindy = first(pop)
    bestloss = bestindy.loss
    distr = zeros(Int, length(bestindy.hits))
    for indy in pop
        if indy.loss < bestloss   # strict < keeps the first (smallest) min-loss individual
            bestloss = indy.loss
            bestindy = indy
        end
        distr .+= indy.hits
    end
    pop, bestindy, distr, model_evaluations(bestindy)
end

function plot_best(strata, plots)
   
    pop, bestindy, distr, ev = get_best(strata)
    lo = getproperty.(setup.interval_targets.intervals, :lo)
    hi = getproperty.(setup.interval_targets.intervals, :hi)
    scatter(setup.inputs[1], setup.noisy_targets, label="target", color=:black)
    plot!(setup.inputs[1], setup.ideal_targets, label="ideal", color=:black)
    plot!(setup.inputs[1], lo, fillrange=hi, fillalpha=0.2, alpha=0, label="interval")
    plt = plot!(setup.inputs[1], ev, label="y", lw=2)
    push!(plots, plt)
    display(plt)
    pop, bestindy, distr, ev
end

function scatter_best(strata, plots)
   
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
    push!(plots, plt)
    display(plt)
    pop, bestindy, distr, ev
end


# Run inside a function so the loop bookkeeping (lasttime, lasteffort, …) are locals;
# at top level a `for` body is soft scope and reassigning these globals is ambiguous.
# Returns `true` when an individual solved every case in the inner band (feasible ⇒
# ratchet), `false` when the global effort cap was hit first (true stop).
function evolve!(strata, setup, effort; effortcap = 10, stoponmaxhits=true, plots=[])
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
            
            if !isempty(plots)
                push!(plots, plots[end]) # gives it some indication of speed in animation
            end

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

function set_tol_ratchet!(strata, setup, tol_vec; shrink=0.9, noise=0.0,
                          floor_abs=1e-8, floor_rel=1e-6)
    noisy = setup.noisy_targets
    intervals = setup.interval_targets.intervals
    n = length(noisy)

    setup.params.tau_outer = maximum(tol_vec)   # freeze outer at just-solved widths

    _, _, _, champ = get_best(strata)           # champion's chosen (raw|scaled) eval series

    
    # find the case furthest to the boundary
    dist = @. min(sup(intervals) - champ, champ - inf.(intervals))
    
    # find the current tolerance vector: intervalvec = noise +/- tol_vec
    tol_vec .-= 1.01 .* dist 

    newintervals = IntervalVector(intervaltype.(noisy .- tol_vec, noisy .+ tol_vec))
    setup.interval_targets.intervals .= newintervals.intervals

    for stratum in strata
        for i in eachindex(stratum)
            stratum[i] = AlignedGP.retarget(stratum[i], setup)
        end
    end
    return tol_vec
end

# Champion-anchored sibling of `anneal`: thread a per-case `tol_vec`. Each ratchet is
# guaranteed to invalidate the incumbent (nbite ≥ 1) or the run has converged to the noise
# floor (nbite == 0 ⇒ stop). `shrink` controls how hard each case is bitten below the
# champion's residual; `noise` is the stop floor.
function anneal_ratchet(strata, setup, effort; effortcap=10, tol=10, shrink=0.9,
                        noise=0.0, floor_abs=1e-8, floor_rel=1e-6, plots=[])
    n = length(setup.noisy_targets)
    tol_vec = fill(float(tol), n)
    _, best_solution, _ = get_best(strata)
    while maximum(tol_vec) > 1e-8
        solved = evolve!(strata, setup, effort; effortcap=effortcap, plots=plots)
        solved || break

        if length(setup.inputs) == 1
            plot_best(strata, plots)
        else
            scatter_best(strata, plots)
        end

        _, best_solution, _ = get_best(strata)

        tol_vec = set_tol_ratchet!(strata, setup, tol_vec; shrink=shrink,
                                             noise=noise, floor_abs=floor_abs, floor_rel=floor_rel)

    end
    best_solution
end


tol = 10.0
#tol=0.2
setup = keijzer1(tol=tol, noise=0.01)
#setup = keijzer4(tol=tol, noise = 0.02);
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

plots = []
#best = anneal(strata, setup, effort, effortcap=11.0, tol=tol, plots=plots)
best = anneal_ratchet(strata, setup, effort, effortcap=11.0, tol=tol, shrink=0.9, noise=0.01, plots=plots)

#animation = @animate for pl in plots; plot(pl); end
#gif(animation, fps=2)