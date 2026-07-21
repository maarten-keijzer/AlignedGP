using AlignedGP
using CairoMakie
using Statistics

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
    pop, bestindy, distr, model_evaluations(bestindy)
end

function plot_best(strata)

    pop, bestindy, distr, ev = get_best(strata)
    lo = getproperty.(setup.interval_targets.intervals, :lo)
    hi = getproperty.(setup.interval_targets.intervals, :hi)
    x  = setup.inputs[1]
    fig = Figure()
    ax  = Axis(fig[1, 1]; xlabel = "x", ylabel = "y")
    band!(ax, x, lo, hi; color = (:blue, 0.2), label = "interval")
    scatter!(ax, x, setup.noisy_targets; color = :black, label = "target")
    lines!(ax, x, setup.ideal_targets; color = :black, label = "ideal")
    lines!(ax, x, ev; label = "y")
    axislegend(ax)
    display(fig)
    pop, bestindy, distr, ev
end

function scatter_best(strata)

    pop, bestindy, distr, ev = get_best(strata)
    lo = getproperty.(setup.interval_targets.intervals, :lo)
    hi = getproperty.(setup.interval_targets.intervals, :hi)
    noisy = setup.noisy_targets
    order = sortperm(noisy)
    xs = noisy[order]
    fig = Figure()
    ax  = Axis(fig[1, 1]; xlabel = "target", ylabel = "model")
    band!(ax, xs, lo[order], hi[order]; color = (:blue, 0.2), label = "interval")
    lines!(ax, xs, xs; color = :black, label = "y=x")
    scatter!(ax, ev, noisy; label = "model")
    axislegend(ax)
    display(fig)
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

function set_eps_tol!(strata, setup)

    pop = [indy for stratum in strata for indy in stratum]
    t = setup.noisy_targets

    tol = Vector{Float64}(undef, length(t))
    scratch = Vector{Float64}(undef, length(pop))    
    for c in eachindex(t)
        for j in eachindex(pop)
            val = abs(pop[j].slope * pop[j].evals[c] + pop[j].intercept - t[c])
            #val = abs(pop[j].evals[c] - t[c])
            if isnan(val) 
                val = Inf 
            end
            scratch[j] = val
        end
        med = median(scratch)
        tol[c] = quantile(scratch, 0.4) #median(abs.(scratch .- med))
        #@show tol[c], med, sum(isfinite.(scratch))
        if !isfinite(tol[c]) 
            @warn "Infinite interval"
        end
    end

    # set bounds
    newbounds = intervaltype.(t .- tol, t .+ tol)
    setup.interval_targets.intervals .= newbounds
    setup.params.tau_outer = maximum(tol)

    for stratum in strata
        for i in eachindex(stratum)
            stratum[i] = AlignedGP.retarget(stratum[i], setup)
        end
    end
end


tol = 100.0
#tol=0.2
#setup = keijzer1(tol=tol, noise=0.0)
setup = keijzer4(tol=tol, noise = 0.0);
#setup = load_pmlb("1027_ESL", tol=tol)
#setup = load_pmlb("706_sleuth_case1202", tol=tol)
#setup = load_pmlb("1096_FacultySalaries", tol=tol)
#setup = load_pmlb("560_bodyfat", tol=tol)
#setup = load_pmlb("601_fri_c1_250_5", tol=tol)
#setup = load_pmlb("583_fri_c1_1000_50", tol=tol)
#setup = pagie(tol=tol, noise = 0.0, unaries=[-,inv], binaries=[+,*])
#setup = feynman("I.11.19", tol=tol)
#setup = feynman("I.18.12", tol=tol)
#setup = feynman("I.13.4", tol=tol)
#setup = feynman("I.34.14")
#setup = feynman("III.19.51", tol=tol)
#setup = feynman("II.35.18", tol=tol)
#setup = feynman("III.7.38", tol=tol)
setup = feynman("I.15.3t", tol=tol)



setup.params.method = Stab
#setup.params.use_residual_lexicase = true

#setup.params.max_complexity=300
# Two-band ratchet starts with the outer floor equal to the inner band (single-band
# search); each ratchet freezes the floor at the tol just solved. The problem
# constructor already set tau_outer == tol, so no two-band pressure until the first shrink.
print("Initializing... ")
strata, effort = initstrata(setup);
println("done")
#set_eps_tol!(strata, setup)          # tighten inner band + retarget all
#plot_best(strata);

# Feasibility-gated ratchet: only tighten the inner band once an individual has solved
# every case in it. The effort cap inside evolve! is the true global stop.
function anneal(tol, strata, setup, effort;effortcap=1000)
    while tol > 1e-8
        solved = evolve!(strata, setup, effort; effortcap=effortcap)
        solved || break                          # effort exhausted before feasible → stop

        if length(setup.inputs) == 1
            plot_best(strata);
        else
            @show "scattering"
            scatter_best(strata)
        end

        setup.params.tau_outer = tol             # freeze outer floor at the solved tol
        tol *= 0.5
        @show tol
        setup.params.tau_outer = tol * 1.5
        set_tol!(strata, setup, tol)   
        #set_eps_tol!(strata, setup)          # tighten inner band + retarget all
    end
end

anneal(tol, strata, setup, effort)