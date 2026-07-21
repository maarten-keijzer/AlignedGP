# Standard (semi-dynamic) ε-lexicase GP — a clean, self-contained baseline.
#
# Deliberately NON-INTRUSIVE: it reuses the problem setup (`ProblemSetup` via
# `load_pmlb`/`keijzer1`/…), `valid_init`, the standard variation operators
# (`standard_crossover`, `sizefair_mutation`), `linear_scale`, and `evaluate` — but it
# does NOT touch the interval machinery, the `Tree` fitness type, or the module's
# `lexicase`/`two_band_lexicase`. Selection is a from-scratch ε-lexicase on continuous
# residuals, with optional parsimony pressure by adding tree size as extra ε-lexicase
# case(s) (see SIZE_CASE_WEIGHT). No params or problem constructors are modified (params
# are read-only; the `tol` argument only feeds the untouched constructor and is ignored).

using AlignedGP
using Statistics: median
using Random: randperm

const POP_SIZE      = 5_000   # generational population size N
const NGEN          = 50    # fixed generation budget (stopping rule)
const INIT_MAX_RANGE = 5:30     # initial trees ramp target range
const MAX_RETRIES   = 20     # variation retries before falling back to a parent copy
const SIZE_TOURNEY  = 3

# Parsimony pressure via size-as-cases: tree size is appended to each individual's error
# vector as one or more extra "cases", so ε-lexicase filters on size with the exact same
# machinery it uses for residuals. Weighting is by DUPLICATION — the size-case appears
# `ncopies` times in every shuffled sweep, so more copies ⇒ earlier/stronger pressure.
#   :off     → 0 copies (byte-identical to the no-parsimony baseline)
#   :uniform → 1 copy   (size is one case among m+1)
#   :logn    → round(log(m)) copies (m = #cases): a minority, so residuals still dominate
const SIZE_CASE_WEIGHT = :uniform   # :off | :uniform | :logn

# Number of duplicate size-cases, given the case count m.
function size_case_copies(m::Int)
    SIZE_CASE_WEIGHT === :off     && return 0
    SIZE_CASE_WEIGHT === :uniform && return 1
    SIZE_CASE_WEIGHT === :logn    && return max(1, ceil(Int, log2(m)))
    error("unknown SIZE_CASE_WEIGHT = $SIZE_CASE_WEIGHT")
end

# A lightweight individual: the tree plus its per-case values and aggregate MSE. The first
# m entries of `errors` are the linearly-scaled absolute residuals; the trailing `ncopies`
# entries (if any) are all the tree's size — the duplicated size-cases. Nothing here mutates
# in place, so individuals may be shared across generations safely (parent-copy fallback).
mutable struct Indy
    root::Node
    errors::Vector{Float64}   # [ |slope·eval[c]+intercept−target[c]| (×m) ; size (×ncopies) ]
    mse::Float64
    pathlen::Int
    Indy(root, errors, mse) = new(root, errors, mse, pathlen_complexity(root))
end

Base.length(indy::Indy) = length(indy.root)
pathlen(indy::Indy) = indy.pathlen
getmse(indy::Indy) = indy.mse

# Build an individual from a raw tree, appending `ncopies` size-cases. Returns `nothing`
# when the tree is non-finite on any RESIDUAL case — the caller rejects it, keeping the
# whole population finite so `linear_scale` and the per-case MAD are always well-defined.
# (Size-cases are always finite and never cause rejection.)
function make_indy(root::Node, setup::ProblemSetup, ncopies::Int)
    evals = evaluate(root, setup.inputs)
    all(isfinite, evals) || return nothing
    slope, intercept, mse = setup.params.use_l2_scaling ? 
        linear_scale(evals, setup.noisy_targets) :
        linear_scale_l1(evals, setup.noisy_targets)
    errors = abs.(slope .* evals .+ intercept .- setup.noisy_targets)
    all(isfinite, errors) || return nothing
    ncopies > 0 && append!(errors, fill(Float64(length(root)), ncopies))
    return Indy(root, errors, mse)
end

# Per-case ε = raw median absolute deviation of that case's values over the WHOLE current
# population. Computed once per generation (the "semi-dynamic" flavor recomputes only the
# elite during the sweep, not ε). All values are finite by construction. This runs over
# every column uniformly, so the trailing size-cases automatically get ε = MAD(sizes).
function case_epsilons(pop::Vector{Indy}, ncases::Int)
    eps = Vector{Float64}(undef, ncases)
    scratch = Vector{Float64}(undef, length(pop))
    for c in 1:ncases
        @inbounds for j in eachindex(pop)
            scratch[j] = pop[j].errors[c]
        end
        med = median(scratch)
        eps[c] = median(abs.(scratch .- med))
    end
    return eps
end

# Semi-dynamic ε-lexicase parent selection: sweep every case in a fresh random order;
# per case the elite is the minimum error over the SURVIVING pool, and survivors are those
# within `elite + ε[case]`. Collapse to one survivor, else break ties at random.
function eps_lexicase(pop::Vector{Indy}, epsilons::Vector{Float64}, rng)
    remaining = collect(1:length(pop))
    for case in randperm(rng, length(epsilons))
        length(remaining) <= 1 && break
        elite = minimum(pop[j].errors[case] for j in remaining)
        threshold = elite + epsilons[case]
        remaining = [j for j in remaining if pop[j].errors[case] <= threshold]
    end
    return length(remaining) == 1 ? pop[remaining[1]] : pop[rand(rng, remaining)]
end

function eps_lexicase_front(pop::Vector{Indy}, epsilons::Vector{Float64}, rng)
    @assert SIZE_CASE_WEIGHT == :uniform # make sure last one is eps weight
    if SIZE_TOURNEY == 1 
        remaining = collect(1:length(pop))
    else
        rp = randperm(rng, length(pop))
        remaining = Int[]
        while length(rp) >= SIZE_TOURNEY 
            best = pop!(rp)
            for _ = 2:SIZE_TOURNEY
                contender = pop!(rp)
                if length(pop[contender]) < length(pop[best])   
                    best = contender
                end
            end
            push!(remaining, best)
        end
    end

    for case in randperm(rng, length(epsilons)-1)
        length(remaining) <= 1 && break
        error_elite = minimum(pop[j].errors[case] for j in remaining)
        error_threshold = error_elite + epsilons[case]
       
        remaining = [j for j in remaining if pop[j].errors[case] <= error_threshold]
    end
    
    return length(remaining) == 1 ? pop[remaining[1]] : pop[rand(rng, remaining)]
end

eps_selection = eps_lexicase_front

# Produce one child via standard crossover (prob `cross_mut_prob`) or size-fair mutation,
# both from ε-lexicase-selected parents. Reject non-finite children and retry; after
# MAX_RETRIES failures, fall back to an (already valid) selected parent.
function make_child(pop::Vector{Indy}, epsilons::Vector{Float64},
                    setup::ProblemSetup, effort, ncopies::Int, rng)
    for _ in 1:MAX_RETRIES
        p1 = eps_selection(pop, epsilons, rng)
        if rand(rng) < setup.params.cross_mut_prob
            p2 = eps_selection(pop, epsilons, rng)
            root = AlignedGP.standard_crossover(p1.root, p2.root, effort, rng)
        else
            root = AlignedGP.sizefair_mutation(p1.root, setup.symboltable, effort, rng)
        end
        child = make_indy(root, setup, ncopies)
        child === nothing || return child
    end
    @warn "Too many retries $MAX_RETRIES"
    return eps_selection(pop, epsilons, rng)   # give up: reuse a valid parent
end

# Ramped initialization: cycle target sizes over INIT_MAX_RANGE. `valid_init` already
# guarantees finite output, so `make_indy` never rejects here (guarded anyway).
function init_pop(setup::ProblemSetup, n::Int, ncopies::Int, rng)
    pop = Indy[]
    while length(pop) < n
        tsize = rand(INIT_MAX_RANGE)
        root, _ = valid_init(setup.symboltable, tsize, setup.inputs, rng)
        child = make_indy(root, setup, ncopies)
        child === nothing && continue
        push!(pop, child)
    end
    return pop
end

function report(gen::Int, best::Indy, pop::Vector{Indy})
    avgsize = sum(length(ind.root) for ind in pop) / length(pop)
    genbest = argmin(ind -> ind.mse, pop)

    println("gen $gen: best mse=$(round(best.mse, sigdigits=5)) " *
            "size=$(length(best.root)) " *
            "gen mse=$(round(genbest.mse, sigdigits=5)) size=$(length(genbest)) " *
            "avgsize=$(round(avgsize, digits=1))")
end

# Generational, no elitism: the whole population is replaced each generation; the best-ever
# individual is tracked for reporting but never reinjected.
function run_eps_lexicase(setup::ProblemSetup, front::Front{Indy}; n::Int=POP_SIZE, ngen::Int=NGEN)
    rng     = setup.rng
    effort  = AlignedGP.EffortStats(0, 0)   # dummy sink for the variation operators
    m       = length(setup.noisy_targets)
    ncopies = size_case_copies(m)
    ncases  = m + ncopies                   # residual cases + duplicated size-cases
    println("SIZE_CASE_WEIGHT=$SIZE_CASE_WEIGHT ⇒ $ncopies size-case(s) over $m residual cases")

    pop  = init_pop(setup, n, ncopies, rng)
    
    merge_with_front!(front, pop)

    best = argmin(ind -> ind.mse, pop)
    report(0, best, pop)



    for gen in 1:ngen
        epsilons = case_epsilons(pop, ncases)
        newpop = Vector{Indy}(undef, n)
        Threads.@threads for k in 1:n
            newpop[k] = make_child(pop, epsilons, setup, effort, ncopies, rng)
        end
        pop = newpop
        
        merge_with_front!(front, pop)

        genbest = argmin(ind -> ind.mse, pop)
        genbest.mse < best.mse && (best = genbest)
        report(gen, best, pop)

        plotfront(front)
    end
    return best
end

using CairoMakie 

function plotfront(front::Front)

    sizes  = complexities_front(front)
    errors = errors_front(front)

    fig = Figure()
    ax  = Axis(fig[1,1]; xlabel = "complexity", ylabel = "mse", xscale = log10)
    CairoMakie.scatter!(ax, sizes, errors; label = "run")
    axislegend(ax)
    display(fig)
    return fig
end

# --- driver (mirrors run_gp.jl's top-level menu) ---------------------------------------
#setup = keijzer1(noise=0.0)
setup = keijzer4(noise=0.0)
#setup = load_pmlb("1027_ESL", tol=tol)
#setup = load_pmlb("560_bodyfat")
#setup = load_pmlb("601_fri_c1_250_5", tol=tol)
#setup = pagie()
#setup = feynman("I.11.19")
#setup = feynman("I.18.12")
#setup = feynman("I.13.4")
#setup = feynman("I.34.14")
#setup = feynman("III.19.51")
#setup = feynman("II.35.18")
#setup = feynman("III.7.38")
setup = feynman("I.15.3t")

front = Front{Indy}(errfn=getmse)
best = run_eps_lexicase(setup, front)
