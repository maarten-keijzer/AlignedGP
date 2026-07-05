@enum OptMethod begin 
    Standard
    Stab
    RecursiveStab
end

export OptMethod, Standard, Stab, RecursiveStab

struct SymbolTable
    nvars::Int
    unaries::Vector{Function}
    binaries::Vector{Function}
end

@kwdef mutable struct GPParams 
    cross_mut_prob::Float64 = 0.8
    max_complexity::Int = 100
    population_size::Int = 5_000

    max_lexicase_comparisons::Int = 100
end

struct ProblemSetup 
    inputs::Vector{Vector{Float64}}
    
    ideal_targets::Vector{Float64}
    noisy_targets::Vector{Float64}
    interval_targets::Vector{CIntervals}

    symboltable::SymbolTable

    params::GPParams

    rng::AbstractRNG
end

export keijzer1, keijzer4, keijzer4_dup

function keijzer4()
    x = collect(0.0:0.1:10.0)
    t = @. x^3 * exp(-x) * cos(x) * sin(x) * (sin(x)^2 * cos(x) - 1)
    tol = 0.01

    ProblemSetup(
        [x],
        t,
        t,
        CIntervals.(CInterval.(t .- tol, t .+ tol)),
        SymbolTable(1, [identity], [+, *, -, /]),
        GPParams(),
        Random.GLOBAL_RNG
    )
end

function keijzer4_dup()
    x = collect(0.0:0.1:10.0)
    t = @. x^3 * exp(-x) * cos(x) * sin(x) * (sin(x)^2 * cos(x) - 1)
    tol1 = 0.01
    tol2 = 0.025
    int1 = CIntervals.(CInterval.(t .- tol1, t .+ tol1))
    int2 = CIntervals.(CInterval.(t .- tol2, t .+ tol2))

    ProblemSetup(
        [vcat(x,x)],
        vcat(t,t),
        vcat(t,t),
        vcat(int1, int2),
        SymbolTable(1, [exp, cos, sin], [+, *, -, /]),
        GPParams(),
        Random.GLOBAL_RNG
    )
end


function keijzer1(unaries=[sqrt, log, exp], binaries=[+,-,*,/]; rng::AbstractRNG=Random.GLOBAL_RNG)
    x = sort(rand(rng, 100))
    t = @. 0.3 * x * sin(2π * x)
    tol = 0.01

    noisy = @. t #+  tol/2 * randn()

    ProblemSetup(
        [x],
        t,
        noisy,
        CIntervals.(CInterval.(noisy .- tol, noisy .+ tol)),
        SymbolTable(1, unaries, binaries),
        GPParams(),
        rng
    )
end

