
mutable struct ProblemSetup 
    inputs::Vector{Vector{Float64}}
    
    ideal_targets::Vector{Float64}
    noisy_targets::Vector{Float64}
    interval_targets::Vector{Interval{Float64, Closed, Closed}}

    symboltable::SymbolTable
    cross_mut_prob::Float64
    max_complexity::Int
    population_size::Int 

    max_lexicase_comparisons::Int

    rng::AbstractRNG
end

function simple_regression()
    x = sort(rand(100))
    t = @. 0.3 * x * sin(2π * x)
    tol = 0.01

    noisy = @. t +  tol/2 * randn()

    ProblemSetup(
        [x], 
        t, 
        noisy, 
        Interval.(noisy .- tol, noisy .+ tol),
        SymbolTable(1, [log, sqrt, exp], [+, *, -, /]),
        0.8,
        200, 
        1500,
        100,
        Random.GLOBAL_RNG
    )
end

