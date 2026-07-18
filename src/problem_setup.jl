@enum OptMethod begin 
    Standard
    Stab
    RecursiveStab
    ConstantStab
end

export OptMethod, Standard, Stab, RecursiveStab, ConstantStab

struct SymbolTable
    nvars::Int
    unaries::Vector{Function}
    binaries::Vector{Function}
end

@kwdef mutable struct GPParams 
    method::OptMethod = RecursiveStab
    cross_mut_prob::Float64 = 0.8
    max_complexity::Int = 150
    population_size::Int = 5_000
    max_lexicase_comparisons::Int = 100
    use_l2_scaling::Bool = true
    constant_stab_probability = 0.2
    use_tournament_stratum = false
    # Outer (secondary) tolerance band half-width for the two-band ratchet.
    # Start it equal to the inner tol; the feature is OFF while tau_outer == tol
    # (secondary band ≡ inner band ⇒ single-band search).
    tau_outer::Float64 = 0.0
    # Experimental: swap PARENT selection from two-band lexicase to residual ε-lexicase
    # (semi-dynamic, MAD-based). OFF by default ⇒ byte-identical to the two-band path.
    # Constant fitting, replacement, and hit/loss bookkeeping are unaffected.
    use_residual_lexicase::Bool = false
end

struct ProblemSetup 
    inputs::Vector{Vector{Float64}}
    
    ideal_targets::Vector{Float64}
    noisy_targets::Vector{Float64}
    interval_targets::IntervalVector

    symboltable::SymbolTable

    params::GPParams

    rng::AbstractRNG
end

export keijzer1, keijzer4, keijzer4_dup, load_pmlb

function keijzer4(;noise = 0.0, tol=0.01, unaries=[sin, exp, log, sqrt], binaries=[+, *, -, /])
    x = collect(0.0:0.1:10.0)
    t = @. x^3 * exp(-x) * cos(x) * sin(x) * (sin(x)^2 * cos(x) - 1)

    noisy = @. t +  noise * randn()

    ProblemSetup(
        [x],
        t,
        noisy,
        IntervalVector(intervaltype.(noisy .- tol, noisy .+ tol)),
        SymbolTable(1, unaries, binaries),
        GPParams(tau_outer = tol),
        Random.GLOBAL_RNG
    )
end

function keijzer4_dup()
    x = collect(0.0:0.1:10.0)
    t = @. x^3 * exp(-x) * cos(x) * sin(x) * (sin(x)^2 * cos(x) - 1)
    tol1 = 0.01
    tol2 = 0.025
    int1 = IntervalVector(intervaltype.(t .- tol1, t .+ tol1))
    int2 = IntervalVector(intervaltype.(t .- tol2, t .+ tol2))

    ProblemSetup(
        [vcat(x,x)],
        vcat(t,t),
        vcat(t,t),
        vcat(int1, int2),
        SymbolTable(1, [exp, cos, sin], [+, *, -, /]),
        GPParams(tau_outer = tol2),
        Random.GLOBAL_RNG
    )
end


# Requires JULIA_PYTHONCALL_EXE to point at a Python with pmlb installed,
# e.g. export JULIA_PYTHONCALL_EXE=/opt/homebrew/bin/python3.11
function load_pmlb(
    dataset::AbstractString;
    unaries=[sqrt, log, exp],
    binaries=[+, -, *, /],
    tol::Float64=0.01,
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    pmlb = pyimport("pmlb")
    df = pmlb.fetch_data(dataset)
    cols = pyconvert(Vector{String}, df.columns.tolist())
    data = pyconvert(Matrix{Float64}, df.values)

    target_col = findfirst(==("target"), cols)
    feature_cols = setdiff(1:length(cols), [target_col])

    targets = data[:, target_col]
    inputs = [data[:, i] for i in feature_cols]

    params = GPParams()
    params.tau_outer = tol
    params.max_lexicase_comparisons = min(params.max_lexicase_comparisons, length(targets))

    ProblemSetup(
        inputs,
        targets,
        targets,
        IntervalVector(intervaltype.(targets .- tol, targets .+ tol)),
        SymbolTable(length(inputs), unaries, binaries),
        params,
        rng
    )
end

function keijzer1(; noise = 0.0, tol=0.01, unaries=[sqrt, log, exp], binaries=[+,-,*,/], rng::AbstractRNG=Random.GLOBAL_RNG)
    x = sort(rand(rng, 100))
    t = @. 0.3 * x * sin(2π * x)

    noisy = @. t +  noise * randn()

    ProblemSetup(
        [x],
        t,
        noisy,
        IntervalVector(intervaltype.(noisy .- tol, noisy .+ tol)),
        SymbolTable(1, unaries, binaries),
        GPParams(tau_outer = tol),
        rng
    )
end

export pagie
function pagie(;noise = 0.0, tol=0.01, unaries=[sin, exp, log, sqrt], binaries=[+, *, -, /])
    x = rand(500) * 10 .- 5
    y = rand(500) * 10 .- 5
    t = @. 1/(1+x^-4) + 1/(1+y^-4)

    noisy = @. t +  noise * randn()

    ProblemSetup(
        [x, y],
        t,
        noisy,
        IntervalVector(intervaltype.(noisy .- tol, noisy .+ tol)),
        SymbolTable(1, unaries, binaries),
        GPParams(tau_outer = tol),
        Random.GLOBAL_RNG
    )
end

