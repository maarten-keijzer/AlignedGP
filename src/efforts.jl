
export compute_effort 

mutable struct EffortStats 
    sum_evals::Int
    sum_stabs::Int
end

compute_effort(effort::EffortStats, n) = effort.sum_evals * n + effort.sum_stabs * n * log(n)

