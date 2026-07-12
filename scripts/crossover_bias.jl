using AlignedGP
using Plots

pathlen_complexity = AlignedGP.pathlen_complexity

function xo(tree1::Node, tree2::Node)
    insert(tree1, tree2[rand(1:length(tree2))], rand(1:length(tree1)))
end


setup = keijzer1([],[+,*,/,+])

s = Float64[]
ss = Float64[]
plot()

pop = [init(setup.symboltable, rand(5:100)) for i in 1:1000];
n = length(pop)
for i in 1:500
    for _ in eachindex(pop)
        pop[rand(1:length(pop))] = xo(pop[rand(1:length(pop))], pop[rand(1:length(pop))])
    end

    if i > 100 
        sumsize = sum(length, pop)
        sumpath = sum(pathlen_complexity.(pop) ./ length.(pop))
        @show i, sumsize/n, sumpath/n
        push!(s, sumsize/n)
        push!(ss, sumpath/n)
    end    
end

scatter(s, ss)

using Polynomials 

poly = fit(ss, s, 2)
sum(abs2, poly.(ss) .- s)

