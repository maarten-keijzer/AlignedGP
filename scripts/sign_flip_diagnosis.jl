using AlignedGP

# Demonstrates the sign-flip / parent-regression scenario identified during diagnosis.
#
# Setup: parent = left * right, two data points, targets [1,3] for both.
# Left structural evals: [7.6, -0.4]  (what you get after left is optimized with c_left=3.6)
# Right structural evals: [0.1, -5.0]  (point-1 is near zero, vulnerable to sign flip)

left_evals      = [7.6, -0.4]
right_structural = [0.1, -5.0]
targets = [CIntervals(CInterval(1.0, 3.0)), CIntervals(CInterval(1.0, 3.0))]

println("=" ^ 60)
println("STEP 1 — right surrogates given current left_evals")
println("=" ^ 60)
right_surrogates = [rightinverse(*, targets[j], left_evals[j]) for j in eachindex(targets)]
println("right_surrogates : ", right_surrogates)

shifted_right = [right_surrogates[j] .- right_structural[j] for j in eachindex(right_surrogates)]
println("valid c_right    : ", shifted_right)
r_right, d_right = max_overlap_region(shifted_right)
println("optimal c_right  : depth=", d_right, "  region=", r_right)
c_right = select_constant(r_right)
println("chosen c_right   : ", c_right)

right_evals_optimized = right_structural .+ c_right
println("right_evals after opt : ", right_evals_optimized)
println("  → sign of right_evals before: ", sign.(right_structural))
println("  → sign of right_evals after : ", sign.(right_evals_optimized))
sign_flips = count(i -> sign(right_structural[i]) != sign(right_evals_optimized[i]), eachindex(right_structural))
println("  → sign flips: ", sign_flips)

println()
println("=" ^ 60)
println("STEP 2 — parent hits BEFORE right optimization")
println("=" ^ 60)
parent_before = left_evals .* right_structural
println("parent_evals     : ", parent_before)
r_p, d_p = max_overlap_region(targets .- parent_before)
println("parent hits (best c): ", d_p, "  region=", r_p)

println()
println("=" ^ 60)
println("STEP 3 — parent hits AFTER right optimization")
println("=" ^ 60)
parent_after = left_evals .* right_evals_optimized
println("parent_evals     : ", parent_after)
r_p2, d_p2 = max_overlap_region(targets .- parent_after)
println("parent hits (best c): ", d_p2, "  region=", r_p2)

if d_p2 < d_p
    println("\n*** PARENT REGRESSION: ", d_p, " → ", d_p2, " hits ***")
elseif d_p2 == d_p
    println("\nParent hits unchanged: ", d_p)
else
    println("\nParent hits improved: ", d_p, " → ", d_p2)
end

println()
println("=" ^ 60)
println("STEP 4 — what c_right WOULD cause a sign flip on point-1")
println("=" ^ 60)
c_flip = -right_structural[1] - 0.05   # just past the zero crossing
right_flipped = right_structural .+ c_flip
println("c_right that flips point-1: ", c_flip)
println("right_evals with that c   : ", right_flipped)
parent_flipped = left_evals .* right_flipped
println("parent_evals              : ", parent_flipped)
r_pf, d_pf = max_overlap_region(targets .- parent_flipped)
println("parent hits (best c)      : ", d_pf, "  region=", r_pf)
println("Is c_flip in optimal region? ", c_flip in r_right)
