# A LOWER CONVEX HULL over (complexity, error), both minimised. The front is a DENSE vector
# kept sorted by complexity ascending — which, being a Pareto hull, is simultaneously error
# descending. It stores ONLY the convex-hull vertices; interior/dominated points are
# discarded outright (not merely hidden from reporting). The hull is always a subset of the
# Pareto staircase. Complexity may be fractional, so lookups use binary search.
#
# Invariant: complexity strictly increasing, error strictly decreasing, and the (complexity,
# error) sequence is strictly convex from below — every stored vertex strictly changes the
# size/error slope (collinear middles are dropped).
#
# Why hull-only storage is exact. The front's point set only ever GROWS (each generation adds
# a batch; nothing underlying is removed). Adding points can only push the lower hull
# downward, so a point that is interior/above the hull today stays interior forever — a
# discarded point never needs resurrection. Thus keeping just the vertex set is identical to
# recomputing the hull over every point ever seen.
#
# Ties: equal-complexity points collapse to the lowest error; an exact (complexity, error)
# duplicate is dropped (the incumbent wins).
struct Front{T}
    front::Vector{T}
    complexityfn::Function
    errfn::Function
    Front{T}(; complexityfn=Base.length, errfn=getmse) where T = new(T[], complexityfn, errfn)
end

complexities_front(front::Front{T}) where T = front.complexityfn.(front.front)
errors_front(front::Front{T}) where T = front.errfn.(front.front)

# Cross product of (o→a) × (o→b) in (complexity, error) space. > 0 ⇒ b turns left of o→a
# (convex vertex from below, keep the middle); <= 0 ⇒ collinear or right turn (drop middle).
_cross(cx, cy, o, a, b) =
    (cx(a) - cx(o)) * (cy(b) - cy(o)) - (cy(a) - cy(o)) * (cx(b) - cx(o))

# Insert a single individual, maintaining the lower convex hull. Binary-search the complexity
# position, reject the individual if it lies on or above the current hull (interior /
# dominated / collinear), otherwise splice it in and pop the now-non-convex vertices on both
# sides. Returns `true` if the individual became a hull vertex, `false` if rejected.
function add_to_front!(frontstruct::Front{T}, indy::T) where T
    front = frontstruct.front
    cx = frontstruct.complexityfn
    cy = frontstruct.errfn
    cross(o, a, b) = _cross(cx, cy, o, a, b)

    c = cx(indy)
    e = cy(indy)

    lo = searchsortedfirst(front, indy; by = cx)     # first index with complexity >= c

    # Reject test: is indy on or above the current hull?
    if lo <= length(front) && cx(front[lo]) == c
        # Exact-complexity incumbent: indy must strictly beat its error to matter.
        cy(front[lo]) <= e && return false
        # Otherwise indy is strictly below that vertex ⇒ a vertex; the dominated incumbent
        # (error > e) is removed by the right-side dominance sweep below.
    elseif lo > 1 && lo <= length(front)
        # Between two vertices: reject if on/above the segment spanning c (<= drops collinear).
        cross(front[lo-1], indy, front[lo]) <= 0 && return false
    elseif lo > 1
        # New maximum complexity: reject if Pareto-dominated by the current last vertex.
        cy(front[lo-1]) <= e && return false
    # else new minimum complexity (or empty front): always a hull vertex.
    end

    # Right side: drop right neighbours dominated by indy (error >= e, a contiguous prefix
    # since error is descending), then convex-pop the rest against indy.
    r = lo
    while r <= length(front) && cy(front[r]) >= e
        r += 1
    end
    while r < length(front) && cross(indy, front[r], front[r+1]) <= 0
        r += 1
    end

    # Left side: all left neighbours have error > e (else indy would have been rejected), so
    # only convex-pop is needed.
    l = lo - 1
    while l >= 2 && cross(front[l-1], front[l], indy) <= 0
        l -= 1
    end

    # Splice: [l+1, r-1] is the contiguous run of dropped vertices; indy takes their place.
    deleteat!(front, (l+1):(r-1))
    insert!(front, l + 1, indy)
    return true
end

# Merge a whole batch into the hull in one shot: sort the batch by (complexity ascending,
# error ascending), two-pointer merge it against the existing hull into one complexity-
# ascending stream, then recompute the lower hull with a combined Pareto-guard + monotone-
# chain sweep. Because the hull is a unique, order-independent function of the point set
# (collinear dropped), this produces exactly the same front as adding the batch one-by-one
# with add_to_front!. Rebuilds `front` in place; returns the Front.
function merge_with_front!(frontstruct::Front{T}, batch::AbstractVector{T}) where T
    front = frontstruct.front
    cx = frontstruct.complexityfn
    cy = frontstruct.errfn
    cross(o, a, b) = _cross(cx, cy, o, a, b)

    # Sort the batch by complexity, breaking ties by error (both ascending).
    sorted = sort(batch; by = x -> (cx(x), cy(x)))

    # Two-pointer merge existing hull + sorted batch into one complexity-ascending stream.
    # On a full (complexity, error) tie the incumbent (old) sorts first so it wins the tie.
    old = front
    merged = Vector{T}(undef, length(old) + length(sorted))
    i = 1; j = 1; k = 1
    while i <= length(old) && j <= length(sorted)
        if (cx(old[i]), cy(old[i])) <= (cx(sorted[j]), cy(sorted[j]))
            merged[k] = old[i]; i += 1
        else
            merged[k] = sorted[j]; j += 1
        end
        k += 1
    end
    while i <= length(old); merged[k] = old[i]; i += 1; k += 1; end
    while j <= length(sorted); merged[k] = sorted[j]; j += 1; k += 1; end

    # Sweep the complexity-ascending stream, recomputing the lower hull:
    #  (1) Pareto guard: skip anything not strictly better than the best error so far (drops
    #      dominated points and collapses equal-complexity duplicates to the lowest error);
    #  (2) convex pop: drop the previous vertex while the turn is non-strictly-convex.
    result = T[]
    for cand in merged
        if !isempty(result) && cy(cand) >= cy(result[end])
            continue
        end
        while length(result) >= 2 && cross(result[end-1], result[end], cand) <= 0
            pop!(result)
        end
        push!(result, cand)
    end

    empty!(front)
    append!(front, result)
    return frontstruct
end
