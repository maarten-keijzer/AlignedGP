"""
    max_overlap_region(intervals) -> (; region, depth)
"""

using IntervalArithmetic: in_interval

function max_overlap_region(intervals::Vector{IntervalType})
  
    isempty(intervals) && return (; region=IntervalType[], depth=0)

    # Encode: (coordinate, 0) = enter, (coordinate, 1) = leave
    events = Vector{Tuple{Float64,Int}}(undef, 2 * length(intervals))
    for (i, iv) in enumerate(intervals)
        events[2i-1] = (iv.lo, 0)
        events[2i]   = (iv.hi, 1)
    end
    sort!(events)

    depth     = 0
    max_depth = 0
    seg_start = 0.0
    components = IntervalType[]

    for (coord, kind) in events
        if kind == 0        # enter
            depth += 1
            if depth > max_depth
                max_depth = depth
                empty!(components)
                seg_start = coord
            elseif depth == max_depth
                seg_start = coord
            end
        else                # leave
            if depth == max_depth
                push!(components, intervaltype(seg_start, coord))
            end
            depth -= 1
        end
    end

    return (; region=components, depth=max_depth)
end

"""
    select_constant(region, rng=Random.GLOBAL_RNG) -> c

Pick the concrete additive constant from a max-overlap region, in priority order:
1. Returns `0` if zero lies in any component.
2. Samples uniformly over the union of components (interval chosen proportional to
   width, then a uniform point within it).
3. Returns `0` if the region is empty.

TODO: handle half-infinite bounds better, we currently use just the finite bound

"""
function select_constantorg(region::Vector{IntervalType}, rng::AbstractRNG=Random.GLOBAL_RNG)
    isempty(region) && return 0.0
    any(in_interval.(Ref(0.0), region)) && return 0.0
    finite_bounds = Float64[]    # finite endpoints of half-infinite intervals
    for iv in region
        lo, hi = iv.lo, iv.hi
        if isfinite(lo) && isfinite(hi)
            if abs(lo) >= 9.0e18 || abs(hi) >= 9.0e18
                # Very large finite bounds: hi-lo may overflow to Inf.
                # Record the bound closest to 0 so the !isfinite(total) fallback works.
                push!(finite_bounds, abs(lo) <= abs(hi) ? lo : hi)
            end
        elseif isfinite(hi)      # CI(-Inf, hi): 0 not in region so hi < 0
            push!(finite_bounds, hi)
        elseif isfinite(lo)      # CI(lo, Inf): 0 not in region so lo > 0
            push!(finite_bounds, lo)
        end
        # CI(-Inf, Inf): unreachable here — 0 ∈ CI(-Inf,Inf) would have returned above
    end
    widths = [iv.hi - iv.lo for iv in region]
    total  = sum(widths)
    if !isfinite(total)
        # Half-infinite intervals present and no integers found: return the finite bound
        return rand(rng, finite_bounds)
    end
    r      = rand(rng) * total
    cumw   = 0.0
    for (iv, w) in zip(region, widths)
        cumw += w
        if r <= cumw
            return iv.lo + rand(rng) * w
        end
    end
    iv = last(region.items)
    return iv.lo + rand(rng) * (iv.hi - iv.lo)
end

include("constant_sampler.jl")
function select_constant(region::Vector{IntervalType}, rng::AbstractRNG=Random.GLOBAL_RNG)
    isempty(region) && return 0.0
    any(in_interval.(Ref(0.0), region)) && return 0.0
    s = ConstantSampler()
    return draw_constant(rng, s, region)
end
 
const TWO_PI = 2π

"""
    fold_stab(carcs::IntervalVector, C=2π) -> (; region, depth)
    fold_stab(arcs::Vector{IntervalType}, C=2π)          # each arc its own case

Circular stab for the additive constant directly below `sin` (Situation 1 of
`specs/sine-inversion.md`). The constant is only defined mod `C`, so every arc is
reduced mod `C` and emitted **twice** — once where it lands, once translated by `-C` —
so a region wrapping past 0 is found whole (§3.3). Full arcs (`w >= C`, i.e. the case
is hit for every constant) are counted in `nfull` rather than materialised.

Depth counts **distinct cases**, not arcs: a sweep bumps depth only when a case's arc
count goes 0→1 and drops it 1→0. A case can legitimately contribute several arcs (rising
+ falling, or several arcs after a composed inversion) that may even *touch* mod `C`
after inner-rounding; counting cases rather than arcs keeps `depth == hits` without
relying on the per-case disjointness invariant (§1 / open-question #4) holding exactly.

Same `(; region, depth)` contract as `max_overlap_region`, with `nfull` folded in.
"""
fold_stab(arcs::Vector{IntervalType}, C::Float64 = TWO_PI) = fold_stab(IntervalVector(arcs), C)

function fold_stab(carcs::IntervalVector, C::Float64 = TWO_PI)
    # (coord, kind, case): kind 0 = enter, 1 = leave. Each non-full arc is emitted at
    # s and s-C, both tagged with the owning case so same-case overlap counts once.
    events = Tuple{Float64,Int,Int}[]
    nfull  = 0
    for i in eachindex(carcs)
        arcs = carcs[i]
        # A full arc (w >= C) means the case is hit for every constant. It then counts
        # once via nfull and contributes no sweep events — including for any *other*
        # (ordinary) arcs of the same case, which would otherwise double-count it.
        if any(arc -> sup(arc) - inf(arc) >= C, arcs)
            nfull += 1
            continue
        end
        for arc in arcs
            w = sup(arc) - inf(arc)
            s = mod(inf(arc), C)
            push!(events, (s, 0, i));     push!(events, (s + w, 1, i))
            push!(events, (s - C, 0, i)); push!(events, (s - C + w, 1, i))
        end
    end
    isempty(events) && return (; region = IntervalType[], depth = nfull)

    sort!(events)                          # coord asc; enter(0) before leave(1) at ties

    active     = zeros(Int, length(carcs)) # open-arc count per case
    depth      = 0                         # number of distinct active cases
    max_depth  = 0
    seg_start  = 0.0
    components = IntervalType[]
    for (coord, kind, case) in events
        if kind == 0                       # enter
            if active[case] == 0           # case becomes active
                depth += 1
                if depth > max_depth
                    max_depth = depth
                    empty!(components); seg_start = coord
                elseif depth == max_depth
                    seg_start = coord
                end
            end
            active[case] += 1
        else                               # leave
            if active[case] == 1           # case becomes inactive
                depth == max_depth && push!(components, intervaltype(seg_start, coord))
                depth -= 1
            end
            active[case] -= 1
        end
    end
    return (; region = components, depth = max_depth + nfull)
end

"""
    circular_hits(carcs, value, C=2π) -> Int

Count how many cases the constant `value` hits, given the **c-space arcs** `carcs`
(the same `(targets - evals)` arcs `fold_stab` folds). Case `i` hits when one of its
arcs — reduced mod `C` and emitted at `s` and `s-C` exactly as in `fold_stab` — covers
`value`, with **at most one hit per case** (the disjointness invariant of §1).

Using `fold_stab`'s own arcs and endpoint arithmetic (rather than recomputing from the
u-space `targets` and subtracting `evals`, which rounds differently) keeps this an
independent recount of *disjointness* without drifting from the fold's depth by a ULP at
a boundary. A hit against any `+kC` translate of the canonical arc is a genuine forward
hit because `sin` is `C`-periodic, so this is the truthful hit count for the constant
directly below `sin`; the linear `compute_hits` would miss translated hits.
"""
function circular_hits(carcs::IntervalVector, value::Real, C::Float64 = TWO_PI)
    hits = 0
    for i in eachindex(carcs)
        for arc in carcs[i]
            lo, hi = inf(arc), sup(arc)
            w = hi - lo
            if w >= C
                hits += 1
                break
            end
            s = mod(lo, C)
            if (s <= value <= s + w) || (s - C <= value <= s - C + w)
                hits += 1
                break
            end
        end
    end
    return hits
end

