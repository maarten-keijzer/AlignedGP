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
    fold_stab(arcs, C=2π) -> (; region, depth)

Circular stab for the additive constant directly below `sin` (Situation 1 of
`specs/sine-inversion.md`). The constant is only defined mod `C`, so every arc is
reduced mod `C` and emitted **twice** — once where it lands, once translated by `-C`.
The duplicate lets the plain linear `max_overlap_region` find a region that wraps past
0 as one unbroken interval (§3.3). Full arcs (`w >= C`, i.e. the case is hit for every
constant) are counted separately in `nfull` rather than materialised, since both copies
of a `≥C` arc would cover the same point.

Same `(; region, depth)` contract as `max_overlap_region`, with `nfull` folded into
`depth`. Widths are taken **before** reduction and carried; endpoints are never reduced
separately (that would split an arc across the seam).
"""
function fold_stab(arcs::Vector{IntervalType}, C::Float64 = TWO_PI)
    out   = IntervalType[]
    nfull = 0
    for iv in arcs
        lo, hi = inf(iv), sup(iv)
        w = hi - lo
        if w >= C
            nfull += 1                       # every constant works for this case
        else
            s = mod(lo, C)
            push!(out, intervaltype(s, s + w))
            push!(out, intervaltype(s - C, s - C + w))
        end
    end
    res = max_overlap_region(out)
    return (; region = res.region, depth = res.depth + nfull)
end

"""
    circular_hits(evals, targets, C=2π) -> Int

Count, for the chosen constant already folded into `evals`, how many cases are hit
**mod `C`**: case `i` hits when some arc of `targets[i]` covers `evals[i]` up to a
multiple of `C`. This is the truthful hit count for a constant directly below `sin`
(a hit against any `+kC` translate of the canonical arc is a genuine forward hit,
because `sin` is `C`-periodic); the linear `compute_hits` would miss translated hits.
"""
function circular_hits(evals::Vector{<:Real}, targets::IntervalVector, C::Float64 = TWO_PI)
    hits = 0
    for i in eachindex(targets)
        v = evals[i]
        for arc in targets[i]
            lo, hi = inf(arc), sup(arc)
            w = hi - lo
            if w >= C || mod(v - lo, C) <= w
                hits += 1
                break
            end
        end
    end
    return hits
end

