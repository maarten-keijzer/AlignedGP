"""
    max_overlap_region(intervals) -> (; region, depth)

Sweep `intervals` (a vector of `CIntervals`) to find the maximum depth `k` and
return every closed range of values contained in exactly `k` of the sub-intervals,
together with `k`. Each `CIntervals` element may carry multiple disjoint ranges;
they are all flattened before the sweep.

Empty `CIntervals` elements (the error/invalid state) are silently skipped. If
all elements are empty or contain only invalid sentinels, returns
`(; region=CIntervals(), depth=0)`.
"""
function max_overlap_region(intervals::Vector{CIntervals})
    flat = flatten(intervals)
    isempty(flat) && return (; region=CIntervals(), depth=0)

    _is_useless(ci) = !isfinite(ci.lo) && ci.lo == ci.hi
    valid = filter(ci -> !_is_invalid(ci) && !_is_useless(ci), flat)
    isempty(valid) && return (; region=CIntervals(), depth=0)

    intervals = valid

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
    components = CInterval[]

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
                push!(components, CInterval(seg_start, coord))
            end
            depth -= 1
        end
    end

    region = CIntervals(components)
    return (; region, depth=max_depth)
end

"""
    select_constant(region, rng=Random.GLOBAL_RNG) -> c

Pick the concrete additive constant from a max-overlap region, in priority order:
1. Returns `0` if zero lies in any component.
2. Samples uniformly over the union of components (interval chosen proportional to
   width, then a uniform point within it).
3. Returns `0` if the region is empty.
"""
function select_constant(region::CIntervals, rng::AbstractRNG=Random.GLOBAL_RNG)
    isempty(region.items) && return 0.0
    0.0 in region && return 0.0
    finite_bounds = Float64[]    # finite endpoints of half-infinite intervals
    for iv in region.items
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
    widths = [iv.hi - iv.lo for iv in region.items]
    total  = sum(widths)
    if !isfinite(total)
        # Half-infinite intervals present and no integers found: return the finite bound
        return rand(rng, finite_bounds)
    end
    r      = rand(rng) * total
    cumw   = 0.0
    for (iv, w) in zip(region.items, widths)
        cumw += w
        if r <= cumw
            return iv.lo + rand(rng) * w
        end
    end
    iv = last(region.items)
    return iv.lo + rand(rng) * (iv.hi - iv.lo)
end
