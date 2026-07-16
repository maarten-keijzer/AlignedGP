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
 
const SEAM_TOL = 1e-12

function max_overlap_region_circular(arcs::Vector{IntervalType}, C::Float64 = 2π)

    isempty(arcs) && return (; region=IntervalType[], depth=0)

    # Encode: (coordinate, 0) = enter, (coordinate, 1) = leave.
    # A wrapping arc splits at the seam: tail [lo, C] plus head [0, hi-C].
    events = Vector{Tuple{Float64,Int}}()
    sizehint!(events, 2 * length(arcs) + 2)
    for iv in arcs
        lo, hi = bounds(iv)
        if hi <= C
            push!(events, (lo, 0));  push!(events, (hi, 1))
        else
            push!(events, (lo, 0));  push!(events, (C, 1))         # tail
            push!(events, (0.0, 0)); push!(events, (hi - C, 1))    # head
        end
    end
    sort!(events)

    depth     = 0
    max_depth = 0
    seg_start = 0.0
    components = IntervalType[]

    for (coord, kind) in events        # ← identical to the linear version
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

    return (; region=merge_seam(components, C), depth=max_depth)
end

# A component straddling the seam leaves the sweep as two pieces —
# one starting at 0, one ending at C. Physically they are one arc.
function merge_seam(components::Vector{IntervalType}, C::Float64)
    length(components) < 2 && return components
    i0 = findfirst(c -> inf(c) <= SEAM_TOL,     components)
    iC = findfirst(c -> sup(c) >= C - SEAM_TOL, components)
    (i0 === nothing || iC === nothing || i0 == iC) && return components

    merged = IntervalType[intervaltype(inf(components[iC]), C + sup(components[i0]))]
    for (k, c) in enumerate(components)
        (k == i0 || k == iC) && continue
        push!(merged, c)
    end
    merged
end

