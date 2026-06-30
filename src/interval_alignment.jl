using Intervals

_invalid_sentinel(::Type{T}) where {T<:Real} = Interval(typemax(T), typemax(T))
_is_invalid(t::Interval) = isinf(first(t)) && first(t) == last(t)

"""
    max_overlap_region(intervals) -> (; region, depth)

Sweep `intervals` to find the maximum depth `k` and return every closed range of
values contained in exactly `k` of the input intervals, together with `k`.

Intervals that are invalid sentinels (degenerate points at ±Inf, produced when a
target is unreachable) are silently filtered before the sweep. If all intervals are
invalid, returns `(; region=IntervalSet([Interval(Inf,Inf)]), depth=0)`.

Empty input returns `(; region=IntervalSet([]), depth=0)`.
"""
function max_overlap_region(intervals::Vector{Interval{T,Closed,Closed}}) where {T<:Real}
    isempty(intervals) && return (; region=IntervalSet(Interval{T,Closed,Closed}[]), depth=0)

    valid = filter(!_is_invalid, intervals)
    isempty(valid) && return (; region=IntervalSet([_invalid_sentinel(T)]), depth=0)

    intervals = valid

    # Encode: (coordinate, 0) = enter, (coordinate, 1) = leave
    events = Vector{Tuple{T,Int}}(undef, 2 * length(intervals))
    for (i, iv) in enumerate(intervals)
        events[2i-1] = (first(iv), 0)
        events[2i]   = (last(iv),  1)
    end
    sort!(events)

    depth     = 0
    max_depth = 0
    seg_start = zero(T)
    components = Interval{T,Closed,Closed}[]

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
                push!(components, Interval(seg_start, coord))
            end
            depth -= 1
        end
    end

    region = IntervalSet(components)
    return (; region, depth=max_depth)
end

"""
    select_constant(region) -> c

Pick the concrete additive constant from a max-overlap region:
- Returns `0` if zero lies in any component.
- Otherwise returns the midpoint of the widest component.
- Returns `0` if the region is empty.
"""
function select_constant(region::IntervalSet{Interval{T,Closed,Closed}}) where {T<:Real}
    isempty(region.items) && return zero(T)
    zero(T) in region && return zero(T)
    widest = argmax(iv -> last(iv) - first(iv), region.items)
    return (first(widest) + last(widest)) / 2
end
