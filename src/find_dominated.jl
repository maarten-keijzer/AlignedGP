
# Find ONE BitVector that is dominated-or-equal by another in the set.
#
# "vecs[i] dominated-or-equal by vecs[j]"  <=>  vecs[i] ⊆ vecs[j]:
# every bit set in i is also set in j.  Chunk test: a & ~b == 0 per UInt64.
#
# Relies on two BitVector guarantees:
#   * `.chunks::Vector{UInt64}` packs the bits, 64 per word;
#   * unused high bits of the final chunk are kept zero,
# so the chunk-wise test is correct as long as all vectors share one length.

@inline function is_subset(a::BitVector, b::BitVector)
    ca, cb = a.chunks, b.chunks
    @inbounds for k in eachindex(ca)
        (ca[k] & ~cb[k]) != 0 && return false   # a has a bit b lacks
    end
    return true
end

popcnt(b::BitVector) = @inbounds mapreduce(count_ones, +, b.chunks; init = 0)

# One-word summary: bit (k-1)&63 set iff chunk k is nonzero (OR-bucketed >64 words).
# Necessary condition for subset: sig(a) & ~sig(b) == 0.  Rejects most misses cheaply.
function signature(b::BitVector)
    s = zero(UInt64)
    @inbounds for (k, w) in enumerate(b.chunks)
        w != 0 && (s |= one(UInt64) << ((k - 1) & 63))
    end
    return s
end


"""
    find_dominated(vecs) -> Union{Tuple{Int,Int}, Nothing}

Return `(i, j)` for the first `vecs[i]` dominated-or-equal by `vecs[j]`
(`vecs[i] ⊆ vecs[j]`), or `nothing` if the set is an antichain.
All bitvectors are assumed to share one length.
"""
function find_dominated(vecs::AbstractVector{BitVector})
    n = length(vecs)
    n < 2 && return nothing

    # Fast path: exact duplicates dominate each other. O(n·w).
    seen = Dict{Vector{UInt64}, Int}()
    @inbounds for i in 1:n
        key = vecs[i].chunks
        haskey(seen, key) && return (i, seen[key])
        seen[key] = i
    end

    pc    = [popcnt(v)    for v in vecs]
    sig   = [signature(v) for v in vecs]
    order = sortperm(pc)                        # ascending popcount

    @inbounds for ai in 1:n
        a    = order[ai]                        # candidate dominated: fewest bits first
        A    = vecs[a]; pcA = pc[a]; sigA = sig[a]
        for bi in n:-1:1                        # candidate dominator: most bits first
            b = order[bi]
            pc[b] < pcA && break               # nothing left can be a superset
            b == a && continue
            (sigA & ~sig[b]) != 0 && continue  # signature prefilter
            is_subset(A, vecs[b]) && return (a, b)
        end
    end
    return nothing
end

# Your layout: Vector{Tree}, each Tree with a `hits::BitVector` field.
# Returned indices index straight into `trees` (order is preserved).
find_dominated_trees(trees::AbstractVector) =
    find_dominated([t.hits for t in trees])


