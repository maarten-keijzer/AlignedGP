

    emit(start, end):            # arc given as a real interval, length ≤ 2π
        s = start mod 2π
        e = s + (end - start)
        if e <= 2π:  return [(s, e)]
        else:        return [(s, 2π), (0, e - 2π)]      # seam split

    invert_sin_interval(ℓ, h):   # -> ordered disjoint intervals ⊆ [0, 2π)
    
        ℓ' = max(ℓ, -1);  h' = min(h, 1)
        
        if ℓ' > h':            return []                 # out of range: unwinnable
        
        if ℓ' == -1 and h' == 1: return [(0, 2π)]        # whole circle

        aL = asin(ℓ');  aH = asin(h')                    # both in [-π/2, π/2], aL ≤ aH

        pieces  = emit(aL,     aH)                        # rising branch
        pieces += emit(π - aH, π - aL)                    # falling branch
        return coalesce_on_line(sort(pieces))

