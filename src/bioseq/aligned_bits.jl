
struct BitsChunk
    abits::UInt64
    bbits::UInt64
    ishead::Bool
    istail::Bool
    remaining::Int
end

@inline function remaining(bc::BitsChunk)
    return bc.remaining
end

@inline function abits(bc::BitsChunk)
    return bc.abits
end

@inline function bbits(bc::BitsChunk)
    return bc.bbits
end

@inline function bit_chunks(bc::BitsChunk)
    return abits(bc), bbits(bc)
end

@inline function ishead(bc::BitsChunk)
    return bc.ishead
end

@inline function istail(bc::BitsChunk)
    return bc.istail
end

function aligned_bits(a::BioSequence{A}, b::BioSequence{A})::Channel{BitsChunk} where {A}
    if length(a) > length(b)
        return Channel((c::Channel{BitsChunk}) -> _aligned_bits(c, b, a), ctype=BitsChunk)
    else
        return Channel((c::Channel{BitsChunk}) -> _aligned_bits(c, a, b), ctype=BitsChunk)
    end
end

function _aligned_bits(c::Channel{BitsChunk},
                      a::BioSequence{A},
                      b::BioSequence{A}) where {A}

    @assert length(a) ≤ length(b)
    nexta = bitindex(a, 1)
    stopa = bitindex(a, endof(a) + 1)
    nextb = bitindex(b, 1)
    stopb = bitindex(b, endof(b) + 1)
    #=
    Note that updating `nextb` or `nexta` by 64, increases the chunk
    index, but the `offset(nextb)` will remain the same.

    The first thing we need to sort out is to correctly align the head of
    sequence / subsequence `a`s data such that the offset of `nexta` is
    essentially reduced to 0.
    With sequence / subsequence `a` aligned, from there, we only need to
    worry about the alignment of sequence / subsequence `b` with respect to `a`.
    =#
    if nexta < stopa && offset(nexta) != 0
        # Here we shift the first data chunks to the right so as the first
        # nucleotide of the seq/subseq is the first nibble / pair of bits.
        x = a.data[index(nexta)] >> offset(nexta)
        y = b.data[index(nextb)] >> offset(nextb)
        # Check there is something to go and get from the next chunk of `b`.
        # Check like so: `64 - offset(nextb) < stopb - nextb`.
        # If this is not true of `b`, then it is certainly not true of `a`.
        if offset(nextb) > offset(nexta) && 64 - offset(nextb) < stopb - nextb
            y |= b.data[index(nextb) + 1] << (64 - offset(nextb))
        end
        #=
        Check if the chunk of `a` we are currently aligning contains the end
        of seq/subseq `a`.
        Check by seeing if the distance to the end of the sequence is smaller
        than the distance to the next word.

        If so, the mask used needs to be defined to account for this:
        `mask(stopa - nexta)`, otherwise the mask simply needs to be
        `mask(64 - offset(nexta))`.

        Finally, move position markers by k, meaning they move to either the
        next chunk, or the end of the sequence if it is in the current integer.
        =#
        rem_word = 64 - offset(nexta)
        rem_stop = stopa - nexta
        rem = ifelse(rem_stop < rem_word, rem_stop, rem_word)
        m = mask(rem)
        put!(c, BitsChunk(x & m, y & m, true, false, rem))
        nexta += rem
        nextb += rem
    end
    if offset(nextb) == 0 # data are aligned with each other.
        while stopa - nexta ≥ 64
            x = a.data[index(nexta)]
            y = b.data[index(nextb)]
            put!(c, BitsChunk(x, y, false, false, stopa - nexta))
            nexta += 64
            nextb += 64
        end
        if nexta < stopa # process the remaining tail.
            x = a.data[index(nexta)]
            y = b.data[index(nextb)]
            rem = stopa - nexta
            m = mask(rem)
            put!(c, BitsChunk(x & m, y & m, false, true, rem))
        end
    elseif nexta < stopa # b is unaligned with a.
        y = b.data[index(nextb)]
        nextb += 64
        while stopa - nexta ≥ 64
            x = a.data[index(nexta)]
            z = b.data[index(nextb)]
            y = y >> offset(nextb) | z << (64 - offset(nextb))
            put!(c, BitsChunk(x, y, false, false, stopa - nexta))
            y = z
            nexta += 64
            nextb += 64
        end
        if nexta < stopa # process the remaining tail.
            x = a.data[index(nexta)]
            y = y >> offset(nextb)
            if 64 - offset(nextb) < stopa - nexta
                y |= b.data[index(nextb)] << (64 - offset(nextb))
            end
            rem = stopa - nexta
            m = mask(rem)
            put!(c, BitsChunk(x & m, y & m, false, true, rem))
        end
    end
end
