# count_sites_bitpar.jl
# =====================
#
# Counting sites in a bitparallel fashion.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

@inline bp_counter_type(::Type{<:Site}, ::Type{<:Alphabet}) = Int
@inline bp_start_counter(::Type{S}, ::Type{A}) where {S,A<:Alphabet} = zero(bp_counter_type(S, A))
@inline bp_update_counter(acc::Int, up::Int) = acc + up
@inline bp_correct_emptyspace(::Type{<:Site}, ::Type{<:Alphabet}) = false
@inline bp_emptyspace_correction(nempty::Int, count::Int) = count - nempty

@inline function nempty(::Type{A}, remaining::Int)::Int where {A<:Alphabet}
    return div(64, bitsof(A)) - div(remaining, bitsof(A))
end

function bitpar_counter(::Type{S}, a::BioSequence{A}, b::BioSequence{A}) where {S<:Site,A<:NucAlphs}
    bits_channel = aligned_bits(a, b)
    counts = bp_start_counter(S, A)
    block = take!(bits_channel)
    x, y = bit_chunks(block)
    counts = bp_update_counter(counts, bp_chunk_count(S, A, x, y))
    if ishead(block) && bp_correct_emptyspace(S, A)
        n_empty = nempty(A, remaining(block))
        counts = bp_emptyspace_correction(n_empty, counts)
    end
    for block in bits_channel
        x, y = bit_chunks(block)
        counts = bp_update_counter(counts, bp_chunk_count(S, A, x, y))
    end
    if istail(block) && bp_correct_emptyspace(S, A)
        n_empty = nempty(A, remaining(block))
        counts = bp_emptyspace_correction(n_empty, counts)
    end
    return counts
end
