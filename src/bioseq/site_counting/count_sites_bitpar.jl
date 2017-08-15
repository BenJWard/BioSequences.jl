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

function bitpar_counter(::Type{S}, a::BioSequence{A}, b::BioSequence{A}) where {S<:Site,A}
    #println("Open channel.")
    bits_channel = aligned_bits(a, b)
    counts = bp_start_counter(S, A)
    #println("Take first block.")
    block = take!(bits_channel)
    #println("Block:")
    #println(block)
    x, y = bit_chunks(block)
    counts = bp_update_counter(counts, bp_chunk_count(S, A, x, y))
    #println("First counts: ", counts)
    #println("Was it an aligned head block: ", ishead(block))
    #println("Correct for emptyspace: ", bp_correct_emptyspace(S, A))
    if ishead(block) && bp_correct_emptyspace(S, A)
        nempty = div(64, bitsof(A)) - div(remaining(block), bitsof(A))
        #println("Empty: ", nempty)
        counts = bp_emptyspace_correction(nempty, counts)
        #println("Corrected counts: ", counts)
    end
    while isopen(bits_channel)
        #println("In the body code block.")
        block = take!(bits_channel)
        #println("Block: ")
        #println(block)
        x, y = bit_chunks(block)
        counts = bp_update_counter(counts, bp_chunk_count(S, A, x, y))
        #println("Counts: ", counts)
    end
    if istail(block)
        #println("In the tail code block.")
        if bp_correct_emptyspace(S, A)
            nempty = div(64, bitsof(A)) - div(remaining(block), bitsof(A))
            counts = bp_emptyspace_correction(nempty, counts)
        end
    end
    return counts
end
