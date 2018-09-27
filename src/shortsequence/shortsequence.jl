# ShortSequence
# =============
#
# A compact short sequence type.
#
# Sometimes applications do not use long sequences, but short ones.
# Such applications include tools that analyse and assemble kmers generated from
# reads, or tools that analyse codons: short 3bp long sequences that correspond to
# amino acids when a gene is translated into its protein product.
# It is conveinient to represent such sequences using single 'words'.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Representation
# --------------
#
# Four kinds of nucleotides are encoded as follows:
#
#   nucleotide | binary
#   ---------- | ------
#       A      |   00
#       C      |   01
#       G      |   10
#     T / U    |   11
#
# NucleicAcids are filled from MSBs to LSBs and right-aligned so that all k-mers
# are lexicographically ordered. For example, the memory layout of "TACG" is:
#   64-bit: 0b 00 00 … 00 11 00 01 10
#    4-mer:                T  A  C  G


abstract type ShortSequence{N} <: BioSequence end

encoded_data(x::ShortSequence{64}) = reinterpret(UInt64, x)
encoded_data(x::ShortSequence{32}) = reinterpret(UInt32, x)
encoded_data(x::ShortSequence{16}) = reinterpret(UInt16, x)
encoded_data(x::ShortSequence{8})  = reinterpret(UInt8, x)

encoded_data_eltype(x::ShortSequence) = eltype(encoded_data(x))


"""
    complement(x::T) where {T <: ShortSequence}

Return the complement of a short sequence type `x`.
"""
function BioSymbols.complement(x::T) where {T <: ShortSequence}
    return T(~encoded_data(x))
end


"""
    reverse(x::T) where {T <: ShortSequence}

Return the reverse of short sequence type variable `x`.
"""
function Base.reverse(x::T) where {T <: ShortSequence}
    bits = encoded_data(x)
    rbits = reversebits(bits, BitsPerSymbol{2}())
    return T(rbits >> (sizeof(bits) * 8 - 2 * length(x)))
end


"""
    reverse_complement(x::ShortSequence)

Return the reverse complement of `x`.
"""
reverse_complement(x::ShortSequence) = complement(reverse(x))


function swap(kmer::T, i, j) where {T <: ShortSequence}
    i = 2 * length(kmer) - 2i
    j = 2 * length(kmer) - 2j
    b = encoded_data(kmer)
    x = ((b >> i) ⊻ (b >> j)) & encoded_data_eltype(kmer)(0x03)
    return T(b ⊻ ((x << i) | (x << j)))
end



