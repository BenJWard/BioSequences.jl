# Conversion & Promotion
# ======================
#
# Conversion methods for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Promotion
# ---------

for alph in (DNAAlphabet, RNAAlphabet)
    @eval function Base.promote_rule(::Type{MutableBioSequence{A}}, ::Type{MutableBioSequence{B}}) where {A<:$alph,B<:$alph}
        return MutableBioSequence{promote_rule(A,B)}
    end
end

# Conversion
# ----------

# Conversion between sequences of different alphabet size.
for A in [DNAAlphabet, RNAAlphabet]

    # Convert from a 4 bit encoding to a 2 bit encoding.
    @eval function Base.convert(::Type{MutableBioSequence{$(A{2})}}, seq::MutableBioSequence{$(A{4})})
        # TODO: make it faster with bit-parallel algorithm
        newseq = MutableBioSequence{$(A{2})}(length(seq))
        for (i, x) in enumerate(seq)
            unsafe_setindex!(newseq, x, i)
        end
        return newseq
    end

    # Convert from a 2 bit encoding to a 4 bit encoding.
    @eval function Base.convert(::Type{MutableBioSequence{$(A{4})}}, seq::MutableBioSequence{$(A{2})})
        newseq = MutableBioSequence{$(A{4})}(length(seq))
        for (i, x) in enumerate(seq)
            unsafe_setindex!(newseq, x, i)
        end
        return newseq
    end
end

# Conversion between DNA and RNA sequences.
for (A1, A2) in [(DNAAlphabet, RNAAlphabet), (RNAAlphabet, DNAAlphabet)], n in (2, 4)
    # NOTE: assumes that binary representation is identical between DNA and RNA
    @eval function Base.convert(::Type{MutableBioSequence{$(A1{n})}},
                                seq::MutableBioSequence{$(A2{n})})
        newseq = MutableBioSequence{$(A1{n})}(seq.data, seq.part, true)
        seq.shared = true
        return newseq
    end
end

# Convert from a DNA or RNA vector to a BioSequence.
function Base.convert(::Type{MutableBioSequence{A}}, seq::AbstractVector{DNA}) where {A<:DNAAlphabet}
    return MutableBioSequence{A}(seq, 1, endof(seq))
end
function Base.convert(::Type{MutableBioSequence{A}}, seq::AbstractVector{RNA}) where {A<:RNAAlphabet}
    return MutableBioSequence{A}(seq, 1, endof(seq))
end
function Base.convert(::Type{AminoAcidSequence}, seq::AbstractVector{AminoAcid})
    return AminoAcidSequence(seq, 1, endof(seq))
end

# Convert from a BioSequence to to a DNA or RNA vector
Base.convert(::Type{Vector}, seq::MutableBioSequence) = collect(seq)
function Base.convert(::Type{Vector{DNA}}, seq::MutableBioSequence{<:DNAAlphabet})
    return collect(seq)
end
function Base.convert(::Type{Vector{RNA}}, seq::MutableBioSequence{<:RNAAlphabet})
    return collect(seq)
end
Base.convert(::Type{Vector{AminoAcid}}, seq::AminoAcidSequence) = collect(seq)

# Covert from a string to a BioSequence and _vice versa_.
function Base.convert(::Type{S}, seq::MutableBioSequence) where {S<:AbstractString}
    return convert(S, [Char(x) for x in seq])
end
Base.convert(::Type{MutableBioSequence{A}}, seq::S) where {S<:AbstractString,A} = MutableBioSequence{A}(seq)
