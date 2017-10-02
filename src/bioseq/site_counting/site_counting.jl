# site_counting.jl
# ================
#
# Site counting framework for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md


# Site types
# ----------

abstract type Site end
abstract type Position <: Site end

"""
A `Certain` site describes a site where both of two aligned sites are not an
ambiguity symbol or a gap.
"""
immutable Certain <: Position end

"""
An `Gap` site describes a site where either of two aligned sites are a
gap symbol '-'.
"""
immutable Gap <: Position end

"""
An `Ambiguous` site describes a site where either of two aligned sites are an
ambiguity symbol.
"""
immutable Ambiguous <: Position end

"""
A `Match` site describes a site where two aligned nucleotides are the
same biological symbol.
"""
immutable Match <: Position end

"""
A `Mismatch` site describes a site where two aligned nucleotides are not the
same biological symbol.
"""
immutable Mismatch <: Position end

include("bitops.jl")
include("position_counter.jl")

# Base.count methods
# ---------------------------

# The main bit-parallel counting method.
# Requires that the alphabet of two sequences is the same.
function Base.count(::Type{P}, a::BioSequence{A}, b::BioSequence{A}) where {P<:Position, A<:NucAlphs}
    state = PositionCounter{P}(0)
    state = bitaligned_do(state, a, b)
    return state.count
end

# Some specific edge cases...
for A in (DNAAlphabet, RNAAlphabet)

    # Specific count methods for some edge cases regarding counting sites in
    # 2 bit encoded DNA and RNA sequences.
    seqtype = BioSequence{A{2}}
    @eval begin
        function Base.count(::Type{Certain}, a::$seqtype, b::$seqtype)
            return min(length(a), length(b))
        end
        Base.count(::Type{Gap}, a::$seqtype, b::$seqtype) = 0
        Base.count(::Type{Ambiguous}, a::$seqtype, b::$seqtype) = 0
    end
end

# Specific Base.count sliding window methods
# ------------------------------------------

@inline bp_counter_type(::Type{<:Position}, ::Type{<:Alphabet}) = Int

function Base.count(::Type{P}, a::BioSequence{A}, b::BioSequence{A}, width::Int, step::Int) where {P<:Position,A<:NucAlphs}
    len = min(length(a), length(b))
    ritr = StepRange(width, step, len)
    width -= 1
    results = Vector{IntervalValue{Int,bp_counter_type(P, A)}}(length(ritr))
    r = 1
    @inbounds for i in ritr
        idx = (i - width):i
        results[r] = IntervalValue(first(idx), last(idx), count(P, a[idx], b[idx]))
        r += 1
    end
    return results
end

# Specific count_pairwise methods
# -------------------------------

@inline diag_val(::Type{Int}) = Int(0)

function count_pairwise(::Type{P}, seqs::Vararg{BioSequence{A},N}) where {P<:Position,A<:NucAlphs,N}
    @assert N >= 2 "At least two sequences are required."
    counts = Matrix{bp_counter_type(P, A)}(N, N)
    for i in 1:N
        counts[i,i] = diag_val(eltype(counts))
        for j in (i+1):N
            counts[i,j] = counts[j,i] = count(P, seqs[i], seqs[j])
        end
    end
    return counts
end


# General, non-specific Base.count methods
# ----------------------------------------

"""
    Base.count{S<:Site}(::Type{S}, a::BioSequence, b::BioSequence)

Count the number of sites of type `S`, between two sequences a and b.
"""
@inline function Base.count(::Type{S}, a::BioSequence, b::BioSequence) where {S<:Site}
    seqs = promote(a, b)
    return count(S, seqs...)
end

"""
    Base.count{S<:Site}(::Type{S}, a::BioSequence, b::BioSequence, width::Int, step::Int)

Count the number of sites of type `S`, between two sequences a and b, using a
sliding window of a given `width`, and `step`.
The `width` and `step` must be provided as number of base pairs.
"""
function Base.count(::Type{S}, a::BioSequence, b::BioSequence, width::Int, step::Int) where {S<:Site}
    seqs = promote(a, b)
    return count(S, seqs..., width, step)
end

"""
    count_pairwise{S<:Site,N}(::Type{S}, seqs::Vararg{BioSequence,N})

Count the number of sites of type `S`, between each possible pair of the sequences
provided (`seqs`).

This returns a symmetric matrix contining the output for each pairwise
operation.
"""
function count_pairwise(::Type{S}, seqs::Vararg{BioSequence,N}) where {S<:Site,N}
    seqs = promote(seqs...)
    return count_pairwise(S, seqs...)
end
