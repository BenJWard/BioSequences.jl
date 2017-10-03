# site_counting.jl
# ================
#
# Site counting framework for biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

include("site_types.jl")
include("bitops.jl")
include("position_counter.jl")

function init_state(::Type{PositionCounter{P}},
    ::Type{A}) where {P <: Position, A <: Alphabet}
    return bp_start_counter(P, A)
end

function head_update_state(::Type{PositionCounter{P}}, state, ::Type{A},
    x::UInt64, y::UInt64, k::Integer) where {P <: Position, A <: Alphabet}

    return _head_update_state(P, state, A, x, y, k, correct_emptyspace(P, A))
end

function _head_update_state(::Type{PositionCounter{P}}, state, ::Type{A},
    x::UInt64, y::UInt64, k::Integer,
    ::Type{CorrectEmptyspace{false}}) where {P <: Position, A <: Alphabet}

    m = mask(k)
    c = bp_update_counter(state, bp_chunk_count(P, A, x & m, y & m))
    return c
end

function _head_update_state(::Type{PositionCounter{P}}, state, ::Type{A},
    x::UInt64, y::UInt64, k::Integer,
    ::Type{CorrectEmptyspace{true}}) where {P <: Position, A <: Alphabet}

    m = mask(k)
    c = bp_update_counter(state, bp_chunk_count(P, A, x & m, y & m))

    #println("Correcting for emptyspace...")
    nempty = elems_per_chunk(A) - elems_per_x(k, A)
    #println("nempty: ", nempty)
    c = bp_emptyspace_correction(nempty, c)
    #println("counts: ", counts)
    return c
end

function update_state(::Type{PositionCounter{P}}, state, ::Type{A}, x::UInt64,
    y::UInt64) where {P <: Position, A <: Alphabet}

    return bp_update_counter(state, bp_chunk_count(P, A, x, y))
end

function tail_update_state(::Type{PositionCounter{P}}, state, ::Type{A},
    x::UInt64, y::UInt64, k::Integer) where {P <: Position, A <: Alphabet}

    return head_update_state(PositionCounter{P}, state, A, x, y, k)
end

# Base.count methods
# ---------------------------

# The main bit-parallel counting method.
# Requires that the alphabet of two sequences is the same.
function Base.count(::Type{P}, a::BioSequence{A}, b::BioSequence{A}) where {P <: Position, A <: Alphabet}
    return bitaligned_do(PositionCounter{P}, a, b)
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
