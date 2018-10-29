
struct EachSkipmerIterator{SK <: Skipmer, UT <: Unsigned, SQ <: BioSequence}
    seq::SQ
    cycle_pos::Vector{UInt8}
    last_unknown::Vector{Int64}
    fkmer::Vector{UT}
    rkmer::Vector{UT}
end

skipmer_type(::Type{EachSkipmerIterator{SK}}) where SK <: Skipmer = SK
skipmer_type(x::EachSkipmerIterator) = skipmer_type(typeof(x))

kmersize(x::EachSkipmerIterator) = kmersize(typeof(x))

firstoffset(x::EachSkipmerIterator) = (kmersize(x) - 1) * 2

function kmer_mask(x::EachSkipmerIterator{SK,UT}) where {SK <: Skipmer, UT <: Unsigned}
    return (UT(1) << (kmersize(SK) * 2)) - 1
end

function firstoffset(x::EachSkipmerIterator)
    return (kmersize(SK) - 1) * 2
end

function EachSkipmerIterator(::Type{SK}, seq::SQ) where {SK <: Skipmer, SQ <: BioSequence}
    last_unknown = Vector{Int64}(undef, cycle_len(SK))
    fkmer = Vector{encoded_data_eltype(SK)}(undef, cycle_len(SK))
    rkmer = Vector{encoded_data_eltype(SK)}(undef, cycle_len(SK))
    cycle_pos = Vector{UInt8}(undef, cycle_len(SK))
    return EachSkipmerIterator{SK, encoded_data_eltype(SK), SQ}(seq, cycle_pos, last_unknown, fkmer, rkmer)
end
