# Predicates & comparisons
# ------------------------

function Base.cmp(seq1::BioSequence, seq2::BioSequence)
    m = lastindex(seq1)
    n = lastindex(seq2)
    for i in 1:min(m, n)
        c = cmp(inbounds_getindex(seq1, i),
                inbounds_getindex(seq2, i))
        if c != 0
            return c
        end
    end
    return cmp(m, n)
end

function Base.:(==)(seq1::BioSequence, seq2::BioSequence)
    return eltype(seq1)    == eltype(seq2) &&
           length(seq1)    == length(seq2) &&
           cmp(seq1, seq2) == 0
end

Base.isless(seq1::BioSequence, seq2::BioSequence) = cmp(seq1, seq2) < 0

function Base.isempty(seq::BioSequence)
    return length(seq) == 0
end