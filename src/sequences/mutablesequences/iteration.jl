# Iteration
# =========
#
# Types and methods for iterating over biological sequences.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

# Ambiguous nucleotides iterator
# ------------------------------

struct AmbiguousNucleicAcidIterator{A<:NucleicAcidAlphabet}
    seq::MutableBioSequence{A}
end

ambiguous_positions(seq::MutableBioSequence) = AmbiguousNucleicAcidIterator(seq)

Base.start(it::AmbiguousNucleicAcidIterator) = find_next_ambiguous(it.seq, 1)
Base.done(it::AmbiguousNucleicAcidIterator, nextpos) = nextpos == 0
function Base.next(it::AmbiguousNucleicAcidIterator, nextpos)
    return nextpos, find_next_ambiguous(it.seq, nextpos + 1)
end

Base.iteratorsize(::AmbiguousNucleicAcidIterator) = Base.SizeUnknown()

function find_next_ambiguous(seq::MutableBioSequence{A}, i::Integer) where {A<:NucleicAcidAlphabet{2}}
    # no ambiguity
    return 0
end

function find_next_ambiguous(seq::MutableBioSequence{A}, from::Integer) where {A<:NucleicAcidAlphabet{4}}
    for i in max(from, 1):endof(seq)
        nt = inbounds_getindex(seq, i)
        if isambiguous(nt)
            return i
        end
    end
    # no ambiguity
    return 0
end
