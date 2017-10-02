


struct PositionCounter{P<:Position} <: BitAlignedDoAlgorithm end

@inline bp_counter_type(::Type{<:Site}, ::Type{<:Alphabet}) = Int

@inline bp_start_counter(::Type{S}, ::Type{A}) where {S,A<:Alphabet} = zero(bp_counter_type(S, A))

@inline bp_update_counter(acc::Int, up::Int) = acc + up

@inline bp_emptyspace_correction(nempty::Int, count::Int) = count - nempty

struct CorrectEmptyspace{v} end

correct_emptyspace(::Type{<:Site}, ::Type{<:Alphabet}) = CorrectEmptyspace{false}

for s in (Match, Gap)
    @eval correct_emptyspace(::Type{$s}, ::Type{A}) where {A<:NucAlphs} = CorrectEmptyspace{true}
end
