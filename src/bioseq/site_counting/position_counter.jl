
struct PositionCounter{P<:Position} <: BitAlignedDoState
    count::Int
end

function head_update_state(state::PositionCounter{P}, ::Type{A},
                                x::UInt64, y::UInt64, k::Integer,
                                m::UInt64) where {P <: Position, A <: Alphabet}

    c = bp_chunk_count(P, A, x & m, y & m)
    c += state.count

    if bp_correct_emptyspace(P, A)
        #println("Correcting for emptyspace...")
        nempty = elems_per_chunk(A) - elems_per_x(k, A)
        #println("nempty: ", nempty)
        c -= nempty
        #println("counts: ", counts)
    end

    return PositionCounter{P}(c)
end

function update_state(state::PositionCounter{P}, ::Type{A}, x::UInt64,
                      y::UInt64) where {P <: Position, A <: Alphabet}

    c = bp_chunk_count(P, A, x, y)
    c += state.count
    return PositionCounter{P}(c)
end

function tail_update_state(state::PositionCounter{P}, ::Type{A},
                           x::UInt64, y::UInt64, k::Integer,
                           m::UInt64) where {P <: Position, A <: Alphabet}

    return head_update_state(state, A, x, y, k, m)
end
