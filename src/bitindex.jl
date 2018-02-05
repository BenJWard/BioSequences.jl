# BitIndex
# --------
#
# Utils for indexing bits in a vector of 64-bit integers (internal use only).
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioSequences.jl/blob/master/LICENSE.md

struct BitIndex{N, W}
    val::Int64
end

_index_shift(i::BitIndex{N, UInt64}) where N = 6
_index_shift(i::BitIndex{N, UInt32}) where N = 5
_index_shift(i::BitIndex{N, UInt16}) where N = 4
_index_shift(i::BitIndex{N, UInt8}) where N = 3
_offset_mask(i::BitIndex{N, W}) where {N, W} = UInt8(8 * sizeof(W)) - 0x01

#         index(i)-1        index(i)        index(i)+1
# ....|................|..X.............|................|....
#                          |<-offset(i)-|
#                      |<--- 64 bits -->|

@inline function bitindex(::Val{N}, ::Type{W}, index) where {N, W}
    return BitIndex{N, W}((index - 1) << trailing_zeros(N))
end

index(i::BitIndex) = (i.val >> _index_shift(i)) + 1
offset(i::BitIndex) = i.val & _offset_mask(i)

Base.:+(i::BitIndex{N, W}, n::Int) where {N, W} = BitIndex{N, W}(i.val + n)
Base.:-(i::BitIndex{N, W}, n::Int) where {N, W} = BitIndex{N, W}(i.val - n)
Base.:-(i1::BitIndex, i2::BitIndex) = i1.val - i2.val
Base.:(==)(i1::BitIndex, i2::BitIndex) = i1.val == i2.val
Base.isless(i1::BitIndex, i2::BitIndex) = isless(i1.val, i2.val)
Base.cmp(i1::BitIndex, i2::BitIndex) = cmp(i1.val, i2.val)
Base.start(i::BitIndex) = 1
Base.done(i::BitIndex, s) = s > 2
Base.next(i::BitIndex, s) = ifelse(s == 1, (index(i), 2), (offset(i), 3))
Base.show(io::IO, i::BitIndex) = print(io, '(', index(i), ", ", offset(i), ')')

# Create a bit mask that fills least significant `n` bits (`n` must be a
# non-negative integer).
bitmask(::Type{T}, n::Integer) where {T} = (one(T) << n) - one(T)
bitmask(::Type{T}, ::Val{N}) where {T, N} = (one(T) << N) - one(T)
bitmask(bidx::BitIndex{N, W}) where {N, W} = bitmask(W, N)

#TODO: Possibly delete later, redunandt with BitIndex redesign.
@inline function bitmask(::Type{A}, ::Type{U}) where {A <: Alphabet, U <: Unsigned}
    ba = bits_per_symbol(A)
    return bitmask(U, ba)
end
bitmask(::Type{A}) where {A <: Alphabet} = bitmask(bits_per_symbol(A))
bitmask(n::Integer) = bitmask(UInt64, n)

@inline function extract_encoded_symbol(bidx::BitIndex, data)
    @inbounds chunk = data[index(bidx)]
    offchunk = chunk >> offset(bidx)
    return offchunk & bitmask(bidx)
end
