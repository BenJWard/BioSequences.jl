
include("bitindex.jl")

@inline function reversebits(x::UInt64, ::BitsPerSymbol{2})
     x = (x & 0x3333333333333333) <<  2 | (x & 0xCCCCCCCCCCCCCCCC) >>  2
     x = reversebits(x, BitsPerSymbol{4}())
     return x
end

@inline function reversebits(x::UInt64, ::BitsPerSymbol{4})
     x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>  4
     x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>  8
     x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >> 16
     x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >> 32
     return x
end

function count_a(x::UInt64)
    xinv = ~x
    return count_ones(((xinv >>> 1) & xinv) & Twiddle.repeatbyte(0x55))
end
count_c(x::UInt64) = count_ones((((~x) >>> 1)    &   x       ) & 0x5555555555555555)
count_g(x::UInt64) = count_ones(((  x  >>> 1)    & (~x      )) & 0x5555555555555555)
count_t(x::UInt64) = count_ones((   x            & ( x >>> 1)) & 0x5555555555555555)