# Bitparallel counting
# --------------------

# Now we overload some internal methods for each S <: Site.
# So as different types of site can be counted in bit-parallel
# manner.

for A in (DNAAlphabet, RNAAlphabet)
    @eval begin

        # Gaps
        @inline function bp_chunk_count(::Type{Gap}, ::Type{$A{4}}, x::UInt64)
            return count_zero_nibbles(x)
        end

        @inline function bp_chunk_count(::Type{Gap}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            # Count the gaps in a, count the gaps in b, subtract the number of shared gaps.
            return count_zero_nibbles(a) + count_zero_nibbles(b) - count_zero_nibbles(a | b)
        end

        # Certain
        @inline function bp_chunk_count(::Type{Certain}, ::Type{$A{4}}, x::UInt64)
            x = enumerate_nibbles(x)
            x = x ⊻ 0x1111111111111111
            return count_zero_nibbles(x)
        end

        @inline function bp_chunk_count(::Type{Certain}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            x = enumerate_nibbles(a) ⊻ 0x1111111111111111
            y = enumerate_nibbles(b) ⊻ 0x1111111111111111
            return count_zero_nibbles(x | y)
        end

        # Ambiguous
        @inline function bp_chunk_count(::Type{Ambiguous}, ::Type{$A{4}}, x::UInt64)
            return count_nonzero_nibbles(enumerate_nibbles(x) & 0xEEEEEEEEEEEEEEEE)
        end

        @inline function bp_chunk_count(::Type{Ambiguous}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            return count_nonzero_nibbles((enumerate_nibbles(a) | enumerate_nibbles(b)) & 0xEEEEEEEEEEEEEEEE)
        end

        # Match
        @inline function bp_chunk_count(::Type{Match}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            return count_zero_nibbles(a ⊻ b)
        end

        @inline function bp_chunk_count(::Type{Match}, ::Type{$A{2}}, a::UInt64, b::UInt64)
            return count_zero_bitpairs(a ⊻ b)
        end

        # Mismatch
        @inline function bp_chunk_count(::Type{Mismatch}, ::Type{$A{4}}, a::UInt64, b::UInt64)
            return count_nonzero_nibbles(a ⊻ b)
        end

        @inline function bp_chunk_count(::Type{Mismatch}, ::Type{$A{2}}, a::UInt64, b::UInt64)
            return count_nonzero_bitpairs(a ⊻ b)
        end
    end
end
