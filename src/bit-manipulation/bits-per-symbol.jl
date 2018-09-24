
"""
The number of bits required to represent a packed symbol in a vector of bits.
"""
struct BitsPerSymbol{N} end

BitsPerSymbol(::A) where A <: NucleicAcidAlphabet{2} = BitsPerSymbol{2}()
BitsPerSymbol(::A) where A <: NucleicAcidAlphabet{4} = BitsPerSymbol{4}()
BitsPerSymbol(::AminoAcidAlphabet) = BitsPerSymbol{8}()
BitsPerSymbol(::CharAlphabet) = BitsPerSymbol{32}()
BitsPerSymbol(::VoidAlphabet) = BitsPerSymbol{0}()

bits_per_symbol(::BitsPerSymbol{N}) where N = N
bits_per_symbol(::A) where A <: Alphabet = bits_per_symbol(BitsPerSymbol(A()))