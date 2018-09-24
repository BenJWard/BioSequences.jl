
abstract type ShortSequence{N} <: BioSequence end

encoded_data(x::ShortSequence{64}) = reinterpret(UInt64, x)
encoded_data(x::ShortSequence{32}) = reinterpret(UInt32, x)
encoded_data(x::ShortSequence{16}) = reinterpret(UInt16, x)
encoded_data(x::ShortSequence{8})  = reinterpret(UInt8, x)

encoded_data_eltype(x::ShortSequence) = eltype(encoded_data(x))
