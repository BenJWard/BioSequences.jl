# Site types
# ----------

abstract type Site end
abstract type Position <: Site end

"""
A `Certain` site describes a site where both of two aligned sites are not an
ambiguity symbol or a gap.
"""
struct Certain <: Position end

"""
An `Gap` site describes a site where either of two aligned sites are a
gap symbol '-'.
"""
struct Gap <: Position end

"""
An `Ambiguous` site describes a site where either of two aligned sites are an
ambiguity symbol.
"""
struct Ambiguous <: Position end

"""
A `Match` site describes a site where two aligned nucleotides are the
same biological symbol.
"""
struct Match <: Position end

"""
A `Mismatch` site describes a site where two aligned nucleotides are not the
same biological symbol.
"""
struct Mismatch <: Position end
