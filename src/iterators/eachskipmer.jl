
struct EachSkipmerIterator{SK <: Skipmer, UT <: Unsigned, SQ <: BioSequence}
    seq::SQ
    cycle_pos::Vector{UInt8}
    last_unknown::Vector{Int64}
    fkmer::Vector{UT}
    rkmer::Vector{UT}
end

@inline Base.eltype(::Type{EachSkipmerIterator{SK, UT, SQ}}) where {SK <: Skipmer, UT <: Unsigned, SQ <: BioSequence} = SK

@inline kmersize(x::EachSkipmerIterator) = kmersize(eltype(x))

@inline firstoffset(x::EachSkipmerIterator) = (kmersize(x) - 1) * 2

@inline function kmer_mask(x::EachSkipmerIterator{SK,UT}) where {SK <: Skipmer, UT <: Unsigned}
    return (UT(1) << (kmersize(SK) * 2)) - 1
end

function EachSkipmerIterator(::Type{SK}, seq::SQ) where {SK <: Skipmer, SQ <: BioSequence}
    @assert span(SK) <= length(seq)
    last_unknown = Vector{Int64}(undef, cycle_len(SK))
    fkmer = Vector{encoded_data_eltype(SK)}(undef, cycle_len(SK))
    rkmer = Vector{encoded_data_eltype(SK)}(undef, cycle_len(SK))
    cycle_pos = Vector{UInt8}(undef, cycle_len(SK))
    return EachSkipmerIterator{SK, encoded_data_eltype(SK), SQ}(seq, cycle_pos, last_unknown, fkmer, rkmer)
end

@inline function init_iterator!(it::EachSkipmerIterator)
    N = cycle_len(eltype(it))
    @inbounds for i in 1:N
        it.cycle_pos[i] = N - i
        it.last_unknown[i] = -1
        it.fkmer[i] = 0
        it.rkmer[i] = 0
    end
end

function Base.iterate(it::EachSkipmerIterator)
    N = cycle_len(eltype(it))
    M = bases_per_cycle(eltype(it))
    S = span(eltype(it))
    
    init_iterator!(it)
    
    for pos in 1:S
        for ni in 1:N
            it.cycle_pos[ni] += 1
            if it.cycle_pos[ni] == N
                it.cycle_pos[ni] = 0
            end
            
            if it.cycle_pos[ni] < M
                println("Sequence position: ", pos, ", Phase: ", ni)
                fbits = extract_encoded_symbol(bitindex(it.seq, pos), encoded_data(it.seq))
                rbits = ~fbits & 0x03
                it.fkmer[ni] = ((it.fkmer[ni] << 2) | fbits) & kmer_mask(it)
                it.rkmer[ni] = (it.rkmer[ni] >> 2) | (UInt64(rbits) << firstoffset(it))
            end
        end
    end
    
    fkmer = first(it.fkmer)
    rkmer = first(it.rkmer)
    
    outkmer = ifelse(fkmer < rkmer, fkmer, rkmer)
    
    return reinterpret(eltype(it), outkmer), (S + 1, UInt(1))
    
end

function Base.iterate(it::EachSkipmerIterator, state::Tuple{UInt, UInt})
    pos = state[1]
    fi  = state[2]
    N = cycle_len(eltype(it))
    M = bases_per_cycle(eltype(it))
    
    if pos > lastindex(it.seq)
        return nothing
    end
    
    # Copied from other method
    for ni in 1:N
        it.cycle_pos[ni] += 1
        if it.cycle_pos[ni] == N
            it.cycle_pos[ni] = 0
        end
            
        if it.cycle_pos[ni] < M
            println("Sequence position: ", pos, ", Phase: ", ni)
            fbits = extract_encoded_symbol(bitindex(it.seq, pos), encoded_data(it.seq))
            rbits = ~fbits & 0x03
            it.fkmer[ni] = ((it.fkmer[ni] << 2) | fbits) & kmer_mask(it)
            it.rkmer[ni] = (it.rkmer[ni] >> 2) | (UInt64(rbits) << firstoffset(it))
        end
    end
    ##

    fi = ifelse(fi == (N + 1), 0x01, fi + 0x01)
    
    fkmer = it.fkmer[fi]
    rkmer = it.rkmer[fi]
    
    outkmer = ifelse(fkmer < rkmer, fkmer, rkmer)
    
    return reinterpret(eltype(it), outkmer), (pos + 1, fi)  
    
end










function Base.iterate(it::BioSequences.EachSkipmerIterator{SK}, state) where {SK<:Skipmer}
    pos = state[1]
    fi = state[2]
    S = Int(BioSequences.span(SK))
    N = BioSequences.cycle_len(SK)
    M = BioSequences.bases_per_cycle(SK)
    while pos <= lastindex(it.seq)
        
        for ni in 1:N
            it.cycle_pos[ni] += 1
            if it.cycle_pos[ni] == S
                it.cycle_pos[ni] = 0
            end
            
            if it.cycle_pos[ni] < M
                println("Sequence position: ", pos, ", Phase: ", ni)
                fbits = BioSequences.twobitnucs[reinterpret(UInt8, it.seq[pos]) + 0x01]
                rbits = ~fbits & 0x03
                it.fkmer[ni] = ((it.fkmer[ni] << 2) | fbits) & BioSequences.kmer_mask(it)
                it.rkmer[ni] = (it.rkmer[ni] >> 2) | (UInt64(rbits) << BioSequences.firstoffset(it))
            end
        end
        
        # If we are at pos, the skip-mer that started at (pos - S) is now done. 
        if pos >= S
            if pos == S
                fi = 0x01
            else
                fi += 0x01
                if fi == (N + 1)
                    fi = 0x01
                end
            end
            if it.last_unknown[fi] + S <= pos
                if it.fkmer[fi] <= it.rkmer[fi]
                    return reinterpret(SK, it.fkmer[fi]), (pos + 1, fi)
                    #push!(skipmers, reinterpret(Skipmer{T, M, N, K}, it.fkmer[fi]))
                else
                    return reinterpret(SK, it.fkmer[fi]), (pos + 1, fi)
                    #push!(skipmers, reinterpret(Skipmer{T, M, N, K}, it.rkmer[fi]))
                end
            end
        end
        
        pos += 1
        
    end
    
    return nothing
end


