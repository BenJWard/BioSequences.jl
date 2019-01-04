struct PairedSeqDatastore
    filename::String
    readsize::UInt64
    chunksize::UInt64
    readpos_offset::UInt64
    size::UInt64
    stream::IOStream
end

# PairedSeqDatastore Header
# =========================
# 
# | Field             | Value  | Type   |
# |:-----------------:|:------:|:------:|
# | Magic number      | 0x0B56 | UInt16 |
# | Datastore type    | 0x0001 | UInt16 |
# | Version number    | 0x0001 | UInt16 |
# | Maximum read size | N/A    | UInt64 |
# | Chunk size        | N/A    | UInt64 |
# | Number of Reads   | N/A    | UInt64 |

const PSD_Version = 0x0001
const PSD_Type    = 0x0001

# Junk may get written in the unused space of the readpair chunk.

function PairedSeqDatastore(rdrx, rdry, outname::String, minsize::UInt64, maxsize::UInt64)
    # Create and allocate the sequence and record objects.
    lread = FASTQ.Record()
    rread = FASTQ.Record()
    lseq = GeneralSequence{DNAAlphabet{2}}(maxsize)
    rseq = GeneralSequence{DNAAlphabet{2}}(maxsize)
    
    chunksize::UInt64 = seq_data_len(DNAAlphabet{2}, maxsize)
    fd = open(outname, "w")
    # Write magic no, datastore type, version no, read size, and chunk size.
    sizepos = write(fd, 0x0B56, PSD_Type, PSD_Version, maxsize, chunksize)
    # Write space for size variable (or number of read pairs).
    readpos = write(fd, UInt64(0))
    
    @info string("readspos: ", readpos)
    @info string("sizepos: ", sizepos)
    
    pairs = discarded = truncated = 0
    
    while !eof(rdrx) && !eof(rdry)
        # Read in the two records.
        read!(rdrx, lread)
        read!(rdry, rread)
        
        llen = FASTQ.seqlen(lread)
        rlen = FASTQ.seqlen(rread)
        # If either read is too short, discard them both.
        if llen < minsize || rlen < minsize
            discarded += 1
            continue
        end
        
        if llen > maxsize
            truncated += 1
            ln = maxsize
        else
            ln = llen
        end
        
        if rlen > maxsize
            truncated += 1
            rn = maxsize
        else
            rn = rlen
        end
        
        pairs += 1
        
        # Copy FASTQ records to sequence variables, thus encoding
        # them in 2 bit format.
        copyto!(lseq, 1, lread, 1, ln)
        copyto!(rseq, 1, rread, 1, rn)
        
        # Write the two sizes reads one after the other: 
        # For each sequence, write the size, and then the datachunk.
        
        write(fd, ln)
        write(fd, lseq.data)
        
        write(fd, rn)
        write(fd, rseq.data)
    end
    
    nreads::UInt64 = pairs * 2
    
    seek(fd, sizepos)
    write(fd, nreads)
    
    close(fd)
    
    @info "Done writing paired sequences to datastore."
    @info string(discarded, " read pairs were discarded due to a too short sequence.")
    @info string(truncated, " reads were truncated to ", maxsize, " base pairs.")
    @info string("Created paired sequence datastore with ", pairs, " sequence pairs.")
    
    stream = open(outname, "r+")
    
    return PairedSeqDatastore(outname, maxsize, chunksize, readpos, nreads, stream)
end

function Base.open(::Type{PairedSeqDatastore}, filename::String)
    fd = open(filename, "r")
    magic = read(fd, UInt16)
    dstype = read(fd, UInt16)
    version = read(fd, UInt16)
    
    @assert magic == 0x0B56
    @assert dstype == PSD_Type
    @assert version == PSD_Version
    
    readsize = read(fd, UInt64)
    chunksize = read(fd, UInt64)
    nreads = read(fd, UInt64)
    readpos_offset = position(fd)
    
    @info magic
    @info dstype
    @info version
    @info readsize
    @info chunksize
    @info nreads
    @info readpos_offset
    
    return PairedSeqDatastore(filename, readsize, chunksize, readpos_offset, nreads, fd)
end