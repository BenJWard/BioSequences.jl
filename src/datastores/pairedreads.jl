struct PairedSeqDatastore
    filename::String
    readsize::UInt64
    chunksize::UInt64
    readpos_offset::UInt64
    size::UInt64
    stream::IOStream
end

function PairedSeqDatastore(rdrx, rdry, outname::String, minsize::UInt64, maxsize::UInt64)
    # Create and allocate the sequence and record objects.
    lread = FASTQ.Record()
    rread = FASTQ.Record()
    lseq = GeneralSequence{DNAAlphabet{2}}(maxsize)
    rseq = GeneralSequence{DNAAlphabet{2}}(maxsize)
    
    chunksize::UInt64 = seq_data_len(DNAAlphabet{2}, maxsize)
    fd = open(outname, "w")
    # Write magic no, version no, datastore type, read size, and chunk size.
    sizepos = write(fd, 0x0B56, 0x0001, 0x0001, maxsize, chunksize)
    # Write space for size variable (or number of read pairs).
    readpos = write(fd, UInt64(0))
    
    @info string("readspos: ", readspos)
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
        
        # Write the two sizes of the sequences one after the other.
        write(fd, ln)
        write(fd, rn)
        
        # Write the data chunk to file.
        write(fd, lseq.data)
        write(fd, rseq.data)
    end
    
    nreads::UInt64 = pairs * 2
    
    seek(fd, sizepos)
    write(fd, nreads)
    
    close(fd)
    
    @info "Done writing paired sequences to datastore."
    @info string(discarded, " read pairs were discarded due to a too short sequence.")
    @info string(truncated, " reads were truncated to ", maxsize, " base pairs.")
    @info string("Created paired sequence datastore with ", size, " sequence pairs.")
    
    stream = open(outname, "r+")
    
    return PairedSeqDatastore(outname, maxsize, chunksize, readpos, nreads, stream)
end