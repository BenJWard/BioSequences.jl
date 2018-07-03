using Documenter, BioSequences, BioSymbols

makedocs(
    format = :html,
    sitename = "BioSequences.jl",
    pages = [
        "Home" => "index.md",

        "User Manual" => [
            "Biological Symbols" => "user_manual/symbols.md",
            "Sequence types" => [
                "Overview" => "user_manual/sequences/biosequence.md",
                "General Sequence" => "user_manual/sequences/generalseq.md",
                "Reference Sequences" => "user_manual/sequences/refseq.md",
                "Nucleic acid k-mers" => "user_manual/sequences/kmer.md"
            ],
            "Indexing sequences" => "user/manual/indexing.md",
            "IO" => [
                "FASTA formatted files" => "user_manual/io/fasta.md",
                "FASTQ formatted files" => "user_manual/io/fastq.md",
                "2bit formatted files" => "user_manual/io/twobit.md"
            ],
            "Searching" => "user_manual/search.md",
            "Sequence Composition" => "user_manual/composition.md",
            "Demultiplexing" => "user_manual/demultiplexer.md",
        ],

        "Developer Notes" => [
            "Biological Symbols" => "dev_docs/symbols.md",
            "Biological Sequences" => [
                "BioSequence" => "dev_docs/sequences/biosequence.md"
            ],
        ]
    ],
    authors = "Ben J. Ward, D.C.Jones, Kenta Sato, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/BioSequences.jl.git",
    julia = "0.6",
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing
)
