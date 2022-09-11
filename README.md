# secondary_rewriter

Some aligners such as minimap2 do not write the SEQ and QUAL fields to
secondary alignments (https://github.com/lh3/minimap2/issues/458) making it
hard to analyze them (for example, SNPs will not be visible in a genome browser
for secondary alignments). This program adds these back, referring to the
primary alignment to get the SEQ and QUAL, and adding them to the secondaries.
This is a two-pass program. The first pass collects the secondary alignments
and the second adds the SEQ and QUAL fields to the secondaries.

## Install

First install rust, probably with rustup https://rustup.rs/

Then

```
cargo install secondary_rewriter
```

## Usage

```
## First pass
samtools view -f256 yourfile.bam > secondaries.txt

## Second pass
samtools view -h yourfile.bam | secondary_rewriter --secondaries secondaries.txt | samtools sort -o out.bam

```

You can package this into a small bash script (supports CRAM) see
[write_secondaries.sh](write_secondaries.sh)

The two-pass strategy works as follows

1. First pass: output ALL secondary alignments (reads with flag 256) to a
   separate file (plaintext or gzip)
2. Second pass: reading secondary alignments into memory, and then scan
   original SAM/BAM/CRAM to add SEQ and QUAL fields encountered during scan to
   the secondary alignments that are stored in a hashmap

It uses a two-pass strategy because otherwise it would have to effectively load
the entire SAM/BAM/CRAM into memory

## Help

```

% secondary_rewriter --help

secondary_rewriter 0.1.2
Adds SEQ and QUAL fields to secondary alignments from the primary alignment

USAGE:
    secondary_rewriter [OPTIONS]

OPTIONS:
    -h, --help                         Print help information
        --output-only-new-data
    -s, --secondaries <SECONDARIES>
    -V, --version                      Print version information

```

`--output-only-new-data` only outputs the secondary alignments with their new
SEQ/QUAL fields (and skips all other data). the default mode without this flag
passes all other alignments through stdout and adds the secondary alignments at
the end

