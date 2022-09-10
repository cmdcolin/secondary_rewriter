# secondary_rewriter

Some aligners such as minimap2 do not write the SEQ and QUAL fields to
secondary alignments making it hard to analyze them. This program adds these
back, referring to the primary alignment to get the SEQ and QUAL, and adding
them to the secondaries. This is a two- (or three- for sorted output) pass
program. The first pass collects the secondary alignments and the second adds
the SEQ and QUAL fields to the secondaries.

## Install

First install rust, probably with rustup https://rustup.rs/

Then

```
cargo install secondary_rewriter
```

## Usage

```
## First pass
samtools view yourfile.bam | secondary_rewriter --pass1 > secondaries.txt

## Second pass
samtools view -h yourfile.bam | secondary_rewriter --pass2 --secondaries secondaries.txt | samtools view -o out.bam

## Third pass, sort the BAM since all the new secondaries are added at the end of the file
samtools sort out.bam -o out.sorted.bam
```

The two-pass strategy works as follows

1. First pass: output ALL secondary alignments (reads with flag 256) to a separate file
2. Second pass: reading secondary alignments into memory, and then scan original SAM/BAM/CRAM to add SEQ and QUAL fields encountered during scan to the secondary alignments that are stored in a hashmap

It uses a two-pass strategy because otherwise it would have to effectively load
the entire SAM/BAM/CRAM into memory

## Note

My second rust project!
