# secondary_rewriter

Some aligners such as minimap2 do not write the SEQ and QUAL fields to
secondary alignments (https://github.com/lh3/minimap2/issues/458) making it
hard to analyze them (for example, SNPs will not be visible in a genome browser
for secondary alignments). This program adds these back, referring to the
primary alignment to get the SEQ and QUAL, and adding them to the secondaries.

## Install

First install rust, probably with rustup https://rustup.rs/

Then

```
cargo install secondary_rewriter
```

## Usage

This small shell script automates the multi-step pipeline (supports BAM or CRAM)

```

#!/bin/bash

# write_secondaries.sh
# usage
# ./write_secondaries.sh <input.cram> <ref.fa> <output.cram> <nthreads default 4> <memory gigabytes, per-thread, default 1G>
# e.g.
# ./write_secondaries.sh input.cram ref.fa output.cram 16 2G

THR=${4:-4}
MEM=${5:-1G}


samtools view -@$THR -h $1 -T $2 -f256 > sec.txt
samtools view -@$THR -h $1 -T $2 -F256 | secondary_rewriter --secondaries sec.txt | samtools sort -@$THR -m $MEM - -o $3

```

## Two-pass strategy

The two-pass strategy works as follows

1. First pass: output ALL secondary alignments (reads with flag 256) to a
   external file
2. Second pass: read secondary alignments from external file into memory,
   and then scan original SAM/BAM/CRAM to add SEQ and QUAL fields on the
   primary alignments to the secondary alignments, and pipe to `samtools sort`
   (needed because all the new secondary reads with BAM/CRAM at the end)

This avoids loading the entire SAM/BAM/CRAM into memory, but does require a
re-sort with `samtools sort`. The sort is the most expensive part of the
process.

## Help

```

secondary_rewriter 0.1.6
Adds SEQ and QUAL fields to secondary alignments from the primary alignment

USAGE:
    secondary_rewriter [OPTIONS]

OPTIONS:
    -h, --help                         Print help information
    -s, --secondaries <SECONDARIES>
    -V, --version                      Print version information

```
