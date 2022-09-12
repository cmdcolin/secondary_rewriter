# secondary_rewriter

Some aligners such as minimap2 do not write the SEQ and QUAL fields to
secondary alignments (https://github.com/lh3/minimap2/issues/458) making it
hard to analyze them (for example, SNPs will not be visible in a genome browser
for secondary alignments). This program adds these back, referring to the
primary alignment to get the SEQ and QUAL, and adding them to the secondaries.
This is a three-pass program. The first pass collects the secondary alignments
into an external file the second adds the SEQ and QUAL fields to the external
file, and then the third pass inserts the secondary alignments in place.

## Install

First install rust, probably with rustup https://rustup.rs/

Then

```
cargo install secondary_rewriter
```

## Usage

```
samtools view file.bam | secondary_rewriter --pass1 > sec.txt
samtools view file.bam | secondary_rewriter --pass2 --secondaries sec.txt > sec2.txt
samtools view -h file.bam | secondary_rewriter --pass3 --secondaries sec2.txt | samtools view - -o out.bam
```

This small shell script automates this (supports CRAM)

```

#!/bin/bash

# write_secondaries.sh
# usage
# ./write_secondaries.sh <input.cram> <ref.fa> <output.cram> <nthreads default 4>
# e.g.
# ./write_secondaries.sh input.cram ref.fa output.cram 16

THR=${4:-4}

samtools view -@$THR -h $1 -T $2 | secondary_rewriter --pass1 > sec.txt
echo Done pass 1
samtools view -@$THR -h $1 -T $2 | secondary_rewriter --pass2 --secondaries sec.txt > sec2.txt
echo Done pass 2
samtools view -@$THR -h $1 -T $2 | secondary_rewriter --pass3 --secondaries sec2.txt.gz | samtools view -@$THR -T $2 - -o $3
echo Done pass 3


```

The three-pass strategy works as follows

1. First pass: output ALL secondary alignments (reads with flag 256) to a
   external file
2. Second pass: reading secondary alignments from external file into memory,
   and then scan original SAM/BAM/CRAM to add SEQ and QUAL fields on the
   primary alignments to the secondary alignments that are stored in a hashmap
3. Third pass: the secondary alignments with the new SEQ and QUAL fields are
   re-inserted in place into the SAM file, and re-encoded. By re-inserting them
   in-place, it avoids a full `samtools sort` on the output which is
   disk/memory/cpu intensive.

This seems laborious, but the three-pass strategy avoids loading the entire
SAM/BAM/CRAM into memory

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
