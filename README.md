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

This small shell script automates the multi-step pipeline (supports BAM or CRAM)

```

#!/bin/bash

# write_secondaries.sh
# usage
# ./write_secondaries.sh <input.bam/cram> <ref.fa> <output.bam/cram> <nthreads default 4>
# e.g.
# ./write_secondaries.sh input.cram ref.fa output.cram 16

THR=${4:-4}

echo "Filtering secondary reads to sec.txt..."
samtools view -@$THR -h $1 -T $2 | secondary_rewriter --pass1 > sec.txt


echo "Applying seq and qual to secondary reads, and sorting by line, outputting to sec2.txt..."
samtools view -@$THR -h $1 -T $2 | secondary_rewriter --pass2 --secondaries sec.txt | LC_ALL=C sort -k1,1n --parallel=$THR > sec2.txt

echo "Creating $3 with seq and qual on secondary reads..."
samtools view -@$THR -h $1 -T $2 | secondary_rewriter --pass3 --secondaries sec2.txt | samtools view -@$THR -T $2 - -o $3

echo "Finished!"
```

The three-pass strategy works as follows

1. First pass: output ALL secondary alignments (reads with flag 256) to a
   external file
2. Second pass: reading secondary alignments from external file into memory,
   and then scan original SAM/BAM/CRAM to add SEQ and QUAL fields on the
   primary alignments to the secondary alignments, outputting a new file with
   `line_num_of_secondary_read SEQ QUAL` and this is then sorted using the
   external program sort
3. Third pass: the original BAM is read line by line in parallel with the
   `line_num_of_secondary_read SEQ QUAL` file, re-inserted in place into the
   SAM file, and re-encoded. By re-inserting them in-place, it avoids a full
   `samtools sort` on the output which may be disk/cpu intensive. Not that step
   2 is still fairly disk/cpu intensive

This seems laborious, but the strategy avoids loading the entire SAM/BAM/CRAM
into memory and a full samtools sort of the entire file

## Help

```

Adds SEQ and QUAL fields to secondary alignments from the primary alignment

USAGE:
    secondary_rewriter [OPTIONS]

OPTIONS:
    -h, --help                         Print help information
        --pass1
        --pass2
        --pass3
    -s, --secondaries <SECONDARIES>
    -V, --version                      Print version information

```
