# secondary_rewriter

Some aligners such as minimap2 do not write the SEQ and QUAL fields to
secondary alignments (https://github.com/lh3/minimap2/issues/458) making it
hard to analyze them (for example, SNPs will not be visible in a genome browser
for secondary alignments). This program adds these back, referring to the
primary alignment to get the SEQ and QUAL, and adding them to the secondaries.
This is a two- (or three- for sorted output) pass program. The first pass
collects the secondary alignments and the second adds the SEQ and QUAL fields
to the secondaries.

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

You can package this into a small bash script

```
#!/bin/bash

# write_secondaries.sh
# example usage
# ./write_secondaries.sh input.cram ref.fa output.cram


samtools view -@3 $1 -T $2 | secondary_rewriter --pass1 > secondaries.txt
samtools view -@3 -h $1 -T $2 | secondary_rewriter --pass2 --secondaries secondaries.txt | samtools view -@3 -T $2 - -o $3
```

The two-pass strategy works as follows

1. First pass: output ALL secondary alignments (reads with flag 256) to a separate file
2. Second pass: reading secondary alignments into memory, and then scan original SAM/BAM/CRAM to add SEQ and QUAL fields encountered during scan to the secondary alignments that are stored in a hashmap

It uses a two-pass strategy because otherwise it would have to effectively load
the entire SAM/BAM/CRAM into memory

## Help

```

% secondary_rewriter --help
secondary_rewriter 0.1.0
Adds SEQ and QUAL fields to secondary alignments which are often missing from minimap2

USAGE:
secondary_rewriter [OPTIONS]

OPTIONS:
-h, --help Print help information
--output-only-new-data
--pass1
--pass2
-s, --secondaries <SECONDARIES>
-V, --version Print version information

```

`--output-only-new-data` only outputs the secondary alignments with their new
SEQ/QUAL fields (and skips all other data). the default mode without this flag
passes all other alignments through stdout and adds the secondary alignments at
the end

## Note

My second rust project!

```

```
