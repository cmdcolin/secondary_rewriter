#!/bin/bash
set -v

# write_secondaries.sh
# usage
# ./write_secondaries.sh <input.cram> <ref.fa> <output.cram> <nthreads default 4> <memory allocated per-thread, default 1G>
# e.g.
# ./write_secondaries.sh input.cram ref.fa output.cram 16 2G 

THR=${4:-4}

## this is per-thread
MEM=${5:-1G}

samtools view -@$THR -h $1 -T $2 -f256 > sec.txt
samtools view -@$THR -h $1 -T $2 -F256 | cargo run -- --secondaries sec.txt | samtools sort -@$THR -m $MEM - -o $3

