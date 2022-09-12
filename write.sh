#!/bin/bash

# write_secondaries.sh
# usage
# ./write_secondaries.sh <input.cram> <ref.fa> <output.cram> <nthreads default 4> <memory gigabytes default 16G>
# e.g.
# ./write_secondaries.sh input.cram ref.fa output.cram 16 32G

THR=${4:-4}
MEM=${5:-16G}


samtools view -@$THR -h $1 -T $2 -f256 > sec.txt
samtools view -@$THR -h $1 -T $2 -F256 | secondary_rewriter --secondaries sec.txt | samtools sort -@$THR -m $MEM - -o $3

