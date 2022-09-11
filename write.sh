#!/bin/bash

# write_secondaries.sh
# usage
# ./write_secondaries.sh <input.cram> <ref.fa> <output.cram> <nthreads default 4>
# e.g.
# ./write_secondaries.sh input.cram ref.fa output.cram 16

THREADS=${4:-4}

samtools view -@$THREADS $1 -f 256 -T $2 | gzip -c > secondaries.txt.gz
samtools view -@$THREADS -h $1 -T $2 | secondary_rewriter --pass2 --secondaries secondaries.txt.gz | samtools sort -@$THREADS --reference $2 -o $3
