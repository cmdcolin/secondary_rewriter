#!/bin/bash

# write_secondaries.sh
# example usage
# ./write_secondaries.sh input.cram ref.fa output.cram


samtools view -@16 $1 -f 256 -T $2 | gzip -c > secondaries.txt.gz
samtools view -@16 -h $1 -T $2 | secondary_rewriter --pass2 --secondaries secondaries.txt.gz | samtools sort -@16 --reference $2 -o $3
