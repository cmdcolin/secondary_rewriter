#!/bin/bash

# write_secondaries.sh
# usage
# ./write_secondaries.sh <input.cram> <ref.fa> <output.cram> <nthreads default 4>
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

