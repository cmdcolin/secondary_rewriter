#!/bin/bash

# write_secondaries.sh
# usage
# ./write_secondaries.sh <input.cram> <ref.fa> <output.cram> <nthreads default 4>
# e.g.
# ./write_secondaries.sh input.cram ref.fa output.cram 16

THR=${4:-4}

echo "Filtering secondary reads to sec.txt..."
samtools view -@$THR -h $1 -T $2 | cargo run -- --pass1 > sec.txt


echo "Applying seq and qual to secondary reads to sec2.txt..."
samtools view -@$THR -h $1 -T $2 | cargo run -- --pass2 --secondaries sec.txt > sec2.txt
echo Done pass 2

echo "Sorting sec2.txt by line number..."
LC_ALL=C sort -k1,1n sec2.txt > sec3.txt


echo "Re-writing file to $3..."
samtools view -@$THR -h $1 -T $2 | cargo run -- --pass3 --secondaries sec3.txt | samtools view -@$THR -T $2 - -o $3

echo "Finished!"

