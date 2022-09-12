#!/bin/bash

# write_secondaries.sh
# usage
# ./write_secondaries.sh <input.cram> <ref.fa> <output.cram> <nthreads default 4>
# e.g.
# ./write_secondaries.sh input.cram ref.fa output.cram 16

THR=${4:-4}

samtools view -@$THR -h $1 -T $2 | cargo run -- --pass1 > sec.txt
echo Done pass 1
samtools view -@$THR -h $1 -T $2 | cargo run -- --pass2 --secondaries sec.txt > sec2.txt
echo Done pass 2
samtools view -@$THR -h $1 -T $2 | cargo run -- --pass3 --secondaries sec2.txt | samtools view -@$THR -T $2 - -o $3
echo Done pass 3

