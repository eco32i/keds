#!/bin/bash

# Trim adapters from 3' ends of the reads
ADAPTER=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
cd ../data

if [ ! -d trimmed ]; then
    mkdir trimmed
fi

# plus channel samples
for fq in `ls *plus*`
do
    cutadapt -a $ADAPTER -o trimmed/${fq} $fq
done

# minus channel samples
for fq in `ls *minus*`:
do
    cutadapt -a $ADAPTER -o trimmed/${fq} $fq
done
cd -
