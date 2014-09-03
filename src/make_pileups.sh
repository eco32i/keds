#!/bin/bash

# Align, sort, index, generate .pileup
# Processes all *.gz files in ../data/trimmed directory

REF=/path/to/ref/index
FASTAREF=/path/to/fasta
INPUT=../data/trimmed
OUTPUT=../results

BARCODES="AGGAAT ACATCG TGGTCA GGACGG CGAAAC GATCTG"

for bc in $BARCODES
do
    fastaref=${FASTAREF}_${bc}.fasta
    refbase=${REF}_${bc}
    bowtie2-build $fastaref $refbase
    for fq in `ls ${INPUT} | grep ${bc}`
    do
        base=$(basename $fq)
        bambase=`echo $base | cut -d'.' -f1`
        gunzip -c ${INPUT}/${fq} | bowtie2 -p 30 -N 1 -x $refbase -U - | samtools view -bhS - > ${OUTPUT}/${bambase}.bam
        samtools sort ${OUTPUT}/${bambase}.bam ${OUTPUT}/${bambase}_sorted
        samtools index ${OUTPUT}/${bambase}_sorted.bam
        samtools mpileup -d 5000000 -f $fastaref ${OUTPUT}/${bambase}_sorted.bam > ${OUTPUT}/${bambase}.pileup &
    done
    rm ../ref/dHSR1/tmp/dHSR*
done
    
