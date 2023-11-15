#!/bin/bash
#SBATCH -J hisat2
#SBATCH -N 1
#SBATCH -n 40

tissue=
hisat2-build final.mask final.mask.fasta
for i in $tissue
do
	mkdir ./$i
	cd ./$i
	hisat2 -p 40 \
	-t -x final.mask \
	-1 transcript.$i.R1.fq.gz \
	-2 transcript.$i.R2.fq.gz \
	-S transcript.$i.sam

	samtools view -bhS \
	-t final.mask.fasta.fai \
	-@ 40 -o transcript.$i.bam \
 	transcript.$i.sam

 	samtools sort -@ 40 \
 	-o transcript.$i.sort.bam \
 	transcript.$i.bam
done
