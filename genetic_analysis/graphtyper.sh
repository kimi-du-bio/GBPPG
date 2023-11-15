#!/bin/bash
#SBATCH -J graphtyper
#SBATCH -N 1
#SBATCH -n 60

graphtyper genotype_sv ref.fasta \
SV.vcf.gz \
--sams=sample.bam.list \
--threads=60 \
--region_file=chr.txt
