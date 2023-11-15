#!/bin/bash
#SBATCH -J LACHESIS
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 200

bwa aln -t 200 polished.fasta \
hic.R1.fastq.gz \
> hic.R1.sai

bwa aln -t 200 polished.fasta \
hic.R2.fastq.gz \
> hic.R2.sai


bwa sampe polished.fasta \
hic.R1.sai \
hic.R2.sai \
hic.R1.fastq.g \
hic.R2.fastq.gz \
> hic.bwa_aln.sam


bash LACHESIS_software/PreprocessSAMs.sh

Lachesis LACHESIS.ini
