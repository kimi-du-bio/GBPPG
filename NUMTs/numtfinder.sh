#!/bin/bash
#SBATCH -J numtfinder
#SBATCH -N 1
#SBATCH -n 30

module load R/4.0.2

python SLiMSuite-1.11.0/dev/numtfinder.py \
seqin=genome.fasta \
mtdna=Duroc_MT.fasta \
circle=T forks=30
