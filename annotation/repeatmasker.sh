#!/bin/bash
#SBATCH -J repeatmasker
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --mem=100gb

module load RepeatMasker/4.1.2
RepeatMasker -pa 10 -species pig -gff -a -inv \
final.fasta -dir .
