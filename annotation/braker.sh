#!/bin/bash
#SBATCH -J braker
#SBATCH -N 1
#SBATCH -n 48

export GENEMARK_PATH=gmes_linux_64_4
ProtHint-2.6.0/bin/prothint.py \
final.mask.fasta \
Vertebrata.fa \
--workdir prohint/ --threads 48

bam=$(ls *.bam | xargs -I {} printf {}"," | sed 's/,$//')

time braker.pl --cores 48 --species=assembly --genome=final.mask.fasta \
     --softmasking --bam=$bam \
     --gff3 --skipAllTraining --useexisting \
     --prot_seq=prohint/prothint_augustus.gff --prg=gth
