#!/bin/bash
#SBATCH -J syri
#SBATCH -N 1
#SBATCH -n 30

minimap2 -t 30 \
-ax asm5 --eqx ref.fa \
query.fa \
> align.sam

syri -c align.sam \
-r ref.fa  \
-q query.fa -k -F S --nc 30
