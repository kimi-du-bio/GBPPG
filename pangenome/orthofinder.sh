#!/bin/bash
#SBATCH -J orthofinder
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 80

conda activate orthofinder
orthofinder -f protein/ \
-M msa -S diamond -t 80 -a 80 -T iqtree
