#!/bin/bash
#SBATCH -J smartie_sv
#SBATCH -N 1
#SBATCH -n 30

snakemake -s Snakefile -w 50 -p -k -j 30
