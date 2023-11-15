#!/bin/bash
#SBATCH -J gapclose
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 200

TGS-GapCloser-release_v1.0.1/TGS-GapCloser.sh \
--scaff  LACHESIS.fasta \
--reads  pacbio_reads.fasta \
--output final \
--racon  racon \
--tgstype pb \
--thread 200
