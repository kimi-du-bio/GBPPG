#!/bin/bash
#SBATCH -J assembly
#SBATCH -p cluster
#SBATCH -N 1
#SBATCH -n 200

module load python/2.7.18_falcon

# run it!
fc_run ./fc_run.cfg
