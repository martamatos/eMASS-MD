#!/bin/bash -l

#SBATCH -A m2009
#SBATCH -t 5:00:00
#SBATCH -p regular

module load g09

cd /global/homes/n/nmih/marta_temp/S7P 

g09 S7P_ideal.gau > S7P_ideal.out
