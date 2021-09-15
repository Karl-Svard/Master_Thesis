#!/bin/bash -l
 
#SBATCH -A snic2021-22-197
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 07:40:00
#SBATCH -J TTbeta_anc_exp
#SBATCH --mail-type=ALL
#SBATCH --mail-user karl.svard970@gmail.com

# Load modules
module load WPS/4.1 # required for scipy to load
module load SciPy-bundle/2019.10-foss-2019b-Python-3.7.4
module load bioinfo-tools
module load msprime


# Commands
python3 anc_expansion_sim.py 
