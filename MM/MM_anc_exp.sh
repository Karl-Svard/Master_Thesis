#!/bin/bash -l

#SBATCH -A snic2021-22-197
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 09:00:00
#SBATCH -J MM_anc_exp
#SBATCH --mail-type=ALL
#SBATCH --mail-user karl.svard970@gmail.com


# load modules
module load WPS/4.1 # required for scipy to load
module load SciPy-bundle/2019.10-foss-2019b-Python-3.7.4
#module load python/2.7.15
#module load SciPy-bundle/2019.03-foss-2019a
module load MPFR/4.0.2-GCCcore-8.3.0
module load bioinfo-tools
module load msprime


# execute python script
python3 cc_anc_exp.py
