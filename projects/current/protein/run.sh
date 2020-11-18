#!/bin/bash
#SBATCH -p nodes
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 15000

module load python3
module load mvapich2-2.2-psm/gcc
cd || exit
cd "/home/student/hort4152/drsasa/current/protein" || exit
srun python3 run_dr_sasa.py
