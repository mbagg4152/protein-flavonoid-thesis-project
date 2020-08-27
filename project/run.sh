#!/bin/sh
#SBATCH -N 4
#SBATCH -J KeggPython
#SBATCH --mem 14000
#SBATCH -n 4
#SBATCH --mail-type=END,FAIL
module load python3
python3 kegg-prog.py