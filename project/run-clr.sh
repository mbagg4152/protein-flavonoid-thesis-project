#!/bin/sh
#SBATCH -N 4
#SBATCH -J KeggPython
#SBATCH --mem 1400
#SBATCH -n 4
#SBATCH --mail-type=END,FAIL
module load python3
rm -rf data
python3 kegg-prog.py