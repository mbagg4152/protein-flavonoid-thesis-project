#!/bin/sh
#SBATCH -N 1
#SBATCH -J kegg-py-script
#SBATCH --mem 20000
#SBATCH -n 1
#SBATCH --mail-type=END,FAIL

cd
cd "/home/student/hort4152/project" || exit
module load python3
python3 kegg-prog.py
