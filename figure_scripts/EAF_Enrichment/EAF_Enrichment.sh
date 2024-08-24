#!/bin/bash
#SBATCH --job-name=RidgeEAFs    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=leslie.smith1@ufl.edu     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --qos=kgraim-b
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=350gb
#SBATCH --time=03:05:00              # Time limit hrs:min:sec
#SBATCH --output=ridgeEafs%j.stdout   # Standard output and error log
#SBATCH --error=ridgeEafs%j.err # Error log


ml R

#Rscript --max-ppsize=500000 ./explore_Gnomad.R
Rscript ./figure_scripts/Ridge_EAFs.R

