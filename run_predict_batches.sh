#!/bin/bash
#SBATCH --job-name=predict_batches    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<user email>     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3gb
#SBATCH --time=120:05:00              # Time limit hrs:min:sec
#SBATCH --output=predict_batches_%j.stdout   # Standard output and error log
#SBATCH --error=predict_batches_%j.err # Error log


ml R
Rscript ./run_predict_batches.R
