#!/bin/bash
#SBATCH --job-name=multi_validation_run    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<user email>     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --qos=<account>
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=95gb
#SBATCH --time=40:05:00              # Time limit hrs:min:sec
#SBATCH --output=multi_validation_%j.stdout   # Standard output and error log
#SBATCH --error=multi_validation_%j.err # Error log


ml R

# run_multi_study_validation.R calls function from R/multi_study_breast_validation.R
# Funtion call takes in the disease and output directory name - results put in results/$DEFINED_DIR
Rscript ./run_multi_study_validation.R breast Multi_Study_Validation
