#!/bin/bash
#SBATCH --job-name=phyloFrame_run_single    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=leslie.smith1@ufl.edu     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=90gb
#SBATCH --time=120:05:00              # Time limit hrs:min:sec
#SBATCH --output=phyloFrame_single_%j.stdout   # Standard output and error log
#SBATCH --error=phyloFrame_single_%j.err # Error log


ml R

Rscript ./run_single.R thyroid TCGA_Thyroid_Gnomad4_cosmic_single_tester TCGA_eur_train_single_cos
Rscript ./run_single.R uterine TCGA_Uterine_Gnomad4_cosmic_single_tester TCGA_eur_train_single_cos
Rscript ./run_single.R breast TCGA_Breast_Gnomad4_cosmic_single_tester TCGA_eur_train_single_cos
