#!/bin/bash
#SBATCH --job-name=phyloFrame_main_run    # Job name
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<user email>     # Where to send mail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=96gb
#SBATCH --time=120:05:00              # Time limit hrs:min:sec
#SBATCH --output=phyloFrame_%j.stdout   # Standard output and error log
#SBATCH --error=phyloFrame_%j.err # Error log


ml R
#run.R calls main function call from R/phyloFrame_driver.R
Rscript ./run.R breast TCGA_Breast_Gnomad4_tester FALSE 
Rscript ./run.R thyroid TCGA_Thyroid_Gnomad4_tester FALSE
Rscript ./run.R uterine TCGA_Uterine_Gnomad4_tester FALSE
#Rscript ./run.R breast TCGA_Breast_Gnomad4_range250v30genes FALSE
