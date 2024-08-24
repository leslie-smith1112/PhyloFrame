
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)

disease <- args[1]
results_path <- args[2]

devtools::load_all()
## - main call for the batches version of phyloFrame on validation set - ##
breast_validation_main_batches(results_path)

## - main call for the single version of phyloFrame on validation set - ##
breast_validation_main_batches_main_single(paste0(results_path,"_single"), "TCGA")
