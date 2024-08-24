
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandArgs(trailingOnly=FALSE))

## only to be used for the TCGA_runs
disease <- args[1]
results_path <- args[2]
training_set_name <- args[3]

devtools::load_all()
if(disease == "breast"){
  expression <- expression_breast
  samples_ancestry <- ancestry_breast
}else if(disease == "thyroid"){
  expression <- expression_thyroid
  samples_ancestry <- ancestry_thyroid
}else if(disease == "uterine"){
  expression <- expression_uterine
  samples_ancestry <- ancestry_uterine
}else{
  stop()
}
#main(disease, results_path, new_batches)
#define the expression matrix
train_samples <- samples_ancestry$patient[samples_ancestry$consensus_ancestry == "eur"]

single_expr_driver(expression,train_samples, disease, results_path, training_set_name)

