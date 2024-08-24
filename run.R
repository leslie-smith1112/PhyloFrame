
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)
# trailingOnly=TRUE means that only your arguments are returned, check:
# print(commandArgs(trailingOnly=FALSE))

disease <- args[1]
results_path <- args[2]
new_batches <- args[3]

devtools::load_all()
main(disease, results_path, new_batches)
