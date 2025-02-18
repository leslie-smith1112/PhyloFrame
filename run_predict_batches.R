#splitting predicted samples into batches
devtools::load_all()

args <- commandArgs(trailingOnly = TRUE)
brca <- args[1] #should just be a directory name: EX: "TCGA_Thyroid_Gnomad4_corrected"
ucec <- args[2]
thca <- args[3]
# Get the batch number for each sample
prep_batches(brca)
prep_batches(ucec)
prep_batches(thca)

predict_batches(brca, "breast")
predict_batches(ucec, "thyroid")
predict_batches(thca, "uterine")
