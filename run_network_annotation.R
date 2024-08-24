
args <- commandArgs(trailingOnly = TRUE)
tissue_network <- args[1]

devtools::load_all()
message("Starting network annotation.")
annotate_network(tissue_network)
message("Annotation complete.")

