
#' Change subtype
#'
#' @param expression
#' @param split_criteria
#'
#' @return
#' @export
#'
#' @examples
subtype_prep <- function(expression){
  #check that there is a subtype column in the expression matrix - this is out binary prediction task
  if(!("subtype" %in% colnames(expression))){
    message("There muse be a subtype column with binary prediction outcome variables.")
    stop()
    #check that the subtype prediction has only 2 subtype - currently only have binary prediction avaoilable
  }else if(length(unique(expression$subtype)) != 2){
    message("There are too many outcome variables. Currently PhyloFrame only supports binary predictions.")
    stop()
    #if the subtypes are not 0 and 1, change them to be 0 and 1
    #here we assume that if one of the subtypes if not a 0 or 1 then the subtypes are not in 0 1 format
  }else if(!(expression$subtype[1] %in% c(0, 1))){ # if the binary variables are not 0 or 1 - make them
    df <- expression %>% dplyr::mutate(subtype_binary = as.numeric(factor(subtype)) - 1)
    #get rid of subtype column and change name of subtype binary
    expr_subtype <- df %>% dplyr::select(-subtype)
    names(expr_subtype)[names(expr_subtype) == "subtype_binary"] <- "subtype"
    expr_subtype$subtype <- as.factor(expr_subtype$subtype)
  }
  return(expr_subtype)
}


## main for running all as a single batch
#expression <- expression_breast

#must have the train samples you want and the network name
##TODO expression matrix must have samples as the rownames
#' Single batch phyloFrame driver
#'
#' @param in_expr Expression matrix with all samples (both train and test) must have samples as rownames
#' @param train_samples A list of the training samples that should be used to train the model
#' @param disease Disease to define the tissue network used
#' @param results_path Path to the directory you want results written to
#' @param training_set_name The name you want in the results file name - this is mostly to accomodate for the batch runs
#'
#' @return
#' @export
#'
#' @examples
single_expr_driver <- function(in_expr,train_samples, disease, results_path, training_set_name){
  #validation set will just be read in already correctlt defined

  eafs  <- load_EAF()
  message("Loaded Enhanced Allele Frequencies")
  if(disease == "breast"){
    tissue_network <- load_breast_network()
  }else if(disease == "thyroid"){
    tissue_network <- load_thyroid_network()
  }else if(disease == "uterine"){
    expression <- expression_uterine
    tissue_network <- load_uterine_network()
  }else{
    message("Please enter a valid disease, currently they are: 1. breast 2. thyroid 3. uterine")
    stop()
  }
  message("Loaded Tissue Network")
  ## notes:
  #need new set up dirs for a single expression
  expression <- subtype_prep(in_expr)
  #split into train and test data (here we soing this by using eur)
  #eur <- samples_ancestry[samples_ancestry$consensus_ancestry == "eur",]$patient
  train_expression <- expression[rownames(expression) %in% train_samples,]
  test_expression <- expression[!(rownames(expression) %in% train_samples),]


  set_up_output_dirs_single_run(results_path) #doesnt need the model or ancestry
  results_dir <- here::here("results",results_path)

  message(paste0("Results for all models will be in results/", results_path))

#the first parameter training_set_name doesnt matter I dont think. Its used for file and directory naming that we dont have in the single expr run
 # the trainin_set_name replaces i here - its for naming the output files
  output <- phyloFrame(NA, train_expression, cancer_type, eafs, tissue_network, training_set_name, results_dir, TRUE)
  bench_model <- output$benchmark
  pf_model <- output$phyloFrame
  curr_samples <- output$samples
  benchmark_dir <- output$benchmark_dir
  phyloFrame_dir <- output$phyloFrame_dir #here::here(results_dir, "phyloFrame")

  ##- test on the ancestry that we trained on -##
  test_model(bench_model, test_expression, "test_set", benchmark_dir, training_set_name, "0", "1")
  test_model(pf_model, test_expression, "test_set", phyloFrame_dir, training_set_name, "0", "1")

  scat_plot_single(disease,results_dir, training_set_name, "test_set")
}
#  write.table(metrics,paste0(out_dir,"/",file_prefix,"_metrics.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

