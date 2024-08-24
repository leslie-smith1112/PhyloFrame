
#' Multi Study Validation batch main
#'
#' @param results_path
#'
#' @return
#' @export
#'
#' @examples
breast_validation_main_batches <- function(results_path){

  ## - validation_raw_expr is original name of expression matrix raw
  ## - was subbed with multi_val_expr_exclude_GSE142102 to try without study predicted wrong often by both benchmark and phyloFrame

  ## - load enhanced allele frequencies - ##
  eafs  <- load_EAF()
  ## - initialize multi-study breast validation test set, this is TCGA has been preprocessed to have common genes as other datset- ##
  # vali without the sub saharan:  breast_validation_training
  expression <- tcga_validation_training
  samples_ancestry <- ancestry_breast
  tissue_network <- load_breast_network()
  cancer_type <- "brca"
  subtype1 <- "Basal"
  subtype2 <- "Luminal"

  set_up_output_dirs(results_path)

  ## - code for if we want to recreate the batches for this run - ##
  #NOTE: if new batches are breasted they are put into {user defined directory}_samples directory in the data-raw section under a

  samples_dir <- here::here("data-raw",paste0(tolower(cancer_type),"_","samples"))

  results_dir <- here::here("results",results_path,"model_runs")
  message(paste0("Results for all models will be in results/", results_path,"/model_runs"))
  # for each ancestry
  # -get the number of batches
  # -read in all mixed samples so that we can test on mixed. Since mixed samples are interwoven with other batches
  # -seperate ancestry samples out so we can test

  eur <- samples_ancestry[samples_ancestry$consensus_ancestry == "eur",]$patient
  afr <- samples_ancestry[samples_ancestry$consensus_ancestry == "afr",]$patient
  eas <- samples_ancestry[samples_ancestry$consensus_ancestry == "eas",]$patient
  admix <- samples_ancestry[samples_ancestry$consensus_ancestry %in% c("eas_admix","afr_admix","eur_admix","admix"),]$patient
  #read in mixed samples - not all samples are used so we need to read in the batches to get samples
  mixed_sample_files <- list.files(here::here("data-raw",paste0(cancer_type,"_samples"),"mixed"), full.names = TRUE)
  #paste path the sample files
  mixed_list <- lapply(mixed_sample_files, readr::read_tsv)
  mixed_combined <- do.call(rbind, mixed_list)
  mixed <- mixed_combined$x

  # get the number of batches (equated to the number of sample batch files we have)
  eur_batch_num <- length(list.files(here::here(samples_dir,"eur")))
  afr_batch_num <- length(list.files(here::here(samples_dir,"afr")))
  eas_batch_num <- length(list.files(here::here(samples_dir,"eas")))
  admix_batch_num <- length(list.files(here::here(samples_dir,"admix")))
  mixed_batch_num <- length(list.files(here::here(samples_dir,"mixed")))

  ## - get ancestry specific expression matrices - ##
  eur_expression <- expression[rownames(expression) %in% eur,]
  afr_expression <- expression[rownames(expression) %in% afr,]
  eas_expression <- expression[rownames(expression) %in% eas,]
  admix_expression <- expression[rownames(expression) %in% admix,]
  mixed_expression <- expression[rownames(expression) %in% mixed,]

  ## - batch corrected matrix: multi_study_breast_val_expr
  # for each ancestry, for each batch of samples, run phyloFrame and test on other ancestries
  curr_ancestry <- "eur"
  for(j in 1:eur_batch_num){
    set_up_model_dir(results_dir, curr_ancestry, j)

    output <- phyloFrame("eur", eur_expression, cancer_type, eafs, tissue_network, j, results_dir, FALSE)
    eur_bench_model <- output$benchmark
    eur_pf_model <- output$phyloFrame
    curr_samples <- output$samples
    benchmark_dir <- output$benchmark_dir
    phyloFrame_dir <- output$phyloFrame_dir

    ##- test on the validation set -##
    test_model(eur_bench_model, all_val_expr, "multi" ,benchmark_dir, "Multi_study_validation", subtype1, subtype2)
    test_model(eur_pf_model,  all_val_expr, "multi" ,phyloFrame_dir, "Multi_study_validation", subtype1, subtype2)

  }


  curr_ancestry <- "afr"
  for(j in 1:afr_batch_num){
    set_up_model_dir(results_dir, curr_ancestry, j)

    output <- phyloFrame("afr", afr_expression, cancer_type, eafs, tissue_network, j, results_dir, FALSE)
    afr_bench_model <- output$benchmark
    afr_pf_model <- output$phyloFrame
    curr_samples <- output$samples
    benchmark_dir <- output$benchmark_dir
    phyloFrame_dir <- output$phyloFrame_dir

    test_model(afr_bench_model,  all_val_expr, "multi" ,benchmark_dir, "Multi_study_validation", subtype1, subtype2)
    test_model(afr_pf_model,  all_val_expr, "multi" ,phyloFrame_dir, "Multi_study_validation", subtype1, subtype2)

  }

  curr_ancestry <- "eas"
  # note that uterine does not have enough samples to make a model for east asian
  for(j in 1:eas_batch_num){
    set_up_model_dir(results_dir, curr_ancestry, j)

    output <- phyloFrame("eas", eas_expression, cancer_type, eafs, tissue_network, j, results_dir, FALSE)
    eas_bench_model <- output$benchmark
    eas_pf_model <- output$phyloFrame
    curr_samples <- output$samples
    benchmark_dir <- output$benchmark_dir
    phyloFrame_dir <- output$phyloFrame_dir

    ##- test on the ancestry that we trained on -##
    test_model(eas_bench_model, all_val_expr, "multi"  ,benchmark_dir, "Multi_study_validation", subtype1, subtype2)
    test_model(eas_pf_model, all_val_expr, "multi"  ,phyloFrame_dir, "Multi_study_validation", subtype1, subtype2)

  }


  curr_ancestry <- "admix"
  for(j in 1:admix_batch_num){
    set_up_model_dir(results_dir, curr_ancestry, j)

    output <- phyloFrame("admix", admix_expression, cancer_type, eafs, tissue_network, j, results_dir, FALSE)
    admix_bench_model <- output$benchmark
    admix_pf_model <- output$phyloFrame
    curr_samples <- output$samples
    benchmark_dir <- output$benchmark_dir
    phyloFrame_dir <- output$phyloFrame_dir

    ##- test on the ancestry that we trained on -##
    test_model(admix_bench_model, all_val_expr, "multi"  ,benchmark_dir, "Multi_study_validation", subtype1, subtype2)
    test_model(admix_pf_model, all_val_expr, "multi"  ,phyloFrame_dir, "Multi_study_validation", subtype1, subtype2)

  }

  curr_ancestry <- "mixed"
  for(j in 1:mixed_batch_num){
    set_up_model_dir(results_dir, curr_ancestry, j)

    output <- phyloFrame("mixed", mixed_expression, cancer_type, eafs, tissue_network, j, results_dir, FALSE)
    mixed_bench_model <- output$benchmark
    mixed_pf_model <- output$phyloFrame
    curr_samples <- output$samples
    benchmark_dir <- output$benchmark_dir
    phyloFrame_dir <- output$phyloFrame_dir

    test_model(mixed_bench_model, all_val_expr, "multi"  ,benchmark_dir, "Multi_study_validation", subtype1, subtype2)
    test_model(mixed_pf_model, all_val_expr, "multi"  ,phyloFrame_dir, "Multi_study_validation", subtype1, subtype2)

  }

  # plot results for each model
  ## - ancestry models we are testing - ##
  ancestries <- c("admix","afr","eas","eur","mixed")
  lapply(ancestries, function(ancestry) {
    scat_plot_multi(disease,results_path, ancestry, samples_dir)
  })
}

####################################### Multi study Validation Set Single Batch #######################################

#' Multi Study Validation single
#'
#' @param results_path
#' @param training_set_name
#'
#' @return
#' @export
#'
#' @examples
breast_validation_main_batches_main_single <- function(results_path, training_set_name){
  ## - validation_raw_expr is original name of expression matrix raw
  ## - was subbed with multi_val_expr_exclude_GSE142102 to try without study predicted wrong often by both benchmark and phyloFrame

  #validation set will just be read in already correctlt defined
  eafs  <- load_EAF()
  message("Loaded Enhanced Allele Frequencies")
  expression <- tcga_validation_training
  samples_ancestry <- ancestry_breast
  ## - for only european training: ##
  #train_samples <- samples_ancestry$patient[samples_ancestry$consensus_ancestry == "eur"]

  train_samples <- rownames(breast_validation_training)
  tissue_network <- load_breast_network()
  cancer_type <- "brca"
  subtype1 <- "Basal"
  subtype2 <- "Luminal"
  message("Loaded Tissue Network")

  train_expression <- expression[rownames(expression) %in% train_samples,]
  #test_expression <- expression[!(rownames(expression) %in% train_samples),]


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

  test_model(bench_model, all_val_expr, "multi"  ,benchmark_dir, "Multi_study_validation", subtype1, subtype2)
  test_model(pf_model, all_val_expr, "multi"  ,phyloFrame_dir, "Multi_study_validation", subtype1, subtype2)


  scat_plot_single(disease,results_path)
}

