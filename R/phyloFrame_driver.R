
main <- function(disease, results_path, new_batches){
  ## - load enhanced allele frequencies - ##
  eafs  <- load_EAF()
  ## - initialize given cancer - ##

  if(disease == "breast"){
    expression <- expression_breast
    samples_ancestry <- ancestry_breast
    tissue_network <- load_breast_network()
    cancer_type <- "brca"
    subtype1 <- "Basal"
    subtype2 <- "Luminal"
  }else if(disease == "thyroid"){
    #load("data/expression_thyroid.rda")
    expression <- expression_thyroid
    #load("data/ancestry_thyroid.rda")
    samples_ancestry <- ancestry_thyroid
    tissue_network <- load_thyroid_network()
    cancer_type <- "thca"
    subtype1 <- "M0"
    subtype2 <- "MX"
  }else if(disease == "uterine"){
    #load("data/expression_uterine.rda")
    expression <- expression_uterine
    #load("data/ancestry_uterine.rda")
    samples_ancestry <- ancestry_uterine
    tissue_network <- load_uterine_network()
    cancer_type <- "ucec"
    subtype1 <- "Endometrioid"
    subtype2 <- "Serous"
  }else{
    message("Please enter a valid disease, currently they are: 1. breast 2. thyroid 3. uterine")
    stop()
  }

  # en.mixture <- 1 #first classification run penalty for base signature
  # run.penalties <- 0 #second run for final signature
  # variable.genes <- 10000

  set_up_output_dirs(results_path)

  ## - code for if we want to recreate the batches for this run - ##
  #NOTE: if new batches are breasted they are put into {user defined directory}_samples directory in the data-raw section under a
  if(new_batches == TRUE){ ##TODO I THINK WE SHOULD CHANGE THIS TO BE IN THE RESULTS DIRECTORY
    samples_dir <- here::here("data-raw",paste0(tolower(results_path),"_","samples"))
    create_batches(samples_dir, cancer_type, expression, subtype1, subtype2)
  }else{
    samples_dir <- here::here("data-raw",paste0(tolower(cancer_type),"_","samples"))
  }
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

  # for each ancestry, for each batch of samples, run phyloFrame and test on other ancestries.

  curr_ancestry <- "eur"
  for(j in 1:eur_batch_num){
    set_up_model_dir(results_dir, curr_ancestry, j)

    output <- phyloFrame("eur", eur_expression, cancer_type, eafs, tissue_network, j, results_dir, FALSE)
    eur_bench_model <- output$benchmark
    eur_pf_model <- output$phyloFrame
    curr_samples <- output$samples
    benchmark_dir <- output$benchmark_dir
    phyloFrame_dir <- output$phyloFrame_dir

    ##- test on the ancestry that we trained on -##
    test_set <- eur_expression[!(rownames(eur_expression) %in% curr_samples$x),]
    if(nrow(test_set) != 0){
      test_model(eur_bench_model, test_set, "eur", benchmark_dir, j, subtype1, subtype2)
      test_model(eur_pf_model, test_set, "eur", phyloFrame_dir, j, subtype1, subtype2)
    }
    test_model(eur_bench_model, afr_expression, "afr", benchmark_dir, j, subtype1, subtype2)
    test_model(eur_pf_model, afr_expression, "afr", phyloFrame_dir, j, subtype1, subtype2)

    test_model(eur_bench_model, eas_expression, "eas", benchmark_dir, j, subtype1, subtype2)
    test_model(eur_pf_model, eas_expression, "eas", phyloFrame_dir, j, subtype1, subtype2)

    test_model(eur_bench_model, admix_expression, "admix", benchmark_dir, j, subtype1, subtype2)
    test_model(eur_pf_model, admix_expression, "admix", phyloFrame_dir, j, subtype1, subtype2)

    mixed_test_expression <- mixed_expression[!(rownames(mixed_expression) %in% curr_samples$x),]
    test_model(eur_bench_model, mixed_test_expression, "mixed", benchmark_dir, j, subtype1, subtype2)
    test_model(eur_pf_model, mixed_test_expression, "mixed", phyloFrame_dir, j, subtype1, subtype2)
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

    ##- test on the ancestry that we trained on -##
    test_set <- afr_expression[!(rownames(afr_expression) %in% curr_samples$x),]
    if(nrow(test_set) != 0){
      test_model(afr_bench_model, test_set, "afr", benchmark_dir, j, subtype1, subtype2)
      test_model(afr_pf_model, test_set, "afr", phyloFrame_dir, j, subtype1, subtype2)
    }
    test_model(afr_bench_model, eur_expression, "eur", benchmark_dir, j, subtype1, subtype2)
    test_model(afr_pf_model, eur_expression, "eur", phyloFrame_dir, j, subtype1, subtype2)

    test_model(afr_bench_model, eas_expression, "eas", benchmark_dir, j, subtype1, subtype2)
    test_model(afr_pf_model, eas_expression, "eas", phyloFrame_dir, j, subtype1, subtype2)

    test_model(afr_bench_model, admix_expression, "admix", benchmark_dir, j, subtype1, subtype2)
    test_model(afr_pf_model, admix_expression, "admix", phyloFrame_dir, j, subtype1, subtype2)

    mixed_test_expression <- mixed_expression[!(rownames(mixed_expression) %in% curr_samples$x),]
    test_model(afr_bench_model, mixed_test_expression, "mixed", benchmark_dir, j, subtype1, subtype2)
    test_model(afr_pf_model, mixed_test_expression, "mixed", phyloFrame_dir, j, subtype1, subtype2)
  }

  curr_ancestry <- "eas"
  # note that uterine does not have enough samples to make a model for east asian
  if(disease != "uterine"){
    for(j in 1:eas_batch_num){
      set_up_model_dir(results_dir, curr_ancestry, j)

      output <- phyloFrame("eas", eas_expression, cancer_type, eafs, tissue_network, j, results_dir, FALSE)
      eas_bench_model <- output$benchmark
      eas_pf_model <- output$phyloFrame
      curr_samples <- output$samples
      benchmark_dir <- output$benchmark_dir
      phyloFrame_dir <- output$phyloFrame_dir

      ##- test on the ancestry that we trained on -##
      test_set <- eas_expression[!(rownames(eas_expression) %in% curr_samples$x),]
      if(nrow(test_set) != 0){
        test_model(eas_bench_model, test_set, "eas", benchmark_dir, j, subtype1, subtype2)
        test_model(eas_pf_model, test_set, "eas", phyloFrame_dir, j, subtype1, subtype2)
      }
      test_model(eas_bench_model, eur_expression, "eur", benchmark_dir, j, subtype1, subtype2)
      test_model(eas_pf_model, eur_expression, "eur", phyloFrame_dir, j, subtype1, subtype2)

      test_model(eas_bench_model, afr_expression, "afr", benchmark_dir, j, subtype1, subtype2)
      test_model(eas_pf_model, afr_expression, "afr", phyloFrame_dir, j, subtype1, subtype2)

      test_model(eas_bench_model, admix_expression, "admix", benchmark_dir, j, subtype1, subtype2)
      test_model(eas_pf_model, admix_expression, "admix", phyloFrame_dir, j, subtype1, subtype2)

      mixed_test_expression <- mixed_expression[!(rownames(mixed_expression) %in% curr_samples$x),]
      test_model(eas_bench_model, mixed_test_expression, "mixed", benchmark_dir, j, subtype1, subtype2)
      test_model(eas_pf_model, mixed_test_expression, "mixed", phyloFrame_dir, j, subtype1, subtype2)
    }
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
    test_set <- admix_expression[!(rownames(admix_expression) %in% curr_samples$x),]
    if(nrow(test_set) != 0){
      test_model(admix_bench_model, test_set, "admix", benchmark_dir, j, subtype1, subtype2)
      test_model(admix_pf_model, test_set, "admix", phyloFrame_dir, j, subtype1, subtype2)
    }
    test_model(admix_bench_model, eur_expression, "eur", benchmark_dir, j, subtype1, subtype2)
    test_model(admix_pf_model, eur_expression, "eur", phyloFrame_dir, j, subtype1, subtype2)

    test_model(admix_bench_model, afr_expression, "afr", benchmark_dir, j, subtype1, subtype2)
    test_model(admix_pf_model, afr_expression, "afr", phyloFrame_dir, j, subtype1, subtype2)

    test_model(admix_bench_model, eas_expression, "eas", benchmark_dir, j, subtype1, subtype2)
    test_model(admix_pf_model, eas_expression, "eas", phyloFrame_dir, j, subtype1, subtype2)

    mixed_test_expression <- mixed_expression[!(rownames(mixed_expression) %in% curr_samples$x),]
    test_model(admix_bench_model, mixed_test_expression, "mixed", benchmark_dir, j, subtype1, subtype2)
    test_model(admix_pf_model, mixed_test_expression, "mixed", phyloFrame_dir, j, subtype1, subtype2)
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

    ##- test on the ancestry that we trained on -##
    test_set <- mixed_expression[!(rownames(mixed_expression) %in% curr_samples$x),]
    if(nrow(test_set) != 0){
      test_model(mixed_bench_model, test_set, "mixed", benchmark_dir, j, subtype1, subtype2)
      test_model(mixed_pf_model, test_set, "mixed", phyloFrame_dir, j, subtype1, subtype2)
    }

    #make sure current trained mixed samples are excluded from ancestry matrices
    eur_cut <- eur_expression[!(rownames(eur_expression) %in% curr_samples$x),]
    test_model(mixed_bench_model, eur_cut, "eur", benchmark_dir, j, subtype1, subtype2)
    test_model(mixed_pf_model, eur_cut, "eur", phyloFrame_dir, j, subtype1, subtype2)

    afr_cut <- afr_expression[!(rownames(afr_expression) %in% curr_samples$x),]
    test_model(mixed_bench_model, afr_cut, "afr", benchmark_dir, j, subtype1, subtype2)
    test_model(mixed_pf_model, afr_cut, "afr", phyloFrame_dir, j, subtype1, subtype2)

    eas_cut <- eas_expression[!(rownames(eas_expression) %in% curr_samples$x),]
    test_model(mixed_bench_model, eas_cut, "eas", benchmark_dir, j, subtype1, subtype2)
    test_model(mixed_pf_model, eas_cut, "eas", phyloFrame_dir, j, subtype1, subtype2)

    admix_cut <- admix_expression[!(rownames(admix_expression) %in% curr_samples$x),]
    test_model(mixed_bench_model, admix_cut, "admix", benchmark_dir, j, subtype1, subtype2)
    test_model(mixed_pf_model, admix_cut, "admix", phyloFrame_dir, j, subtype1, subtype2)
  }

  # plot results for each model
  if(disease != "uterine"){
    ancestries <- c("admix","afr","eas","eur","mixed")
  }else{
    ancestries <- c("admix","afr","eur","mixed")
  }

  lapply(ancestries, function(ancestry) {
    scat_plot(disease,results_path, ancestry, samples_dir)
  })
}
