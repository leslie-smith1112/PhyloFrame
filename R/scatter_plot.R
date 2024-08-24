
#' Read in ancestry scores for a given ancestry model
#'
#' @param model_num Character List: a character list of 1 the total number of models that were trained for the primary ancestry, so we know how many results metrics to expect
#' @param prim_ancestry Single string value: the ancestry of the model we are plotting metrics for, used to read in files
#' @param test_ancestry Single string value: the test ancestry; the data the primary ancestry model was tested on used to read in files
#'
#' @return A data frame with the mcc and roc for phyloFrame and the benchmark for the given ancestry
#' @export
#'
#' @examples get_master_dat(list("1","2","3"), "afr", "admixed")
get_master_dat <- function(model_num, prim_ancestry,test_ancestry, results_dir){
  # we need this function to track the ancestry of the files we read in
  model_num <- as.character(model_num)
  #combinations <- expand.grid(anc = test_ancestry, model = model_num)
  # model_list <- apply(combinations, 1, function(x) paste0("model_",x['model'],"/",test_ancestry,"_", x['model'], "_metrics.tsv"))
  # head(model_list)
  ## tempo for old results
  #model_list <- apply(combinations, 1, function(x) paste0("model_",x['model'],"/",test_ancestry,"_",x['model'], "_metrics.tsv"))
  model_list <- paste0("model_",model_num,"/",test_ancestry,"_",model_num, "_metrics.tsv")
  #model_list <- paste0("model_", model_num, "/multi_Multi_study_validation_metrics.tsv")
  #head(model_list)
  #model_1admixed_metrics.tsv
  benchmark_files <- paste0(here::here("results", results_dir, "model_runs", "benchmark",prim_ancestry),"/" ,model_list)
  phyloFrame_files <- paste0(here::here("results", results_dir, "model_runs","phyloFrame",prim_ancestry),"/", model_list)

  bm_files <- lapply(benchmark_files, readr::read_tsv)
  pf_files <- lapply(phyloFrame_files, readr::read_tsv)

  bm_scores <- do.call("rbind", bm_files)
  pf_scores <- do.call("rbind", pf_files)

 results_bm_roc <- bm_scores[bm_scores$.metric %in% c("roc_auc"),]
 results_bm_mcc <- bm_scores[bm_scores$.metric %in% c("mcc"),]
 results_bm_precision <- bm_scores[bm_scores$.metric %in% c("precision"),]
 results_bm_recall <- bm_scores[bm_scores$.metric %in% c("recall"),]
 results_bm_roc$ancestry <- test_ancestry
 results_bm_roc$model <- "benchmark"

 results_pf_roc <- pf_scores[pf_scores$.metric %in% c("roc_auc"),]
 results_pf_mcc <- pf_scores[pf_scores$.metric %in% c("mcc"),]
 results_pf_precision <- pf_scores[pf_scores$.metric %in% c("precision"),]
 results_pf_recall <- pf_scores[pf_scores$.metric %in% c("recall"),]
 results_pf_roc$ancestry <- test_ancestry
 results_pf_roc$model <- "phyloFrame"

 all <- data.frame("ancestry" = results_bm_roc$ancestry, "benchmark_roc" = results_bm_roc$.estimate,
                   "phyloFrame_roc" = results_pf_roc$.estimate, "benchmark_mcc" = results_bm_mcc$.estimate, "phyloFrame_mcc" = results_pf_mcc$.estimate,
                   "benchmark_precision" = results_bm_precision$.estimate, "phyloFrame_precision" = results_pf_precision$.estimate,
                   "benchmark_recall" = results_bm_recall$.estimate, "phyloFrame_recall" = results_pf_recall$.estimate)
 return(all)
}

#' Scatter Plot Metrics
#'
#' @param disease Single string value: the disease we are currently plotting. Only used in plot title.
#' @param results_dir Single string value: the name of the directory that contains the other run results.
#' @param prim_ancestry Single string value: the ancestry name for the model whose results we are currently plotting.
#'
#' @return NULL. Function writes plots to file in "results/{results_dir}/plots"
#' @export
#'
#' @examples scat_plot("brca", "brca_run", "afr")
scat_plot <- function(disease, results_dir, prim_ancestry, samples_dir){
  #for each model of each ancesry we want to plot their metrics. So here we check to see how many models the
  #current ancestry has, if it is 1 or 0 they do not have any other batches to test on so we dont include
  #in the ancestries we will read in
  test_ancestries <- c("admix","afr","eas","eur","mixed")
  model_num <- 1:length(list.files(here::here(samples_dir,prim_ancestry)))

  if(length(model_num) <= 1){
    test_ancestries <- test_ancestries[test_ancestries != prim_ancestry]
  }
  model_num <- as.character(model_num)

  results <- lapply(test_ancestries, function(test_ancestry) {
    get_master_dat(model_num, prim_ancestry, test_ancestry, results_dir)
  })

  to_plot  <- do.call(rbind, results)
## roc
  p <- ggplot2::ggplot(to_plot, ggplot2::aes(x=benchmark_roc, y=phyloFrame_roc, color=ancestry)) + ggplot2::geom_point(alpha = 8/10) + ggplot2::ggtitle(paste0(toupper(prim_ancestry)," MODEL ROC FOR ",toupper(disease))) +
    ggplot2::scale_color_manual(values = c("afr" = "#306cc3",
                                  "admix"="#eaa512",
                                  "eas"="#4f8826",
                                  "eur" = "#c33039",
                                  "mixed" = "#822cb6")) + ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::geom_abline() + ggplot2::geom_point(ggplot2::aes(size = 2))
  #p
  if(dir.exists(here::here("results",results_dir,"plots")) == FALSE){
    dir.create(here::here("results",results_dir,"plots"))
  }
  png(paste0(here::here("results",results_dir,"plots"),"/",prim_ancestry, "_model_roc.png"),width = 800, height = 600)
  print(p)
  dev.off()
## mcc
  p <- ggplot2::ggplot(to_plot, ggplot2::aes(x=benchmark_mcc, y=phyloFrame_mcc, color=ancestry)) + ggplot2::geom_point(alpha = 8/10) + ggplot2::ggtitle(paste0(toupper(prim_ancestry)," MODEL MCC FOR ",toupper(disease))) +
    ggplot2::scale_color_manual(values = c("afr" = "#306cc3",
                                  "admix"="#eaa512",
                                  "eas"="#4f8826",
                                  "eur" = "#c33039",
                                  "mixed" = "#822cb6")) + ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::geom_abline() + ggplot2::geom_point(ggplot2::aes(size = 2))
  #p
  png(paste0(here::here("results",results_dir,"plots"),"/",prim_ancestry, "_model_mcc.png"),width = 800, height = 600)
  print(p)
  dev.off()

  ## precision
  p <- ggplot2::ggplot(to_plot, ggplot2::aes(x=benchmark_precision, y=phyloFrame_precision, color=ancestry)) + ggplot2::geom_point(alpha = 8/10) + ggplot2::ggtitle(paste0(toupper(prim_ancestry)," MODEL PRECISION FOR ",toupper(disease))) +
    ggplot2::scale_color_manual(values = c("afr" = "#306cc3",
                                           "admix"="#eaa512",
                                           "eas"="#4f8826",
                                           "eur" = "#c33039",
                                           "mixed" = "#822cb6")) + ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::geom_abline() + ggplot2::geom_point(ggplot2::aes(size = 2))
  #p
  png(paste0(here::here("results",results_dir,"plots"),"/",prim_ancestry, "_model_precision.png"),width = 800, height = 600)
  print(p)
  dev.off()

  ## recall
  p <- ggplot2::ggplot(to_plot, ggplot2::aes(x=benchmark_recall, y=phyloFrame_recall, color=ancestry)) + ggplot2::geom_point(alpha = 8/10) + ggplot2::ggtitle(paste0(toupper(prim_ancestry)," MODEL RECALL FOR ",toupper(disease))) +
    ggplot2::scale_color_manual(values = c("afr" = "#306cc3",
                                           "admix"="#eaa512",
                                           "eas"="#4f8826",
                                           "eur" = "#c33039",
                                           "mixed" = "#822cb6")) + ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::geom_abline() + ggplot2::geom_point(ggplot2::aes(size = 2))
  #p
  png(paste0(here::here("results",results_dir,"plots"),"/",prim_ancestry, "_model_recall.png"),width = 800, height = 600)
  print(p)
  dev.off()


}


################################### - Sub-Saharan Validation Plots - ########################################
#' Function to get master metrics dataframe for plotting
#'
#' @param model_num
#' @param prim_ancestry
#' @param test_ancestry
#' @param results_dir
#'
#' @return
#' @export
#'
#' @examples
get_master_dat_sub_sah <- function(model_num,prim_ancestry,test_ancestry, results_dir){
  # we need this function to track the ancestry of the files we read in
  model_list <- paste0("model_",model_num,"/",test_ancestry,"_metrics.tsv")
  head(model_list)
  #model_1admixed_metrics.tsv
  benchmark_files <- paste0(here::here("results", results_dir, "model_runs", "benchmark",prim_ancestry),"/" ,model_list)
  phyloFrame_files <- paste0(here::here("results", results_dir, "model_runs","phyloFrame",prim_ancestry),"/", model_list)

  bm_files <- lapply(benchmark_files, readr::read_tsv)
  pf_files <- lapply(phyloFrame_files, readr::read_tsv)

  bm_scores <- do.call("rbind", bm_files)
  pf_scores <- do.call("rbind", pf_files)

  bm_scores$ancestry <- test_ancestry
  bm_scores$model <- "benchmark"

  pf_scores$ancestry <- test_ancestry
  pf_scores$model <- "phyloFrame"


  all <- data.frame("ancestry" = bm_scores$ancestry, "benchmark_recall" = bm_scores$Recall, "phyloFrame_recall" = pf_scores$Recall,
                    "benchmark_precision" = bm_scores$Precision, "phyloFrame_precision" = pf_scores$Precision)
  return(all)
}

# write.table(to.write, paste0(directory,"/", out_file,"_metrics.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

#' Scatter Plot Main for SubSaharan Breast Cancer Validation Set
#'
#' @param disease
#' @param results_dir
#' @param prim_ancestry
#' @param samples_dir
#'
#' @return
#' @export
#'
#' @examples
scat_plot_sub_sah <- function(disease, results_dir, prim_ancestry, samples_dir){
  #for each model of each ancesry we want to plot their metrics. So here we check to see how many models the
  #current ancestry has, if it is 1 or 0 they do not have any other batches to test on so we dont include
  #in the ancestries we will read in
  test_ancestries <- c("African_American","African_Ethiopian","African_Ghanaian")
  model_num <- 1:length(list.files(here::here(samples_dir,prim_ancestry)))
  model_num <- as.character(model_num)

  ## - for every test ancestry get the results for each model run of the primary ancestry and bind the data frames into a single dataframe - ##
  results <- lapply(test_ancestries, function(test_ancestry) {
    get_master_dat_sub_sah(model_num, prim_ancestry, test_ancestry, results_dir)
  })

  to_plot  <- do.call(rbind, results)
  if(dir.exists(here::here("results",results_dir,"plots")) == FALSE){
    dir.create(here::here("results",results_dir,"plots"))
  }

  # ## precision
  # p <- ggplot2::ggplot(to_plot, ggplot2::aes(x=benchmark_precision, y=phyloFrame_precision, color=ancestry)) + ggplot2::geom_point(alpha = 8/10) + ggplot2::ggtitle(paste0(toupper(prim_ancestry)," MODEL PRECISION FOR ",toupper(disease))) +
  #   ggplot2::scale_color_manual(values = c("African_American"="#86ecf1",
  #                                          "African_Ethiopian"="#341ef2",
  #                                          "African_Ghanaian" = "#1e95f2")) + ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::geom_abline() + ggplot2::geom_point(ggplot2::aes(size = 2))
  # #p
  # png(paste0(here::here("results",results_dir,"plots"),"/",prim_ancestry, "_model_precision.png"),width = 800, height = 600)
  # print(p)
  # dev.off()

  ## recall
  p <- ggplot2::ggplot(to_plot, ggplot2::aes(x=benchmark_recall, y=phyloFrame_recall, color=ancestry)) + ggplot2::geom_point(alpha = 8/10) + ggplot2::ggtitle(paste0(toupper(prim_ancestry)," MODEL RECALL FOR ",toupper(disease))) +
    ggplot2::scale_color_manual(values = c("African_American"="#86ecf1",
                                           "African_Ethiopian"="#341ef2",
                                           "African_Ghanaian" = "#1e95f2")) + ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::geom_abline() + ggplot2::geom_point(ggplot2::aes(size = 2))
  #p
  png(paste0(here::here("results",results_dir,"plots"),"/",prim_ancestry, "_model_recall.png"),width = 800, height = 600)
  print(p)
  dev.off()


}

########################################### MULTI STUDY VALIDATION SCATTER PLOT ##########################################
scat_plot_multi <- function(disease, results_dir, prim_ancestry, samples_dir){
  #for each model of each ancesry we want to plot their metrics. So here we check to see how many models the
  #current ancestry has, if it is 1 or 0 they do not have any other batches to test on so we dont include
  #in the ancestries we will read in
  test_ancestries <- c("multi")
  model_num <- 1:length(list.files(here::here(samples_dir,prim_ancestry)))
  #model_num <- "Multi_study_validation"
  if(length(model_num) <= 1){
    test_ancestries <- test_ancestries[test_ancestries != prim_ancestry]
  }
  model_num <- as.character(model_num)

  results <- lapply(test_ancestries, function(test_ancestry) {
    get_master_dat(model_num, prim_ancestry, test_ancestry, results_dir)
  })

  to_plot  <- do.call(rbind, results)
  #all_dat <- tidyr::gather(to_plot, model_metric, value, benchmark_roc:phyloFrame_recall)

  roc <- data.frame( "Metric" = rep("ROC", nrow(to_plot)),"Benchmark" = to_plot$benchmark_roc, "PhyloFrame" = to_plot$phyloFrame_roc)
  mcc <- data.frame( "Metric" = rep("MCC", nrow(to_plot)),"Benchmark" = to_plot$benchmark_mcc, "PhyloFrame" = to_plot$phyloFrame_mcc)
  recall <- data.frame( "Metric" = rep("RECALL", nrow(to_plot)),"Benchmark" = to_plot$benchmark_recall, "PhyloFrame" = to_plot$phyloFrame_recall)
  precision <- data.frame( "Metric" = rep("PRECISION", nrow(to_plot)),"Benchmark" = to_plot$benchmark_precision, "PhyloFrame" = to_plot$phyloFrame_precision)

  all_dat <- rbind(roc,mcc,recall,precision)

  #all <- data.frame("
                    # "benchmark_precision" = results_bm_precision$.estimate, "phyloFrame_precision" = results_pf_precision$.estimate,
                    # "benchmark_recall" = results_bm_recall$.estimate, "phyloFrame_recall" = results_pf_recall$.estimate)
  ## roc
  p <- ggplot2::ggplot(all_dat, ggplot2::aes(x=Benchmark, y=PhyloFrame, color=Metric)) + ggplot2::geom_point(alpha = 8/10) + ggplot2::ggtitle(paste0(toupper(prim_ancestry)," MODEL METRICS FOR ",toupper(disease))) +
    ggplot2::scale_color_manual(values = c("ROC" = "#f5b154",
                                           "MCC"="#7e0ca3",
                                           "RECALL"="#dfda3f",
                                           "PRECISION" = "#12e9f0")) + ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::geom_abline() + ggplot2::geom_point(ggplot2::aes(size = 2))
  #p
  if(dir.exists(here::here("results",results_dir,"plots")) == FALSE){
    dir.create(here::here("results",results_dir,"plots"))
  }
  png(paste0(here::here("results",results_dir,"plots"),"/",prim_ancestry, "_model_roc.png"),width = 800, height = 600)
  print(p)
  dev.off()


}

