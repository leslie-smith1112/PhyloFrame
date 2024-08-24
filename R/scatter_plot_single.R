
#' Scatter Plot Metrics for Single batch
#'
#' @param disease Single string value: the disease we are currently plotting. Only used in plot title.
#' @param results_dir Single string value: the name of the directory that contains the other run results.
#' @param prim_ancestry Single string value: the ancestry name for the model whose results we are currently plotting.
#' @param training_set_name The string name used to name the model results
#' @param set_name another name used in the file naming of the results
#'
#' @return NULL. Function writes plots to file in "results/{results_dir}/plots"
#' @export
#'
#' @examples scat_plot("brca", "brca_run", "afr")
scat_plot_single <- function(disease, results_dir, training_set_name, test_set_name){
#model_TCGA_eur_train_metrics.tsv
  message("Starting creation of metric scatter plots.")
  benchmark_results <- readr::read_tsv(here::here(results_dir,"benchmark",paste0(test_set_name,"_",training_set_name,"_metrics.tsv")))
  phyloFrame_results <- readr::read_tsv(here::here(results_dir,"phyloFrame",paste0(test_set_name,"_",training_set_name,"_metrics.tsv")))

  metric_dat <- data.frame("metric" = c("roc_auc","accuracy","sens","spec","precision","recall","f_meas","kap","mcc"),
                          "phyloFrame" = c(phyloFrame_results$.estimate[phyloFrame_results$.metric == "roc_auc"],
                                           phyloFrame_results$.estimate[phyloFrame_results$.metric == "accuracy"],
                                           phyloFrame_results$.estimate[phyloFrame_results$.metric == "sens"],
                                           phyloFrame_results$.estimate[phyloFrame_results$.metric == "spec"],
                                           phyloFrame_results$.estimate[phyloFrame_results$.metric == "precision"],
                                           phyloFrame_results$.estimate[phyloFrame_results$.metric == "recall"],
                                           phyloFrame_results$.estimate[phyloFrame_results$.metric == "f_meas"],
                                           phyloFrame_results$.estimate[phyloFrame_results$.metric == "kap"],
                                           phyloFrame_results$.estimate[phyloFrame_results$.metric == "mcc"]),
                           "benchmark" = c(benchmark_results$.estimate[benchmark_results$.metric == "roc_auc"],
                                           benchmark_results$.estimate[benchmark_results$.metric == "accuracy"],
                                           benchmark_results$.estimate[benchmark_results$.metric == "sens"],
                                           benchmark_results$.estimate[benchmark_results$.metric == "spec"],
                                           benchmark_results$.estimate[benchmark_results$.metric == "precision"],
                                           benchmark_results$.estimate[benchmark_results$.metric == "recall"],
                                           benchmark_results$.estimate[benchmark_results$.metric == "f_meas"],
                                           benchmark_results$.estimate[benchmark_results$.metric == "kap"],
                                           benchmark_results$.estimate[benchmark_results$.metric == "mcc"]))


  p <- ggplot2::ggplot(metric_dat, ggplot2::aes(x=benchmark, y=phyloFrame, color=metric)) + ggplot2::geom_point(alpha = 8/10) + ggplot2::ggtitle(paste0(toupper(training_set_name)," MODEL METRICS FOR ",toupper(disease))) +
    ggplot2::scale_color_manual(values = c("roc_auc" = "#eb0f1e",
                                           "accuracy"="#0fd3eb",
                                           "sens"="#25bb86",
                                           "spec" = "#e534c3",
                                           "precision" = "#7b3ae4",
                                           "recall" = "#b1129a",
                                           "f_meas" = "#f6ee1e",
                                           "kap" = "#eb6f07",
                                           "mcc" = "#04139a")) + ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::geom_abline() + ggplot2::geom_point(ggplot2::aes(size = 2))
  #p
  if(dir.exists(paste0(results_dir,"/plots")) == FALSE){
    dir.create(paste0(results_dir,"/plots"))
  }
  png(paste0(results_dir,"/plots","/",training_set_name, "_model_metrics.png"),width = 800, height = 600)
  print(p)
  dev.off()
}



get_master_dat_sub_sah_single <- function(model_num, prim_ancestry,test_ancestry, results_dir){
  # we need this function to track the ancestry of the files we read in


  return(all)
}


#' Scatter Plot for Validation Single Runs
#'
#' @param disease
#' @param results_dir
#'
#' @return
#' @export
#'
#' @examples
scat_plot_sub_sah_single <- function(disease, results_dir){
  #for each model of each ancesry we want to plot their metrics. So here we check to see how many models the
  #current ancestry has, if it is 1 or 0 they do not have any other batches to test on so we dont include
  #in the ancestries we will read in
  test_ancestries <- c("African_American","African_Ethiopian","African_Ghanaian")
  bm_AA <-
  met_files <- paste0(test_ancestries,"_metrics.tsv")

  benchmark_files <- paste0(here::here("results", results_dir, "benchmark"),"/" ,met_files)
  phyloFrame_files <- paste0(here::here("results", results_dir,"phyloFrame"),"/", met_files)

  bm_files <- lapply(benchmark_files, readr::read_tsv)
  bm_files[[1]]$Ancestry <- "African_American"
  bm_files[[2]]$Ancestry <- "African_Ethiopian"
  bm_files[[3]]$Ancestry <- "African_Ghanaian"

  pf_files <- lapply(phyloFrame_files, readr::read_tsv)
  pf_files[[1]]$Ancestry <- "African_American"
  pf_files[[2]]$Ancestry <- "African_Ethiopian"
  pf_files[[3]]$Ancestry <- "African_Ghanaian"

  bm_scores <- do.call("rbind", bm_files)
  pf_scores <- do.call("rbind", pf_files)

  bm_scores$model <- "benchmark"

  pf_scores$model <- "phyloFrame"


  to_plot <- data.frame("ancestry" = bm_scores$Ancestry, "benchmark_recall" = bm_scores$Recall, "phyloFrame_recall" = pf_scores$Recall,
                    "benchmark_precision" = bm_scores$Precision, "phyloFrame_precision" = pf_scores$Precision)
  ## - for every test ancestry get the results for each model run of the primary ancestry and bind the data frames into a single dataframe - ##

  if(dir.exists(here::here("results",results_dir,"plots")) == FALSE){
    dir.create(here::here("results",results_dir,"plots"))
  }


  ## recall
  p <- ggplot2::ggplot(to_plot, ggplot2::aes(x=benchmark_recall, y=phyloFrame_recall, color=ancestry)) + ggplot2::geom_point(alpha = 8/10) + ggplot2::ggtitle(paste0("SUB SAHARAN SINGLE MODEL RECALL FOR ",toupper(disease))) +
    ggplot2::scale_color_manual(values = c("African_American"="#86ecf1",
                                           "African_Ethiopian"="#341ef2",
                                           "African_Ghanaian" = "#1e95f2")) + ggplot2::xlim(0,1) + ggplot2::ylim(0,1) + ggplot2::geom_abline() + ggplot2::geom_point(ggplot2::aes(size = 2))
  #p
  png(paste0(here::here("results",results_dir,"plots"),"/","sub_sah_validation_model_recall.png"),width = 800, height = 600)
  print(p)
  dev.off()


}


