

##### seperate into batches to get predictions #####
#directory <- "TCGA_Breast_Gnomad4_corrected"


#' Seperate sample batches and get metrics
#'
#' @param directory  Directory for the results of all the sample predictions
#'
#' @return
#' @export
#'
#' @examples
predict_batches <- function(directory, disease){
  # ver <- unique(dat$version) ## for both phyloFrame and benchmark
  # for(i in ver){ #for both phyloFrame and benchmark
  #   subset_ver <- dat[dat$version == i,]
  dat <- readr::read_tsv(here::here("results",directory, "model_runs", "sample_batches_and_predictions.tsv"))
  out_dir <- here::here("results",directory, "model_runs","batch_predictions")
  dir.create(out_dir)
  master_dat <- data.frame(matrix(nrow = 0, ncol = 22))

  training_dat <- unique(dat$model_train_data)
  for(j in training_dat){#for every type of training data (ancestry)
    #break it down by data the model was trained on - ex first we will take all european trained models
    subset_train <- dat[dat$model_train_data == j,]
    #grab the number of models for this specific training data
    model_num<- unique(subset_train$model_num)
    for(h in model_num){
      #break it down to a specific model's predictions
      subset_num <- subset_train[subset_train$model_num == h,]
      #from there we will break it down by ancestry test data
      anc_batch <- unique(subset_num$sample_ancestry_batch)
      dir.create(paste0(out_dir,"/",j,"/",h), recursive = TRUE)
      for(k in anc_batch){
        subset_anc_batch <- subset_num[subset_num$sample_ancestry_batch == k,]
        #from there we will break it down by specific batch predictions for that test ancestry
        batch_num <- unique(subset_anc_batch$sample_batch_num)
        for(n in batch_num){
          final_subset <- subset_anc_batch[subset_anc_batch$sample_batch_num == n,]
          ## now we should have a specifci batch
          pf_final <- final_subset[final_subset$version == "phyloFrame",]
          bm_final <- final_subset[final_subset$version == "benchmark",]
          if(disease == "breast"){
            pf_mets <- get_metrics_brca(pf_final)
            bm_mets <- get_metrics_brca(bm_final)
          }else if(disease == "thyroid"){
            pf_mets <- get_metrics_thca(pf_final)
            bm_mets <- get_metrics_thca(bm_final)
          }else if(disease == "uterine"){
            pf_mets <- get_metrics_ucec(pf_final)
            bm_mets <- get_metrics_ucec(bm_final)
          }else{
            message("Error parsing disease prediction file.")
          }
          colnames(pf_mets) <- c("pf_auc", "pf_sensitivity","pf_specificity",
                                 "pf_accuracy", "pf_precision", "pf_recall",
                                 "pf_fmeas","pf_kapp", "pf_mcc")
          colnames(bm_mets) <- c("bm_auc", "bm_sensitivity","bm_specificity",
                                 "bm_accuracy", "bm_precision", "bm_recall",
                                 "bm_fmeas","bm_kapp", "bm_mcc")

          all_metrics <- cbind(pf_mets, bm_mets)
          all_metrics$train_data <- j
          all_metrics$model_num <- h
          all_metrics$ancestry_batch <- k
          all_metrics$batch_num <- n

          out_file <- paste0(out_dir,"/",j,"/",h, "/",k,"_",n,".tsv")
          #dir.create(paste0(out_dir,"/",i,"/",j,"/",h))

          write.table(all_metrics, out_file, sep = "\t", col.names = TRUE, row.names = FALSE)
          master_dat[nrow(master_dat) + 1, ] <- all_metrics
        }
      }
    }
  }
  #}
  colnames(master_dat) <- c("pf_auc", "pf_sensitivity","pf_specificity",
                            "pf_accuracy", "pf_precision", "pf_recall",
                            "pf_fmeas","pf_kapp", "pf_mcc","bm_auc", "bm_sensitivity","bm_specificity",
                            "bm_accuracy", "bm_precision", "bm_recall",
                            "bm_fmeas","bm_kapp", "bm_mcc", "train_data","model_num","ancestry_batch","batch_num")
  write.table(master_dat,paste0(out_dir,"/master_predictions.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE )
}


#' Breast Batch specific metrics
#'
#' @param results these are the sample predictions from all runs of a disease
#'
#' @return return the metrics for benchmark and phyloFrame for a batch
#' @export
#'
#' @examples
get_metrics_brca <- function(results){
  results$subtype <- factor(x = results$subtype, levels = c("Basal", "Luminal"))
  results$.pred_class <- factor(x = results$.pred_class, levels = c("Basal", "Luminal"))
  result_metrics <- dplyr::bind_cols(results$subtype, results$.pred_subtype1, results$.pred_subtype2, results$.pred_class)
  colnames(result_metrics) <- c("subtype", ".pred_subtype1", ".pred_subtype2", ".pred_class")
  result_metrics$subtype <- factor(x = result_metrics$subtype, levels = c("Basal", "Luminal"))

  auc <- yardstick::roc_auc(result_metrics, subtype, .pred_subtype1)
  senss <- yardstick::sens(results,  subtype, .pred_class)
  specc <- yardstick::spec(results,  subtype,  .pred_class)
  acc <- yardstick::accuracy(results,  subtype, .pred_class)
  prec <- yardstick::precision(results,  subtype, .pred_class)
  re <- yardstick::recall(results, subtype, .pred_class)
  f <- yardstick::f_meas(results, subtype, .pred_class)
  kapp <- yardstick::kap(results,  subtype, .pred_class)
  mccc <- yardstick::mcc(results, subtype, .pred_class)
  #metrics <- rbind(auc, acc, senss, specc, prec, re, f, kapp, mccc)
  metrics <- data.frame("auc" = auc$.estimate, "senss" = senss$.estimate, "specc" = specc$.estimate,
                        "acc" = acc$.estimate, "prec" = prec$.estimate, "re" = re$.estimate, "f" = f$.estimate,
                        "kapp" = kapp$.estimate, "mcc" = mccc$.estimate)
  return(metrics)
}
# result_metrics <-  bind_cols(test_expression$subtype, results[,4], results[,5], results$.pred_class)
# colnames(result_metrics) <- c("subtype", paste0(".pred_",subtype_1), paste0(".pred_",subtype_2), ".pred_class")
# confusion <- yardstick::conf_mat(results, truth = subtype,estimate = .pred_class)
# conf_matrix <- as.data.frame(confusion$table)
# write.table(conf_matrix, paste0(out_dir,"/",file_prefix,"_confusion_matrix.tsv"), sep  = "\t")


#' Thyroid Batch specific metrics
#'
#' @param results these are the sample predictions from all runs of a disease
#'
#' @return return the metrics for benchmark and phyloFrame for a batch
#' @export
#'
#' @examples
get_metrics_thca <- function(results){
  results$subtype <- factor(x = results$subtype, levels = c("M0", "MX"))
  results$.pred_class <- factor(x = results$.pred_class, levels = c("M0", "MX"))
  result_metrics <- dplyr::bind_cols(results$subtype, results$.pred_subtype1, results$.pred_subtype2, results$.pred_class)
  colnames(result_metrics) <- c("subtype", ".pred_subtype1", ".pred_subtype2", ".pred_class")
  result_metrics$subtype <- factor(x = result_metrics$subtype, levels = c("M0", "MX"))

  auc <- yardstick::roc_auc(result_metrics, subtype, .pred_subtype1)
  senss <- yardstick::sens(results,  subtype, .pred_class)
  specc <- yardstick::spec(results,  subtype,  .pred_class)
  acc <- yardstick::accuracy(results,  subtype, .pred_class)
  prec <- yardstick::precision(results,  subtype, .pred_class)
  re <- yardstick::recall(results, subtype, .pred_class)
  f <- yardstick::f_meas(results, subtype, .pred_class)
  kapp <- yardstick::kap(results,  subtype, .pred_class)
  mccc <- yardstick::mcc(results, subtype, .pred_class)
  #metrics <- rbind(auc, acc, senss, specc, prec, re, f, kapp, mccc)
  metrics <- data.frame("auc" = auc$.estimate, "senss" = senss$.estimate, "specc" = specc$.estimate,
                        "acc" = acc$.estimate, "prec" = prec$.estimate, "re" = re$.estimate, "f" = f$.estimate,
                        "kapp" = kapp$.estimate, "mcc" = mccc$.estimate)
  return(metrics)
}

#' Uterine Batch specific metrics
#'
#' @param results these are the sample predictions from all runs of a disease
#'
#' @return return the metrics for benchmark and phyloFrame for a batch
#' @export
#'
#' @examples
get_metrics_ucec <- function(results){
  results$subtype <- factor(x = results$subtype, levels = c("Endometrioid", "Serous"))
  results$.pred_class <- factor(x = results$.pred_class, levels = c("Endometrioid", "Serous"))
  result_metrics <- dplyr::bind_cols(results$subtype, results$.pred_subtype1, results$.pred_subtype2, results$.pred_class)
  colnames(result_metrics) <- c("subtype", ".pred_subtype1", ".pred_subtype2", ".pred_class")
  result_metrics$subtype <- factor(x = result_metrics$subtype, levels = c("Endometrioid", "Serous"))

  auc <- yardstick::roc_auc(result_metrics, subtype, .pred_subtype1)
  senss <- yardstick::sens(results,  subtype, .pred_class)
  specc <- yardstick::spec(results,  subtype,  .pred_class)
  acc <- yardstick::accuracy(results,  subtype, .pred_class)
  prec <- yardstick::precision(results,  subtype, .pred_class)
  re <- yardstick::recall(results, subtype, .pred_class)
  f <- yardstick::f_meas(results, subtype, .pred_class)
  kapp <- yardstick::kap(results,  subtype, .pred_class)
  mccc <- yardstick::mcc(results, subtype, .pred_class)
  #metrics <- rbind(auc, acc, senss, specc, prec, re, f, kapp, mccc)
  metrics <- data.frame("auc" = auc$.estimate, "senss" = senss$.estimate, "specc" = specc$.estimate,
                        "acc" = acc$.estimate, "prec" = prec$.estimate, "re" = re$.estimate, "f" = f$.estimate,
                        "kapp" = kapp$.estimate, "mcc" = mccc$.estimate)
  return(metrics)
}








