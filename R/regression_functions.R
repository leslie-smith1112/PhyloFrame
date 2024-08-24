#Functions for running the elastic net model and testing on other ancestries

#' Elastic Net Run
#'
#' @param train_expression The expression matrix you would like to train on. Expects samples as rows and genes as columns.
#' @param out_dir The name of the directory where you would like your output to be written to. This will be written inside the results directory. For example if I pass in "breast_cancer" my results will be in ~/results/breast_cancer/.
#' @param model_num The number of the model currently running for the ancestry. This is primarily for file naming.
#' @param lr_mixture The mixture we want for the model. In phyloFrame this is initially set to 1. In the second pass it is set to 0.
#' @param seed Randomization seed. Default is 4831.
#'
#' @return The trained model.
#' @export
#'
#' @examples lasso_model <- run_elastic_net(training_set, out_dir, i, 1) (Here i is defined in the for loop) Seed uses default.
run_elastic_net <- function(train_expression, out_dir, model_num, lr_mixture, seed = 4831){
  #Resource: https://www.tidymodels.org/start/case-study/
  set.seed(seed)
  ## - split samples into training and test sets - ##
  #OLD:
  split <- rsample::initial_split(train_expression, strata = subtype)
  training_set <- rsample::training(split)
  test_set <- rsample::testing(split)
  validation_set <- rsample::validation_split(training_set,
                                                strata = subtype,
                                                prop = 0.80)
  write.table(training_set, paste0(out_dir,"/model_",model_num,"_training_set.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
  write.table(test_set, paste0(out_dir,"/model_",model_num,"_test_set.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)


  #NEW:
  # split <- rsample::initial_validation_split(train_expression, prop =  c(0.6, 0.2),strata = subtype)
  # training_set <- rsample::training(split)
  # test_set <- rsample::testing(split)
  # validation_set <- rsample::validation_set(split)
  #val <- rsample::vfold_cv(validation_set, v = nrow(validation_set), repeats = 1,strata = subtype)

  ## - make sure there is no intersection between the training and test set - ##
  if(length(intersect(rownames(training_set), rownames(test_set))) != 0){
    message("ERROR: There is overlap between the training and test sets.")
    stop()
  }

  ## - define the model workflow, we will tune the penalty - ##
  lr_model <-
    parsnip::logistic_reg(penalty = tune(), mixture = lr_mixture) %>%
    parsnip::set_engine("glmnet")
  lr_recipe <-
    recipes::recipe(subtype ~ ., data = training_set)
  lr_workflow <- workflows::workflow() %>%
    workflows::add_model(lr_model) %>%
    workflows::add_recipe(lr_recipe)

  ## - create our grid for tuning - ##
  lr_tune_grid <- tibble::tibble(penalty = 10^seq(-4,-1,length.out = 30))

  ## - train and save the results of the model on the validation set - ##
  init_model_fit <-
    lr_workflow %>%
    tune::tune_grid(validation_set,
                    grid = lr_tune_grid,
                    control = tune::control_grid(save_pred = TRUE),
                    metrics = yardstick::metric_set(yardstick::roc_auc))



  penalties <- init_model_fit %>% tune::collect_metrics()
  top_score <- penalties[penalties$mean == max(penalties$mean),]
  final_penalty <- top_score[sample(1:nrow(top_score),1),]$penalty

  lr_final_mod <-
    parsnip::logistic_reg(penalty = final_penalty, mixture = lr_mixture) %>%
    parsnip::set_engine("glmnet", importance = "impurity")

  lr_final_workflow <-
    lr_workflow %>%
    workflows::update_model(lr_final_mod)

  final_fit <- tune::last_fit(lr_final_workflow,split)
  final_prediction <- tune::collect_metrics(final_fit)
  write.table(final_prediction,paste0(out_dir,"/model_",model_num,"_metrics.tsv"),sep = "\t", col.names = TRUE, row.names = FALSE)

  ## - write genes used for prediction to file - ##
  gene_scores <- final_fit %>%
    workflows::extract_fit_engine() %>%
    vip::vi()

  gene_scores <- gene_scores[gene_scores$Importance > 0,]
  write.table(gene_scores, paste0(out_dir,"/model_",model_num,"_all_sig.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)

  my_model <- tune::extract_workflow(final_fit)
  readr::write_rds(my_model, paste0(out_dir,"/model_",model_num,"_EN_model.rds"))
  ## - save model coefficients - ##
  save_model <- broom::tidy(my_model)
  write.table(save_model, paste0(out_dir,"/model_",model_num,"_model_coefficients.tsv"),sep = "\t", col.names = TRUE, row.names = FALSE)

  return(my_model)

}

#' Model Metrics on Testing Data
#'
#' @param my_model The trained model.
#' @param test_expression The expression matrix we are testing on.
#' @param anc The current ancestry whose expression matrix we are testing on, only used for file naming.
#' @param out_dir The directory name we will be writing to. (This is the directory within results directory that we defined initially).
#' @param model_num The number of the model we are testing - only used for naming purposes
#' @param subtype_1 One of the binary subtypes we predicted.
#' @param subtype_2 The other subtype we predicted.
#'
#' @return None
#' @export
#'
#' @examples test_model(benchmark_model, afr_expression, "afr", benchmark_dir, j, subtype1, subtype2) (Here j is defined in the for loop)
test_model <- function(my_model, test_expression, anc, out_dir, model_num, subtype_1, subtype_2){
  file_prefix <- paste0(anc,"_", model_num)
  pred_prob <- predict(my_model, test_expression, type = "prob")
  pred_class <- predict(my_model, test_expression, type = "class")
  results <- test_expression %>% dplyr::select(subtype) %>% bind_cols(pred_class, pred_prob)
  results <- tibble::rownames_to_column(results, "sample_id")
  write.table(results, paste0(out_dir,"/",file_prefix,"_results.tsv"),sep = "\t", col.names = TRUE, row.names = FALSE)

  result_metrics <-  bind_cols(test_expression$subtype, results[,4], results[,5], results$.pred_class)
  colnames(result_metrics) <- c("subtype", paste0(".pred_",subtype_1), paste0(".pred_",subtype_2), ".pred_class")
  auc <- yardstick::roc_auc(result_metrics, subtype, paste0(".pred_",subtype_1))

  confusion <- yardstick::conf_mat(results, truth = subtype,estimate = .pred_class)
  conf_matrix <- as.data.frame(confusion$table)
  write.table(conf_matrix, paste0(out_dir,"/",file_prefix,"_confusion_matrix.tsv"), sep  = "\t")

  senss <- yardstick::sens(results,  subtype, .pred_class)
  specc <- yardstick::spec(results,  subtype,  .pred_class)
  acc <- yardstick::accuracy(results,  subtype, .pred_class)
  prec <- yardstick::precision(results,  subtype, .pred_class)
  re <- yardstick::recall(results, subtype, .pred_class)
  f <- yardstick::f_meas(results, subtype, .pred_class)
  kapp <- yardstick::kap(results,  subtype, .pred_class)
  mccc <- yardstick::mcc(results, subtype, .pred_class)
  metrics <- rbind(auc, acc, senss, specc, prec, re, f, kapp, mccc)
  write.table(metrics,paste0(out_dir,"/",file_prefix,"_metrics.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
}

#my_model, test_expression, anc, out_dir, model_num, subtype_1, subtype_2
## MODEL METRICS FOR SUBSAHARAN AFRICAN VALIDATION SET:
test_model_sub_sah_dataset <- function(my_model, test_expression, directory, out_file, subtype1, subtype2)
{
  pred_prob <- predict(my_model, test_expression, type = "prob")
  pred_class <- predict(my_model, test_expression, type = "class")
  results <- test_expression %>% dplyr::select(subtype) %>% bind_cols(pred_class, pred_prob)
  results <- tibble::rownames_to_column(results, "sample_id")
  write.table(results, paste0(directory,"/", out_file,"_results.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
  false_negative <- nrow(results[results$.pred_class == "Luminal",])
  true_positive <- nrow(results[results$.pred_class == "Basal",])
  ## - compute by hand because we have no luminal sampels - ##
  rec <- true_positive/(true_positive + false_negative)
  prec <- true_positive/(true_positive + 0)
  #spec cannot calculate
  # cannot calculate: mcc = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (0 + 0) * (0 + FN))

  to.write <- data.frame("Recall"= rec, "Precision" = prec)
  write.table(to.write, paste0(directory,"/", out_file,"_metrics.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
}
