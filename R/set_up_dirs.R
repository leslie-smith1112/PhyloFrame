#' Set up the output directories
#'
#' @param The name of the directory you would like your results to be written in. This will be in the results directory.
#'
#' @return None
#' @export
#'
#' @examples set_up_output_dirs(path_to_desired_result_directory)
#'
set_up_output_dirs <- function(path){
  if(dir.exists(here::here("results")) == FALSE){
    dir.create(here::here("results"))
  }
  if(dir.exists(here::here("results",path)) == FALSE){
    dir.create(here::here("results",path))
    #dir.create(here::here("results",path,"samples"))
    dir.create(here::here("results",path,"model_runs"))
    dir.create(here::here("results",path, "model_runs", "benchmark"))
    dir.create(here::here("results",path,"model_runs","benchmark","afr"))
    dir.create(here::here("results",path,"model_runs","benchmark","admix"))
    dir.create(here::here("results",path,"model_runs","benchmark","eas"))
    dir.create(here::here("results",path,"model_runs","benchmark","eur"))
    dir.create(here::here("results",path,"model_runs","benchmark","mixed"))
    dir.create(here::here("results",path,"model_runs","phyloFrame"))
    dir.create(here::here("results",path,"model_runs","phyloFrame","afr"))
    dir.create(here::here("results",path,"model_runs","phyloFrame","admix"))
    dir.create(here::here("results",path,"model_runs","phyloFrame","eas"))
    dir.create(here::here("results",path,"model_runs","phyloFrame","eur"))
    dir.create(here::here("results",path,"model_runs","phyloFrame","mixed"))
    dir.create(here::here("results",path,"model_runs","phyloFrame","base"))
    dir.create(here::here("results",path,"model_runs","phyloFrame","base","afr"))
    dir.create(here::here("results",path,"model_runs","phyloFrame","base","admix"))
    dir.create(here::here("results",path,"model_runs","phyloFrame","base","eas"))
    dir.create(here::here("results",path,"model_runs","phyloFrame","base","eur"))
    dir.create(here::here("results",path,"model_runs","phyloFrame","base","mixed"))
#    dir.create(here::here("results",path,"ancestry_genes"))
  }
}

#' Set up model specific results direcories.
#'
#' @param primary_dir The directory to write results to.
#' @param curr_ancestry The current ancestry being run.
#' @param i The current model number we are working on
#'
#' @return None
#' @export
#'
#' @examples set_up_model_dir(out_dir, ancestry, model_num)
set_up_model_dir <- function(primary_dir, curr_ancestry,i){
  dir.create(here::here(primary_dir,"benchmark",curr_ancestry,paste0("model_",i)))
  #benchmark_dir <- here::here(primary_dir,"benchmark",curr_ancestry,paste0("model_",i))
  dir.create(here::here(primary_dir,"phyloFrame",curr_ancestry,paste0("model_",i)))
  #phyloFrame_dir <- here::here(primary_dir,"phyloFrame",curr_ancestry,paste0("model_",i))
  dir.create(here::here(primary_dir,"phyloFrame","base",curr_ancestry,paste0("model_",i)))
  #base_genes_dir <- here::here(primary_dir,"phyloFrame","base",curr_ancestry,paste0("model_",i))
}



set_up_output_dirs_single_run <- function(path){
  if(dir.exists(here::here("results")) == FALSE){
    dir.create(here::here("results"))
  }
  if(dir.exists(here::here("results",path)) == FALSE){
    dir.create(here::here("results",path))
    #dir.create(here::here("results",path, "samples"))
    #dir.create(here::here("results",path,"model_runs"))
    # dir.create(here::here("results",path,"model_runs", "benchmark"))
    # dir.create(here::here("results",path,"model_runs", "phyloFrame"))
    # dir.create(here::here("results",path,"model_runs", "phyloFrame","base"))

    #removing the model_runs directory
    dir.create(here::here("results",path, "benchmark"))
    dir.create(here::here("results",path, "phyloFrame"))
    dir.create(here::here("results",path, "phyloFrame","base"))
  }
}





