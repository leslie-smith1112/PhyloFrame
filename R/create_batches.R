
#' Split by Subtype
#'
#' @param disease The disease the samples are from. The sample batches will be written to a "disease_samples" directory.
#' @param ancestry The current ancestry we are creating batches for.
#' @param anc_samples The list of samples for the current ancestry we are creating batches for.
#' @param smallest_batch The smallest number of samples we have for ancestry, this will be used to calculate how many batches we need for an ancestry.
#' @param subtype1_all Data.frame with all samples and their subtype,we will use this to grab samples of the current ancestry
#' @param subtyp2_all Data.frame with all samples and their subtype, we will use this to grab samples of the current ancestry
#'
#' @return None.
#' @export
#'
#' @examples split_subtypes("breast","eur",eur_samples,40,"Basal","Luminal)
split_subtypes <- function(disease, ancestry, anc_samples, smallest_batch, subtype1_all, subtype2_all){
  ## - Helper function should be called once for each ancestry we want to split into batches - ##

  ## - calculate the number of samples for this ancestry - ##
  count <- length(anc_samples)
  ## - calculate the number of batches that will be needed - ##
  num_batches <- floor(count/smallest_batch)
  ## - get total number of samples in each subtype for the ancstry - ##
  sub1_sample_list <- subtype1_all$sample_id[subtype1_all$sample_id %in% anc_samples]
  sub2_sample_list <- subtype2_all$sample_id[subtype2_all$sample_id %in% anc_samples]
  sub1_count <- length(sub1_sample_list)
  sub2_count <- length(sub2_sample_list)

  ## - calculate how many sampels will be in each batch from subtypes - ##
  sub1_batch <- floor(sub1_count/num_batches)
  sub2_batch <- floor(sub2_count/num_batches)
  ## - calculate how many samples will be left over - these samples will be added individaull to batches until theuy are gone
  sub1_leftover <- sub1_count - (sub1_batch * num_batches)
  sub2_leftover <- sub2_count - (sub2_batch * num_batches)

  ## - for as many batches as we have, randomly select samples with the correct numebr of samples without replacement
## COMMENT HERE
  total_samples = 0
  for( i in 1:num_batches){
    if(sub1_leftover != 0){
      batch_size <- sub1_batch + 1
      sub1_leftover <- sub1_leftover - 1
    }else{
      batch_size <- sub1_batch
    }
    add_batch <- sample(x = sub1_sample_list, size = batch_size, replace = FALSE)
    ## - remove samples from the samples list so they are not added to any other sampels - ##
    sub1_sample_list <- sub1_sample_list[!(sub1_sample_list %in% add_batch)]

    if(sub2_leftover != 0){
      batch_size <- sub2_batch + 1
      sub2_leftover <- sub2_leftover - 1
    }else{
      batch_size <- sub2_batch
    }
    add_batch_2 <- sample(x = sub2_sample_list, size = batch_size, replace = FALSE)
    ## - remove samples from the samples list so they are not added to any other sampels - ##
    sub2_sample_list <- sub2_sample_list[!(sub2_sample_list %in% add_batch_2)]

    curr_batch <- c(add_batch, add_batch_2)
    total_samples = total_samples + length(curr_batch)
    curr_batch <- as.data.frame(curr_batch)

    out_dir <- here::here("data-raw",paste0(disease,"_samples"),ancestry)
    if(dir.exists(out_dir) == FALSE){
      dir.create(out_dir)
    }
    write.table(curr_batch,paste0(out_dir,"/samples",i,".tsv"),sep = "\t", row.names = FALSE, col.names = TRUE)
  }

  batch_info <- list("ancestry" = ancestry, "total_samples_assigned" = total_samples,
                     "num_batches" = num_batches, "num_sub1_batch" = sub1_batch, "num_sub2_batch" = sub2_batch,
                     "sub1_leftover" = sub1_leftover, "sub2_leftover" = sub2_leftover)

  return(batch_info)
}

#' Create Batches
#'
#' @param outdir the directory where you want the sample batches for each ancestry to be writen
#' @param sample_ancestry data.frame with the samples for the current diease and their ancestry for only the samples we use in the study
#' @param expression The expression matrix for the idease
#' @param subtype1 The name of subtype 1
#' @param subtype2 The name of subtype 2
#'
#' @return None
#' @export
#'
#' @examples create_batches(samples_dir, cancer_type, expression, subtype1, subtype2)
create_batches <- function(disease, samples_ancestry, expression, subtype1, subtype2 ){
  ## - this function is run for a single disease at a time - ##
  ## - samples ancestry is already cut down to samples we are using in our study for this disease - ##
  ## - split samples into associated acnestries - ##
  eur <- samples_ancestry[samples_ancestry$consensus_ancestry == "eur",]$patient
  afr <- samples_ancestry[samples_ancestry$consensus_ancestry == "afr",]$patient
  eas <- samples_ancestry[samples_ancestry$consensus_ancestry == "eas",]$patient
  admix <- samples_ancestry[samples_ancestry$consensus_ancestry %in% c("eas_admix","afr_admix","eur_admix","admix"),]$patient
  ## - get count for each ancestry, we need the count for finding the smallest ancestry size to split into correct number of batches later on - ##
  eur_count <- length(eur)
  afr_count <- length(afr)
  eas_count <- length(eas)
  admix_count <- length(admix)

  ## - for downsampling europeans in mixed batch - ##
  ## - find that max number of samples to determine how many european sampels to add to mixed batches - ##
  max_anc_samples <- max(afr_count, eas_count, admix_count)
  ## - random sampling wihtout replacement of eur samples for mixed batches - ##
  eur_selected <- sample(eur, max_anc_samples, replace = FALSE)


  ## - we keep all samples and all_count to use to make mixed batches later - ##
  all_samples <- c(eur_selected, afr, eas, admix)
  all_count <- length(eur_selected) + afr_count + eas_count + admix_count

  ## - for this disease, find the ancestry with the smallest number of samples so we know how many samples to add to our batches - ##
  smallest_batch <- min(eur_count, afr_count, eas_count, admix_count)

  cut_dat <- data.frame(rownames(expression), expression$subtype)
  colnames(cut_dat) <- c("sample_id","subtype")

  ## - seperate out subtypes so that each batch has ~even number of each - ##
  ## - currently this has all subtypes for all ancestries - ##
  subtype1_dat <- cut_dat[cut_dat$subtype == subtype1,]
  subtype2_dat <- cut_dat[cut_dat$subtype == subtype2,]

  sample_dir <- here::here("data-raw",paste0(disease,"_samples"))
  if(dir.exists(sample_dir) == FALSE){
    dir.create(sample_dir)
  }

  ## - disease only used for naming directories - ##
  disease <- tolower(disease)
  #return values added for testing
  eur_stats <- split_subtypes(disease,"eur",eur,smallest_batch,subtype1_dat,subtype2_dat)
  afr_stats <- split_subtypes(disease,"afr",afr,smallest_batch,subtype1_dat,subtype2_dat)
  eas_stats <- split_subtypes(disease,"eas",eas, smallest_batch,subtype1_dat,subtype2_dat)
  admix_stats <- split_subtypes(disease, "admixed",admix,smallest_batch,subtype1_dat,subtype2_dat)
  mixed_stats <- split_subtypes(disease, "mixed",all_samples, smallest_batch,subtype1_dat,subtype2_dat)

  value_check <- list("EUR" = eur_stats, "AFR" = afr_stats, "EAS" = eas_stats, "ADMIX" = admix_stats,
                      "MIXED" = mixed_stats)
  # message("EUR total samples: ", value_check$EUR$total_samples,", # of batches: ",
  #         value_check$EUR$num_batches,", # of subtype1 in each batch: ",
  #         value_check$EUR$sub1_batch, ", # of subtype2 in each batch: ",
  #         value_check$EUR$sub2_batch, ", samples leftover from each subtype: ",
  #         sub1_leftover," ",sub2_leftover)
  return(value_check)
}

