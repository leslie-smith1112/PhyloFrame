## ADXMIXTURE PLOT ##

## - Results directory naems:
TCGA_Breast_Gnomad4_corrected
ancestry
#eur model was tested on eur, eas, and afr
models <- 1:17
ancestries <- c("eur","eas","afr","admix")
combinations <- expand.grid(Model = models, Ancestry = ancestries)
## - example file name in the eur ancestry model results: eur_14_results.tsv
result_pf <- here::here("results","TCGA_Breast_Gnomad4_corrected","model_runs","phyloFrame","eur",paste0("model_",combinations$Model),paste0(combinations$Ancestry,"_",combinations$Model,"_results.tsv"))
result_pf
## function to read in all metrics - defined below
pf_dat <- get_admixture_performance(result_pf, ancestry, "PF_Eur")
pf_dat$Model <- "PhyloFrame"
pf_dat$percent_correct <- as.numeric(pf_dat$percent_correct)

result_bm <- here::here("results","TCGA_Breast_Gnomad4_corrected","model_runs","benchmark","eur",paste0("model_",combinations$Model),paste0(combinations$Ancestry,"_",combinations$Model,"_results.tsv"))
result_bm
bm_dat <- get_admixture_performance(result_bm, ancestry, "BM_Eur")
bm_dat$Model <- "Benchmark"
bm_dat$percent_correct <- as.numeric(bm_dat$percent_correct)

## - seperate primary ancestry we are looking at and take only those with that ancestry as primary admixture
afr_dat <- rbind(pf_dat[(pf_dat$consensus_ancestry == "afr" & pf_dat$admixture >= 0.50),],bm_dat[(bm_dat$consensus_ancestry == "afr" & bm_dat$admixture >= 0.50),])
eur_dat <- rbind(pf_dat[(pf_dat$consensus_ancestry == "eur" & pf_dat$admixture >= 0.50),],bm_dat[(bm_dat$consensus_ancestry == "eur" & bm_dat$admixture >= 0.50),])
eas_dat <- rbind(pf_dat[(pf_dat$consensus_ancestry == "eas" & pf_dat$admixture >= 0.50),],bm_dat[(bm_dat$consensus_ancestry == "eas" & bm_dat$admixture >= 0.50),])

####### PLOTS #######
ggplot(data= eur_dat, aes(x=admixture, y=percent_correct, group=Model, color = Model, fill = Model))+  geom_smooth(method = loess, level = 0.95) +geom_point() + scale_y_continuous(limits=c(0, 100)) +
  scale_x_continuous(limits=c(0.5, 1.00)) +
  geom_hline(yintercept=50,linetype=2) + theme_minimal() +
  scale_color_manual(values = c(
    "PhyloFrame" = "#CA4136",
    "Benchmark" = "#f5a2a5"
  )) +
  scale_fill_manual(values = c("PhyloFrame" = "#CA4136",
                               "Benchmark" = "#f5a2a5"))
ggsave(here::here("figures","BRCA_Admixture_Eur.png"),plot = last_plot(),width = 7, height = 5)
nrow(eur_dat)/2

ggplot(data= afr_dat, aes(x=admixture, y=percent_correct, group=Model, color = Model, fill = Model))+  geom_smooth(method = loess, level = 0.95) +geom_point() + scale_y_continuous(limits=c(0, 100)) +
  scale_x_continuous(limits=c(0.5, 1.00)) +
  geom_hline(yintercept=50,linetype=2) + theme_minimal() +
  scale_color_manual(values = c("PhyloFrame" = "#1F619E",
                                "Benchmark" = "#9fccfa")) +
  scale_fill_manual(values = c("PhyloFrame" = "#1F619E",
                               "Benchmark" = "#9fccfa"))

ggsave(here::here("figures","BRCA_Admixture_Afr.png"),plot = last_plot(),width = 7, height = 5)
nrow(afr_dat)/2

ggplot(data= eas_dat, aes(x=admixture, y=percent_correct, group=Model, color = Model, fill = Model))+  geom_smooth(method = loess, level = 0.95) +geom_point() + scale_y_continuous(limits=c(0, 100)) +
  scale_x_continuous(limits=c(0.5, 1.00)) +
  geom_hline(yintercept=50,linetype=2) + theme_minimal() +
  scale_color_manual(values = c(
    "PhyloFrame" = "#496849",
    "Benchmark" = "#63b34c"
  )) +
  scale_fill_manual(values = c("PhyloFrame" = "#496849",
                               "Benchmark" = "#63b34c"))

ggsave(here::here("figures","BRCA_Admixture_Eas.png"),plot = last_plot(),width = 7, height = 5)
dim(eas_dat)
nrow(eas_dat)/2
#' Admixture Plot
#'
#' @param file_list List of the model results to plot
#' @param ancestry_info Ancestry information for samples
#' @return datafrae with ancestry and performance for every sample for every admixture
#' @export
#'
#' @examples
get_admixture_performance <- function(file_list, ancestry_info, BMorPF){
  pf_eur_model <- lapply(file_list, readr::read_tsv)
  pf_dat <- do.call(rbind, pf_eur_model)
  pf_dat[1:5,1:5]

  ## add in the ancestry info - merging by sample names
  pf_all <- merge(pf_dat, ancestry_info, by.x = "sample_id", by.y = "patient")
  dim(pf_all)
  head(pf_all)
  pf_all$model <- BMorPF
  pf_all$percent_correct <- NA
  ## - dataframe to keep only one instance of each sample, their ancestry, and how many times they were predicted correctly
  pf_new <- data.frame(matrix(ncol = 9, nrow = 0))

  colnames(pf_new) <- c("sample_id","consensus_ancestry","Admixture % AFR","Admixture % AMR",
                        "Admixture % EAS", "Admixture % EUR", "Admixture % SAS", "Model", "percent_correct")
  ## - lose "admix" ending if there is one
  pf_all$consensus_ancestry <- substr(pf_all$consensus_ancestry, 1, 3)

  ## - for every unique sample, find the percent of time it was correctly predicted
  samples <- unique(pf_all$sample_id)
  for(i in 1:length(samples)){
    temp.dat <- pf_all[pf_all$sample_id == samples[i],] #get all predictions of sample
    total <- nrow(temp.dat)
    correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,]) # find the number of times it was correcly predicted
    percent <- (correct/total) * 100
    pf_all[pf_all$sample_id == samples[i],]$percent_correct <- percent # add prediction percent to sample info
    new_list <- c(temp.dat$sample_id[1], temp.dat$consensus_ancestry[1], temp.dat$`Admixture % AFR`[1], # create new row with all wanted infomation
                  temp.dat$`Admixture % AMR`[1], temp.dat$`Admixture % EAS`[1], temp.dat$`Admixture % EUR`[1],
                  temp.dat$`Admixture % SAS`[1], temp.dat$model[1],percent)
    pf_new[nrow(pf_new) + 1,]<-new_list
  }
  #convert the admixture % to numbers
  j <- c("Admixture % EAS", "Admixture % AFR", "Admixture % EUR", "Admixture % SAS", "Admixture % AMR")
  pf_new[, j] <- apply(pf_new[, j], 2, function(x) as.numeric(x))

  ## - create dataframe for plotting so that each ancestry % has its own row for each sample - then we can take specific admixture slices

  dat_names <- c("sample_id", "consensus_ancestry", "admixture", "model",
                 "percent_correct", "ancestry_percent")
  afr <- pf_new %>% dplyr::select(-`Admixture % AMR`, -`Admixture % EAS`, -`Admixture % EUR`,-`Admixture % SAS`)
  afr$ancestry_slice <- "afr"
  colnames(afr) <- dat_names
  dim(afr)
  eas <- pf_new %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`,-`Admixture % EUR`,-`Admixture % SAS`)
  eas$ancestry_slice <- "eas"
  colnames(eas) <- dat_names
  eur <- pf_new %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`, -`Admixture % EAS`,-`Admixture % SAS`)
  eur$ancestry_slice <- "eur"
  colnames(eur) <- dat_names
  sas <- pf_new %>% dplyr::select(-`Admixture % AMR`, -`Admixture % AFR`, -`Admixture % EAS`,-`Admixture % EUR`)
  sas$ancestry_slice <- "sas"
  colnames(sas) <- dat_names
  amr <- pf_new %>% dplyr::select(-`Admixture % SAS`, -`Admixture % AFR`,-`Admixture % EAS`,-`Admixture % EUR`)
  amr$ancestry_slice <- "amr"
  colnames(amr) <- dat_names

  all_ancestry <- rbind(afr, eas, eur, sas, amr)
  all_ancestry <- all_ancestry %>% tidyr::drop_na(admixture)
  return(all_ancestry)
}






