############## SUPP FIG 5 ######
require(readr)  # for read_csv()
require(dplyr) #data frame handling
library(tidyr)
library(stringi)
library(ggplot2)
####################################  DEFINE DISEASE AND MODEL NUMBERS FOR EACH ANCESTRY ####################################

## THYROID CALLS
disease <- "thyroid"
directory <- "TCGA_Thyroid_Gnomad4_corrected"
sample_predictions(directory, "admix", disease, FALSE) # FALSE = only has one model
sample_predictions(directory, "afr", disease, FALSE)
sample_predictions(directory, "eas", disease, TRUE)
sample_predictions(directory, "eur", disease, TRUE)

## BREAST CALLS
disease <- "Breast"
directory <- "TCGA_Breast_Gnomad4_corrected"
sample_predictions(directory, "admix", disease, FALSE)
sample_predictions(directory, "afr", disease, TRUE)
sample_predictions(directory, "eas", disease, FALSE)
sample_predictions(directory, "eur", disease, TRUE)


## UTERINE CALLS
disease <- "Uterine"
directory <- "TCGA_Uterine_Gnomad4_corrected"
sample_predictions(directory, "admix", disease, FALSE)
sample_predictions(directory, "afr", disease, TRUE)
sample_predictions(directory, "eur", disease, TRUE)

## for the given disease and the given directory
# read in the final predictions file that has all model predictions
# merge file with ancestry
#for every ancestry model (substset of the merged filee)
# for every sample get hte predcitions and their ancestry

sample_predictions <- function(directory, models, disease, normal){
  predictions<- readr::read_tsv(here::here("results",directory,"model_runs","final_results.tsv"))
  colnames(predictions) <- c("sample_id", "subtype",".pred_class", ".pred_sub1",".pred_sub2","train_data","model_num","test_data","version")
  predictions <- predictions[predictions$sample_id != "sample_id",]

  anc <- ancestry |> dplyr::select(patient, consensus_ancestry)
  anc$ancestry <- plyr::mapvalues(anc$consensus_ancestry, c("admix", "afr", "afr_admix", "eas", "eas_admix", "eur", "eur_admix"), c("admix", "afr", "admix", "eas", "admix","eur","admix"))
  all_pred <- merge(predictions, anc, by.x = "sample_id", by.y = "patient")

  sample_list <- unique(all_pred$sample_id)
  calculate_correctness(sample_list, models,all_pred, disease, normal)
  #return(dat)
  # lapply(models,function(model){
  #   calculate_correctness(sample_list, models,all_pred, disease)
  # })

}


calculate_correctness <- function(samples, model, master_predictions, disease, normal){
  ## PHYLOFRAME
  anc_predictions <- master_predictions[(master_predictions$train_data == model) & (master_predictions$version == "phyloFrame"),]
  new_pf <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(new_pf ) <- c("Sample_id", "Ancestry", "Model", "Percent_correct", "Version")
  for(i in 1:length(samples)){
    ## - get a single samples predictions
    temp.dat <- anc_predictions[(anc_predictions$sample_id == sample_list[i]),] #get all predictions of sample
    total <- nrow(temp.dat)
    correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
    percent <- (correct/total) * 100
    new.list <- c(temp.dat$sample_id[1], paste0(temp.dat$ancestry[1],"_PhyloFrame"),
                  temp.dat$train_data[1],percent, "PhyloFrame" )
    new_pf[nrow(new_pf) + 1,] <-new.list
  }

  ## BENCHMARK
  anc_predictions <- master_predictions[(master_predictions$train_data == model) & (master_predictions$version == "benchmark"),]
  new_bm <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(new_bm ) <- c("Sample_id", "Ancestry", "Model", "Percent_correct", "Version")
  for(i in 1:length(samples)){
    ## - get a single samples predictions
    temp.dat <- anc_predictions[(anc_predictions$sample_id == sample_list[i]),] #get all predictions of sample
    total <- nrow(temp.dat)
    correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
    percent <- (correct/total) * 100
    new.list <- c(temp.dat$sample_id[1], paste0(temp.dat$ancestry[1],"_Benchmark"),
                  temp.dat$train_data[1],percent, "Benchmark" )
    new_bm[nrow(new_bm) + 1,] <-new.list
  }
  all_dat <- rbind(new_bm, new_pf)
  all_dat$Percent_correct <- as.numeric(all_dat$Percent_correct)

  if(normal == TRUE){
    ################################# NORMAL MODEL
    all_dat$Ancestry <- factor(all_dat$Ancestry, levels = c("admix_PhyloFrame","admix_Benchmark",
                                                            "afr_PhyloFrame","afr_Benchmark",
                                                            "eas_PhyloFrame","eas_Benchmark",
                                                            "eur_PhyloFrame","eur_Benchmark"))
    p <- ggplot(all_dat, aes(x=Ancestry, y=Percent_correct, fill = Ancestry, color = Ancestry)) + geom_boxplot() +
      scale_fill_manual(values = c("#FDC652","#ffe943","#1F619E", "#9fccfa","#496849","#63b34c","#CA4136","#f5a2a5"
      )) + scale_color_manual(values = c("#ca9e43","#cab836","#153e64", "#7b9dbf","#2b3d2b","#427533","#7c2a23","#b97d7f"))+ theme_bw() +
      ggtitle(paste0(toupper(model) ,"_MODEL_", toupper(disease))) +
      ylim(0, 100) +
      theme(axis.text=element_text(size=15),
            axis.title=element_text(size=15,face="bold")) #; ggsave('THYR_MCC.png',width=5.1,height=4)
    ggsave(here::here("figures","TCGA_Sample_Sample_Predictions",paste0(disease,"_",model,".png")),plot = last_plot(),width = 7, height = 5)

    ggsave(here::here("figures","TCGA_Sample_Sample_Predictions",paste0(disease,"_",model,".svg")),plot = last_plot(),width = 7, height = 5)

    p
  }else if(normal == FALSE & model == "admix"){
    ################################## ADMIXED SINGLE MODEL

    all_dat <- all_dat[!(is.na(all_dat$Sample_id)),]
    all_dat$Ancestry <- factor(all_dat$Ancestry, levels = c("afr_PhyloFrame","afr_Benchmark",
                                                            "eas_PhyloFrame","eas_Benchmark",
                                                            "eur_PhyloFrame","eur_Benchmark"))

    p <- ggplot(all_dat, aes(x=Ancestry, y=Percent_correct, fill = Ancestry, color = Ancestry)) + geom_boxplot() +
      scale_fill_manual(values = c("#1F619E", "#9fccfa","#496849","#63b34c","#CA4136","#f5a2a5"
      )) + scale_color_manual(values = c("#153e64", "#7b9dbf","#2b3d2b","#427533","#7c2a23","#b97d7f"))+ theme_bw()+
      ggtitle(paste0(toupper(model) ,"_MODEL_", toupper(disease))) +
      theme(axis.text=element_text(size=10),
            axis.title=element_text(size=10,face="bold")) #; ggsave('THYR_MCC.png',width=5.1,height=4)
    ggsave(here::here("figures","TCGA_Sample_Sample_Predictions",paste0(disease,"_",model,".png")),plot = last_plot(),width = 7, height = 5)
    ggsave(here::here("figures","TCGA_Sample_Sample_Predictions",paste0(disease,"_",model,".svg")),plot = last_plot(),width = 7, height = 5)

    p
  }else if(normal == FALSE & model == "afr"){
    ################################# AFR SINGLE MODEL
    all_dat <- all_dat[!(is.na(all_dat$Sample_id)),]
    all_dat$Ancestry <- factor(all_dat$Ancestry, levels = c("admix_PhyloFrame","admix_Benchmark",
                                                            "eas_PhyloFrame","eas_Benchmark",
                                                            "eur_PhyloFrame","eur_Benchmark"))

    p <- ggplot(all_dat, aes(x=Ancestry, y=Percent_correct, fill = Ancestry, color = Ancestry)) + geom_boxplot() +
      scale_fill_manual(values = c("#FDC652","#ffe943","#496849","#63b34c","#CA4136","#f5a2a5"
      )) + scale_color_manual(values = c("#ca9e43","#cab836","#2b3d2b","#427533","#7c2a23","#b97d7f"))+ theme_bw()+
    ggtitle(paste0(toupper(model) ,"_MODEL_", toupper(disease))) +
      theme(axis.text=element_text(size=15),
            axis.title=element_text(size=15,face="bold")) #; ggsave('THYR_MCC.png',width=5.1,height=4)
    ggsave(here::here("figures","TCGA_Sample_Sample_Predictions",paste0(disease,"_",model,".png")),plot = last_plot(),width = 7, height = 5)
    ggsave(here::here("figures","TCGA_Sample_Sample_Predictions",paste0(disease,"_",model,".svg")),plot = last_plot(),width = 7, height = 5)

    p
  }else if(normal == FALSE & model == "eas"){
    ################################# EAS SINGLE MODEL
    all_dat <- all_dat[!(is.na(all_dat$Sample_id)),]
    all_dat$Ancestry <- factor(all_dat$Ancestry, levels = c("admix_PhyloFrame","admix_Benchmark",
                                                            "afr_PhyloFrame","afr_Benchmark",
                                                            "eur_PhyloFrame","eur_Benchmark"))

    p <- ggplot(all_dat, aes(x=Ancestry, y=Percent_correct, fill = Ancestry, color = Ancestry)) + geom_boxplot() +
      scale_fill_manual(values = c("#FDC652","#ffe943","#1F619E", "#9fccfa","#CA4136","#f5a2a5"
      )) + scale_color_manual(values = c("#ca9e43","#cab836","#153e64", "#7b9dbf","#7c2a23","#b97d7f"))+ theme_bw()+
      ggtitle(paste0(toupper(model) ,"_MODEL_", toupper(disease))) +
      theme(axis.text=element_text(size=15),
            axis.title=element_text(size=15,face="bold")) #; ggsave('THYR_MCC.png',width=5.1,height=4)
    ggsave(here::here("figures","TCGA_Sample_Sample_Predictions",paste0(disease,"_",model,".png")),plot = last_plot(),width = 7, height = 5)
    ggsave(here::here("figures","TCGA_Sample_Sample_Predictions",paste0(disease,"_",model,".svg")),plot = last_plot(),width = 7, height = 5)

    p
  }else{
    message("Issue plotting ", model)
  }
}



##### BRIEF ADD FOR CHECKING UTERINE SEROUS SAMPLES ####
# clinical <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/diseases/uterine/ucec_tcga_pan_can_atlas_2018_clinical_data.tsv")
# clin <- clinical %>% dplyr::select(`Sample ID`,`Tumor Type`)
# table(clin$`Tumor Type`)
# keeping <- to.return %>% dplyr::select(sample_id,percent_correct)
# test <- merge(keeping, clin, by.x = "sample_id", by.y = "Sample ID")
# head(test)
# table(test$percent_correct,test$`Tumor Type`)

# estimated_ancestry <- readxl::read_xlsx("/home/leslie.smith1/blue_kgraim/leslie.smith1/phyloFrame/all_ancestry_runs/tcga_estimated_ancestry.xlsx",skip = 1)
# samples.ancestry <- estimated_ancestry[estimated_ancestry$tumor_type == "UCEC",]
# samples.ancestry <- samples.ancestry %>% dplyr::select(patient, consensus_ancestry)
# samples.ancestry$patient <- paste0(samples.ancestry$patient, "-01")
# ancestry.ucec <- merge(samples.ancestry, clin, by.x = "patient", by.y = "Sample ID")
# NOTE: EUR model each sample predicted N times, 760 Endometrial predicted 100 a few at 70-90, Serous: between 10-45 % correctly predicted
# ADM one model endometrical 50 wrong 674 right, serous 29 wrong 151 right
#AFr endometrial 757 haf 100% right, 6 had 0 % right and 9 had 50%, Serous had 71 had 0% right, 60 had 100% right, 79 had 50% right
#MIXED endometiral had 756 correct and only 16 below that, serous had 49 wrong and all over tha map basically, only 5 got 100%
