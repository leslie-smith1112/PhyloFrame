


# Pull the data for each cancer, combine into 1 matrix
dat.brca <- read.table(here::here("results","TCGA_Breast_Gnomad4_corrected", "model_runs","batch_predictions","master_predictions.tsv"), sep='\t', header=T, stringsAsFactors=F)
dat.ucec <- read.table(here::here("results","TCGA_Uterine_Gnomad4_corrected", "model_runs","batch_predictions","master_predictions.tsv"), sep='\t', header=T, stringsAsFactors=F)
dat.thyr <- read.table(here::here("results","TCGA_Thyroid_Gnomad4_corrected", "model_runs","batch_predictions","master_predictions.tsv"), sep='\t', header=T, stringsAsFactors=F)


dat.brca$Cancer <- 'BRCA'
dat.ucec$Cancer <- 'UCEC'
dat.thyr$Cancer <- 'THYR'

## here train-data is the model training data, model_num is the model for the trainined data, ancestry_batch is the test data and batch_num is the batch for the test data
dat.brca <- dat.brca[(dat.brca$train_data != "mixed" & dat.brca$ancestry_batch != "mixed"),]
dat.ucec <- dat.ucec[(dat.ucec$train_data != "mixed" & dat.ucec$ancestry_batch != "mixed"),]
dat.thyr <- dat.thyr[(dat.thyr$train_data != "mixed" & dat.thyr$ancestry_batch != "mixed"),]
#master dat with all diseases


#color by ancestry batch

#for each ancestry in ancestries, subset the master dat for this disease to get the test samples of that ancestry trained model and plot them

#brca_gathered <- pivot_longer(brca_mcc, c('pf_mcc','bm_mcc'), names_to='Model', values_to='Mcc')

plot_ancestry <- function(ancestry,dat){
  ancestry_performance <- dat[dat$train_data == ancestry,]
  cols <- c(eur="#CA4136",admix="#FDC652",eas="#496849",afr="#1F619E") # Keep colors same across cancer plots
  dat <- dat[dat$ancestry != "mixed",]
  ## mcc
  p <- ancestry_performance %>%
    ggplot(aes(x=bm_mcc, y=pf_mcc, color=ancestry_batch)) +
    #scale_colour_manual(values=cols) +
    #scale_fill_manual(values = cols) +
    geom_point(alpha = 0.6, size = 5) +
    geom_smooth(color = "grey20", fill = "grey40",method = lm, linewidth = 2) +
    geom_abline(linetype = 2, intercept=0)+
    lims( x=c(-1,1), y=c(-1,1) )+ scale_colour_manual(values=cols)+
    labs(title = "95% Confidence Interval") +
    theme_bw() +
    #theme(plot.title = element_text(face = "bold",hjust = 0.5))
    ggplot2::theme(axis.text=element_text(size=20),
                   axis.title=element_text(size=20,face="bold")) + ggplot2::ggtitle(paste0(toupper(ancestry)," MODEL MCC FOR ",toupper(disease))) + ggplot2::geom_abline()


  p


  ggsave(here::here("figures","Model_Scatter_Plots",paste0(ancestry_performance[1,"Cancer"],ancestry,"_mcc_CI_LM.png")),plot = last_plot(),width = 9, height = 7)
  ggsave(here::here("figures","Model_Scatter_Plots",paste0(ancestry_performance[1,"Cancer"],ancestry,"_mcc_CI_LM.svg")),plot = last_plot(),width = 9, height = 7)

  #ROC
  # ggplot(ancestry_performance, aes(color=ancestry_batch,x=bm_auc, y=pf_auc)) + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols) +
  #   theme_bw() +
  #   theme(axis.text=element_text(size=20),
  #         axis.title=element_text(size=20,face="bold")) +geom_point(alpha = 0.6, size = 5) + ggtitle(paste0(toupper(ancestry),"_",ancestry_performance[1,"Cancer"], "_AUC"))
  #
  # ggsave(here::here("figures","Model_Scatter_Plots",paste0(ancestry_performance[1,"Cancer"],ancestry,"_auc.png")),plot = last_plot(),width = 9, height = 7)
  # ggsave(here::here("figures","Model_Scatter_Plots",paste0(ancestry_performance[1,"Cancer"],ancestry,"_auc_.svg")),plot = last_plot(),width = 9, height = 7)
}

## define ancestry_performance[1,"Cancer"]
disease_dat <- dat.brca #dat.brca, dat.ucec, dat.thyr

# get list of trained ancestries for this disease model
#ancestries <- unique(disease_dat$train_data)
ancestries <- c("admix","afr")
lapply(ancestries, function(ancestry){
  plot_ancestry(ancestry,disease_dat)
})

