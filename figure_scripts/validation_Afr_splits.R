## laod project for validation metadata
devtools::load_all()
## dataframe with all predictions and metadata


dat <- all_vali_preds
# for when we want ot subset to a specific model
#dat <- mega_validation_meta
#dat <- mega_validation_meta[mega_validation_meta$train_data == "eur",]

dat$ancestry_assigned_by_us[ dat$ancestry_assigned_by_us=='None' ] <- 'Unknown'
dat$ancestry_assigned_by_study[ dat$ancestry_assigned_by_study=='Other' ] <- NA

dat$PF_pred_correct <- dat$pf_pred_class==dat$subtype
dat$BM_pred_correct <- dat$bm_pred_class==dat$subtype

dat2 <- dat[,c('sample_id','PF_pred_correct','BM_pred_correct','ancestry_assigned_by_study','ancestry_assigned_by_us','model_num')]
dat3 <- xtabs(PF_pred_correct ~ sample_id + model_num, dat2)
dat4 <- xtabs(BM_pred_correct ~ sample_id + model_num, dat2)
tmp <- dat2[,c('sample_id','ancestry_assigned_by_study','ancestry_assigned_by_us')]
tmp$all <- paste0( tmp$sample_id, tmp$ancestry_assigned_by_study, tmp$ancestry_assigned_by_us )
tmp2 <- tmp[!duplicated(tmp$all),]
tmp2 <- tmp2[tmp2$all != "NANANA",]
num_models  <- length(unique(paste0(dat$train_data,dat$model_num)))
tmp2$PF_pct <- apply(dat3, 1, sum)/num_models * 100
tmp2$BM_pct <- apply(dat4, 1, sum)/num_models * 100
tmp2$afr_ancestry <- plyr::mapvalues(tmp2$ancestry_assigned_by_study, c("African American","African Ethiopian", "African Ghanaian",
                                                                        "American Indian","Asian","Caucasian","Chinese", "Colombian",
                                                                        "European American", "Hispanic","Italian", "Kenyan","Mexican",
                                                                        "NA","Other","Taiwanese"), c("West_African","East_African",
                                                                                                     "West_African", "Non_African",
                                                                                                     "Non_African","Non_African",
                                                                                                     "Non_African","Non_African",
                                                                                                     "Non_African","Non_African",
                                                                                                     "Non_African","East_African",
                                                                                                     "Non_African","Non_African",
                                                                                                     "Non_African","Non_African"))
# afr_preds <- mega_validation_meta$sample_id[mega_validation_meta$train_data == "afr"]
#
# dat <- dat[dat$sample_id == mega_validation_meta$sample_id[mega_validation_meta$train_data == "afr"],]

afr_sub <- tidyr::gather(tmp2,key ="ancestry", "percent", PF_pct, BM_pct)


                                                                                                                                                                                                              "Non_African","Non_African"))
library(ggplot2)
## subset

no_na <- subset(afr_sub, !is.na(afr_ancestry))
no_na$afr_ancestry <- factor(no_na$afr_ancestry, levels = c("Non_African", "East_African", "West_African"))
ggplot(no_na) + geom_boxplot(aes(x=afr_ancestry, y=percent, fill=ancestry)) +
  #scale_fill_manual(values=c(European="firebrick2", African="dodger blue", Asian="forest green", American="goldenrod1", na.value="grey10")) +
  theme_minimal(base_size=24) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(fill='Continental Ancestry', y='% Models Correct', x='Race/Ethnicity Label') +
  geom_hline(aes(yintercept=median(tmp2$PF_pct))) +
  geom_hline(aes(yintercept=median(tmp2$BM_pct)), linetype = 2)
ggsave('Boxplot_PhyloFrame_Validation_Performance.png', width=10, height=7)
ggplot(tmp2) + geom_boxplot(aes(x=ancestry_assigned_by_study, y=BM_pct, fill=ancestry_assigned_by_us)) +
  scale_fill_manual(values=c(European="firebrick2", African="dodger blue", Asian="forest green", American="goldenrod1", na.value="grey10")) +
  theme_minimal(base_size=24) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(fill='Continental Ancestry', y='% Models Correct', x='Race/Ethnicity Label') +
  geom_hline(aes(yintercept=median(BM_pct)), linetype = 2) +
  geom_hline(aes(yintercept=median(PF_pct)))
ggsave('Boxplot_Benchmark_Validation_Performance.png', width=10, height=7)

t.test(tmp2$PF_pct, tmp2$BM_pct, paired = TRUE, alternative = "two.sided")
t.test(tmp2$PF_pct, tmp2$BM_pct, paired = TRUE, alternative = "greater")

t.test(tmp2[ tmp2$ancestry_assigned_by_us=='American','PF_pct'], tmp2[ tmp2$ancestry_assigned_by_us=='American','BM_pct'], paired = TRUE, alternative = "greater")
t.test(tmp2[ tmp2$ancestry_assigned_by_us=='Asian','PF_pct'], tmp2[ tmp2$ancestry_assigned_by_us=='Asian','BM_pct'], paired = TRUE, alternative = "greater")
t.test(tmp2[ tmp2$ancestry_assigned_by_us=='African','PF_pct'], tmp2[ tmp2$ancestry_assigned_by_us=='African','BM_pct'], paired = TRUE, alternative = "greater")
t.test(tmp2[ tmp2$ancestry_assigned_by_us=='European','PF_pct'], tmp2[ tmp2$ancestry_assigned_by_us=='European','BM_pct'], paired = TRUE, alternative = "greater")

median(tmp2$PF_pct)
median(tmp2$BM_pct)

#Specific comparisons of african
##
## create log file
tmp <- file.path(here::here("files"), "test.log")

t.test(tmp2[ tmp2$afr_ancestry=='East_African','PF_pct'], tmp2[ tmp2$afr_ancestry=='East_African','BM_pct'], paired = TRUE, alternative = "greater")
t.test(tmp2[ tmp2$ancestry_assigned_by_us=='European','PF_pct'], tmp2[ tmp2$ancestry_assigned_by_us=='European','BM_pct'], paired = TRUE, alternative = "greater")
t.test(tmp2[ tmp2$ancestry_assigned_by_us=='European','PF_pct'], tmp2[ tmp2$ancestry_assigned_by_us=='European','BM_pct'], paired = TRUE, alternative = "greater")




## Checking African ancestry
# East vs West: East = Kenyan, Ethiopian ; West = African American, Ghanaian
tmp4 <- tmp2[ tmp2ancestry_assigned_by_us=='African',]
tmp2$AFReVw <- tmp4$ancestry_assigned_by_study %in% c('Kenyan','Ethiopian')

# All non-AFR vs AFR
tmp2$nonAFR <- tmp2$ancestry_assigned_by_us=='European'
t.test(tmp2[ tmp2$nonAFR==TRUE,'PF_pct'], tmp2[ tmp2$nonAFR==TRUE,'BM_pct'], paired = TRUE, alternative = "greater")


# Make the barplot of continental ancestry vs subtype
tmp3 <- dat[,c('sample_id','ancestry_assigned_by_study','ancestry_assigned_by_us','subtype')]
tmp3$all <- paste0( tmp3$sample_id, tmp3$ancestry_assigned_by_study, tmp3$ancestry_assigned_by_us )
tmp3 <- tmp3[!duplicated(tmp3$all),]

ggplot(tmp3) + geom_bar(aes(x=ancestry_assigned_by_us, fill=subtype)) + geom_text(aes(label = after_stat(count), x=ancestry_assigned_by_us), stat = "count", vjust = 1.5, colour = "white") + scale_fill_manual(values=c("grey50", "grey70")) + theme_minimal()
ggsave('Barplot_ancestry_subtype.png', width=7,height=7)

## look at only EUR trained models (above is all models)
dat <- read.table('all_predictions_with_metadata.tsv', sep='\t', header=T)
dat <- dat[ dat$train_data=='eur', ]

dat$ancestry_assigned_by_us[ dat$ancestry_assigned_by_us=='None' ] <- 'Unknown'

dat$PF_pred_correct <- dat$pf_pred_class==dat$subtype
dat$BM_pred_correct <- dat$bm_pred_class==dat$subtype

dat2 <- dat[,c('sample_id','PF_pred_correct','BM_pred_correct','ancestry_assigned_by_study','ancestry_assigned_by_us','model_num')]
dat3 <- xtabs(PF_pred_correct ~ sample_id + model_num, dat2)
dat4 <- xtabs(BM_pred_correct ~ sample_id + model_num, dat2)
tmp <- dat2[,c('sample_id','ancestry_assigned_by_study','ancestry_assigned_by_us')]
tmp$all <- paste0( tmp$sample_id, tmp$ancestry_assigned_by_study, tmp$ancestry_assigned_by_us )
tmp2 <- tmp[!duplicated(tmp$all),]
num_models  <- length(unique(paste0(dat$train_data,dat$model_num)))
tmp2$PF_pct <- apply(dat3, 1, sum)/num_models * 100
tmp2$BM_pct <- apply(dat4, 1, sum)/num_models * 100

library(ggplot2)
ggplot(tmp2) + geom_boxplot(aes(x=ancestry_assigned_by_study, y=PF_pct, fill=ancestry_assigned_by_us)) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + geom_hline(aes(yintercept=median(PF_pct)))
ggplot(tmp2) + geom_boxplot(aes(x=ancestry_assigned_by_study, y=BM_pct, fill=ancestry_assigned_by_us)) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + geom_hline(aes(yintercept=median(BM_pct)))

#ggplot(tmp2) + geom_violin(aes(x=ancestry_assigned_by_study, y=PF_pct, colour=ancestry_assigned_by_us, fill=ancestry_assigned_by_us)) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + geom_hline(aes(yintercept=median(tmp2$PF_pct)))
#ggplot(tmp2) + geom_violin(aes(x=ancestry_assigned_by_study, y=BM_pct, colour=ancestry_assigned_by_us, fill=ancestry_assigned_by_us)) + theme(axis.text.x = element_text(angle = 45, hjust=1)) + geom_hline(aes(yintercept=median(tmp2$BM_pct)))

t.test(tmp2$PF_pct, tmp2$BM_pct, paired = TRUE, alternative = "greater")

median(tmp2$PF_pct)
median(tmp2$BM_pct)

# EUR = red, AFR = blue, EAS = green, AMR = yellow
scale_fill_manual(values=c(European="red", African="blue", Asian="green", American="yellow"))
