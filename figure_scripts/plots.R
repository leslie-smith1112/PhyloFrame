
# Recall within each training data ancestry
sapply( unique(dat$Training_Data), function(x) {print(x); t.test(dat[ dat$Training_Data==x,'PF_Recall' ], dat[ dat$Training_Data==x,'BM_Recall'] ) } )
sapply( unique(dat$Training_Data), function(x) {print(x); t.test(dat[ dat$Training_Data==x,'PF_Mcc' ], dat[ dat$Training_Data==x,'BM_Mcc'] ) } )
sapply( unique(dat$Training_Data), function(x) {print(x); t.test(dat[ dat$Training_Data==x,'PF_Precision' ], dat[ dat$Training_Data==x,'BM_Precision'] ) } )

# Pull the data for each cancer, combine into 1 matrix
dat.brca <- read.table(here::here("results","TCGA_Breast_Gnomad4_corrected", "model_runs","batch_predictions","master_predictions.tsv"), sep='\t', header=T, stringsAsFactors=F)
dat.ucec <- read.table(here::here("results","TCGA_Uterine_Gnomad4_corrected", "model_runs","batch_predictions","master_predictions.tsv"), sep='\t', header=T, stringsAsFactors=F)
dat.thyr <- read.table(here::here("results","TCGA_Thyroid_Gnomad4_corrected", "model_runs","batch_predictions","master_predictions.tsv"), sep='\t', header=T, stringsAsFactors=F)


dat.brca$Cancer <- 'BRCA'
dat.ucec$Cancer <- 'UCEC'
dat.thyr$Cancer <- 'THCA'

## here train-data is the model training data, model_num is the model for the trainined data, ancestry_batch is the test data and batch_num is the batch for the test data
dat.brca <- dat.brca[(dat.brca$train_data != "mixed" & dat.brca$ancestry_batch != "mixed"),]
dat.ucec <- dat.ucec[(dat.ucec$train_data != "mixed" & dat.ucec$ancestry_batch != "mixed"),]
dat.thyr <- dat.thyr[(dat.thyr$train_data != "mixed" & dat.thyr$ancestry_batch != "mixed"),]


dat <- rbind( dat.thyr, dat.brca, dat.ucec)

# Plots for the supp figure
cols <- c(eur="#CA4136",mixed='peachpuff4',admix="#FDC652",eas="#496849",afr="#1F619E") # Keep colors same across cancer plots
library(ggplot2)
# scale_fill_manual(values = c("#1F619E", "#9fccfa","#FDC652","#ffe943","#496849","#63b34c",,"#f5a2a5", "#a986e6","#e1d0fe"
# )) + scale_color_manual(values = c("#153e64", "#7b9dbf","#ca9e43","#cab836","#2b3d2b","#427533","#7c2a23","#b97d7f", "#5b31a5","#26075c"))+ theme_bw()
# p

# BRCA
ggplot(dat.brca, aes(color=train_data,x=bm_recall, y=pf_recall)) + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) +geom_point(alpha = 0.6, size = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","BRCA_Recall.png"),plot = last_plot(),width = 7, height = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","BRCA_Recall.svg"),plot = last_plot(),width = 7, height = 5)

ggplot(dat.brca, aes(color=train_data,x=bm_mcc, y=pf_mcc)) + geom_point(alpha = 0.6, size = 5) + geom_abline(intercept=0) + lims( x=c(-1,1), y=c(-1,1) ) + scale_colour_manual(values=cols) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))
ggsave(here::here("figures","Metrics_Figure","most_recent","BRCA_Mcc.png"),plot = last_plot(),width = 7, height = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","BRCA_Mcc.svg"),plot = last_plot(),width = 7, height = 5)

ggplot(dat.brca, aes(color=train_data,x=bm_precision, y=pf_precision)) + geom_point(alpha = 0.6, size = 5) + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) ) + scale_colour_manual(values=cols) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) #+  geom_point();#ggsave('BRCA_precision.png',width=5.1,height=4)
ggsave(here::here("figures","Metrics_Figure","most_recent","BRCA_Precision.png"),plot = last_plot(),width = 7, height = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","BRCA_Precision.svg"),plot = last_plot(),width = 7, height = 5)

ggplot(dat.brca, aes(color=train_data,x=bm_auc, y=pf_auc)) + geom_point(alpha = 0.6, size = 5) + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")); #ggsave('BRCA_AUC.png',width=5.1,height=4)
ggsave(here::here("figures","Metrics_Figure","most_recent","BRCA_Auc.png"),plot = last_plot(),width = 7, height = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","BRCA_Auc.svg"),plot = last_plot(),width = 7, height = 5)

## - we decided not to include these
# ggplot(dat.brca, aes(color=train_data,x=bm_kapp, y=pf_kapp)) + geom_point() + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols) +
#   theme(axis.text=element_text(size=20),
#         axis.title=element_text(size=20,face="bold")) +  geom_point(alpha = 0.5, size = 5);# ggsave('BRCA_AUC.png',width=5.1,height=4)
# ggsave(here::here("figures","Metrics_Figure","BRCA_Kapp.png"),plot = last_plot(),width = 7, height = 5)
#
# ggplot(dat.brca, aes(color=train_data,x=bm_fmeas, y=pf_fmeas)) + geom_point() + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols) +
#   theme(axis.text=element_text(size=20),
#         axis.title=element_text(size=20,face="bold")) +  geom_point(alpha = 0.5, size = 5);# ggsave('BRCA_AUC.png',width=5.1,height=4)
# ggsave(here::here("figures","Metrics_Figure","BRCA_Fmeas.png"),plot = last_plot(),width = 7, height = 5)
#


# UCEC
ggplot(dat.ucec, aes(color=train_data,x=bm_recall, y=pf_recall)) + geom_point(alpha = 0.6, size = 5) + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) # ggsave('UCEC_recall.png',width=5.1,height=4)
ggsave(here::here("figures","Metrics_Figure","most_recent","UCEC_Recall.png"),plot = last_plot(),width = 7, height = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","UCEC_Recall.svg"),plot = last_plot(),width = 7, height = 5)

ggplot(dat.ucec, aes(color=train_data,x=bm_mcc, y=pf_mcc)) + geom_point(alpha = 0.6, size = 5) + geom_abline(intercept=0) + lims( x=c(-1,1), y=c(-1,1) )+ scale_colour_manual(values=cols) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) #ggsave('UCEC_MCC.png',width=5.1,height=4)
ggsave(here::here("figures","Metrics_Figure","most_recent","UCEC_Mcc.png"),plot = last_plot(),width = 7, height = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","UCEC_Mcc.svg"),plot = last_plot(),width = 7, height = 5)

ggplot(dat.ucec, aes(color=train_data,x=bm_precision, y=pf_precision)) + geom_point(alpha = 0.6, size = 5) + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))#ggsave('UCEC_precision.png',width=5.1,height=4)
ggsave(here::here("figures","Metrics_Figure","most_recent","UCEC_Precision.png"),plot = last_plot(),width = 7, height = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","UCEC_Precision.svg"),plot = last_plot(),width = 7, height = 5)

ggplot(dat.ucec, aes(color=train_data,x=bm_auc, y=pf_auc)) + geom_point(alpha = 0.6, size = 5) + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) #ggsave('UCEC_AUC.png',width=5.1,height=4)
ggsave(here::here("figures","Metrics_Figure","most_recent","UCEC_Auc.png"),plot = last_plot(),width = 7, height = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","UCEC_Auc.svg"),plot = last_plot(),width = 7, height = 5)


## - we decided not to include these
# ggplot(dat.ucec, aes(color=train_data,x=bm_kapp, y=pf_kapp)) + geom_point() + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols)+
#   theme(axis.text=element_text(size=20),
#         axis.title=element_text(size=20,face="bold")) +  geom_point(alpha = 0.5, size = 5)#ggsave('UCEC_precision.png',width=5.1,height=4)
# ggsave(here::here("figures","Metric_Figure","UCEC_Kapp.png"),plot = last_plot(),width = 7, height = 5)
#
# ggplot(dat.ucec, aes(color=train_data,x=bm_fmeas, y=pf_fmeas)) + geom_point() + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols) +
#   theme(axis.text=element_text(size=20),
#         axis.title=element_text(size=20,face="bold")) +  geom_point(alpha = 0.5, size = 5); #ggsave('UCEC_AUC.png',width=5.1,height=4)
# ggsave(here::here("figures","Metric_Figure","UCEC_Fmeas.png"),plot = last_plot(),width = 7, height = 5)

# THYR
ggplot(dat.thyr, aes(color=train_data,x=bm_recall, y=pf_recall)) + geom_point(alpha = 0.6, size = 5) + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols)+
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) #; ggsave('THYR_recall.png',width=5.1,height=4)
ggsave(here::here("figures","Metrics_Figure","most_recent","THCA_Recall.png"),plot = last_plot(),width = 7, height = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","THCA_Recall.svg"),plot = last_plot(),width = 7, height = 5)


ggplot(dat.thyr, aes(color=train_data,x=bm_mcc, y=pf_mcc)) + geom_point(alpha = 0.6, size = 5) + geom_abline(intercept=0) + lims( x=c(-1,1), y=c(-1,1) )+ scale_colour_manual(values=cols) +
  theme_bw() +
  theme(axis.text=element_text(size=30),
        axis.title=element_text(size=30,face="bold")) #; ggsave('THYR_MCC.png',width=5.1,height=4)
ggsave(here::here("figures","Metrics_Figure","most_recent","THCA_Mcc.png"),plot = last_plot(),width = 7, height = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","THCA_Mcc.svg"),plot = last_plot(),width = 7, height = 5)

ggplot(dat.thyr, aes(color=train_data,x=bm_precision, y=pf_precision)) + geom_point(alpha = 0.6, size = 5) + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) #; ggsave('THYR_precision.png',width=5.1,height=4)
ggsave(here::here("figures","Metrics_Figure","most_recent","THCA_Precision.png"),plot = last_plot(),width = 7, height = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","THCA_Precision.svg"),plot = last_plot(),width = 7, height = 5)


ggplot(dat.thyr, aes(color=train_data,x=bm_auc, y=pf_auc)) + geom_point(alpha = 0.6, size = 5) + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) #; ggsave('THYR_AUC.png',width=5.1,height=4)
ggsave(here::here("figures","Metrics_Figure","most_recent","THCA_Auc.png"),plot = last_plot(),width = 7, height = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","THCA_Auc.svg"),plot = last_plot(),width = 7, height = 5)


## - we decided not to include these
# ggplot(dat.thyr, aes(color=train_data,x=bm_kapp, y=pf_kapp)) + geom_point() + geom_abline(intercept=0) + lims( x=c(-1,1), y=c(-1,1) )+ scale_colour_manual(values=cols) +
#   theme(axis.text=element_text(size=20),
#         axis.title=element_text(size=20,face="bold")) +  geom_point(alpha = 0.5, size = 5)#; ggsave('THYR_precision.png',width=5.1,height=4)
# ggsave(here::here("figures","Metric_Figure","THCA_Kapp.png"),plot = last_plot(),width = 7, height = 5)
# ggplot(dat.thyr, aes(color=train_data,x=bm_fmeas, y=pf_fmeas)) + geom_point() + geom_abline(intercept=0) + lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols) +
#   theme(axis.text=element_text(size=20),
#         axis.title=element_text(size=20,face="bold")) +  geom_point(alpha = 0.5, size = 5)#; ggsave('THYR_AUC.png',width=5.1,height=4)
# ggsave(here::here("figures","Metric_Figure","THCA_Fmeas.png"),plot = last_plot(),width = 7, height = 5)


# Save t-test results to file
# All model metrics
sink(here::here("results","files","figure_stats.log"), append = TRUE, split = TRUE)
# Splitting by cancer type
cat("=============================\n")
cat("Metrics - All models and all ancestries, split by cancer type\n")
cat("=============================\n")
cat('Recall'); sapply( unique(dat$Cancer), function(x) {print(x); t.test(dat[ dat$Cancer==x,'pf_recall' ], dat[ dat$Cancer==x,'bm_recall'], alternative = "greater", paired = TRUE )$p.value } )
cat('MCC'); sapply( unique(dat$Cancer), function(x) {print(x); t.test(dat[ dat$Cancer==x,'pf_mcc' ], dat[ dat$Cancer==x,'bm_mcc'], alternative = "greater", paired = TRUE  )$p.value } )
#breast is constant
cat('Precision'); sapply( c("THYR","UCEC"), function(x) {print(x); t.test(dat[ dat$Cancer==x,'pf_precision' ], dat[ dat$Cancer==x,'bm_precision'], alternative = "greater", paired = TRUE  )$p.value } )
cat('AUC'); sapply( unique(dat$Cancer), function(x) {print(x); t.test(dat[ dat$Cancer==x,'pf_auc' ], dat[ dat$Cancer==x,'bm_auc'], alternative = "greater", paired = TRUE  )$p.value } )

sink()


cat("=============================\n")
cat("Metrics across all models, all cancers and ancestries together\n")
cat("=============================\n")
t.test( dat$PF_Recall, dat$BM_Recall ) # Recall across ALL tests
t.test( dat$PF_Precision, dat$BM_Precision ) # Precision across ALL tests
t.test( dat$PF_Mcc, dat$BM_Mcc ) # Mcc across ALL tests
t.test( dat$PF_Auc, dat$BM_Auc ) # Auc across ALL tests

# Splitting by ancestries
cat("=============================\n")
cat("Metrics - All models and all cancers, split by ancestry\n")
cat("=============================\n")
cat('Recall'); sapply( unique(dat$Training_Data), function(x) {print(x); t.test(dat[ dat$Training_Data==x,'PF_Recall' ], dat[ dat$Training_Data==x,'BM_Recall'] )$p.value } )
cat('MCC'); sapply( unique(dat$Training_Data), function(x) {print(x); t.test(dat[ dat$Training_Data==x,'PF_Mcc' ], dat[ dat$Training_Data==x,'BM_Mcc'] )$p.value} )
cat('Precision'); sapply( unique(dat$Training_Data), function(x) {print(x); t.test(dat[ dat$Training_Data==x,'PF_Precision' ], dat[ dat$Training_Data==x,'BM_Precision'] )$p.value } )
cat('AUC'); sapply( unique(dat$Training_Data), function(x) {print(x); t.test(dat[ dat$Training_Data==x,'PF_Auc' ], dat[ dat$Training_Data==x,'BM_Auc'] )$p.value } )






sink('t.tests_split_by_cancer_type.txt')

# Ancestries together - for supp figure
cat("=============================\n")
cat("BRCA metrics across all models, all ancestries together\n")
cat("=============================\n")
t.test( dat.brca$PF_Recall, dat.brca$BM_Recall ) # Recall across ALL tests
t.test( dat.brca$PF_Precision, dat.brca$BM_Precision ) # Precision across ALL tests
t.test( dat.brca$PF_Mcc, dat.brca$BM_Mcc ) # Mcc across ALL tests
t.test( dat.brca$PF_Auc, dat.brca$BM_Auc ) # Auc across ALL tests

cat("=============================\n")
cat("THYR metrics across all models, all ancestries together\n")
cat("=============================\n")
t.test( dat.thyr$PF_Recall, dat.thyr$BM_Recall ) # Recall across ALL tests
t.test( dat.thyr$PF_Precision, dat.thyr$BM_Precision ) # Precision across ALL tests
t.test( dat.thyr$PF_Mcc, dat.thyr$BM_Mcc ) # Mcc across ALL tests
t.test( dat.thyr$PF_Auc, dat.thyr$BM_Auc ) # Auc across ALL tests

cat("=============================\n")
cat("UCEC metrics across all models, all ancestries together\n")
cat("=============================\n")
t.test( dat.ucec$PF_Recall, dat.ucec$BM_Recall ) # Recall across ALL tests
t.test( dat.ucec$PF_Precision, dat.ucec$BM_Precision ) # Precision across ALL tests
t.test( dat.ucec$PF_Mcc, dat.ucec$BM_Mcc ) # Mcc across ALL tests
t.test( dat.ucec$PF_Auc, dat.ucec$BM_Auc ) # Auc across ALL tests


# Splitting by ancestries
cat("=============================\n")
cat("BRCA - All models, split by ancestry\n")
cat("=============================\n")
cat('Recall'); sapply( unique(dat.brca$Training_Data), function(x) {print(x); t.test(dat.brca[ dat.brca$Training_Data==x,'PF_Recall' ], dat.brca[ dat.brca$Training_Data==x,'BM_Recall'] )$p.value } )
cat('MCC'); sapply( unique(dat.brca$Training_Data), function(x) {print(x); t.test(dat.brca[ dat.brca$Training_Data==x,'PF_Mcc' ], dat.brca[ dat.brca$Training_Data==x,'BM_Mcc'] )$p.value} )
cat('Precision'); sapply( unique(dat.brca$Training_Data), function(x) {print(x); t.test(dat.brca[ dat.brca$Training_Data==x,'PF_Precision' ], dat.brca[ dat.brca$Training_Data==x,'BM_Precision'] )$p.value } )
cat('AUC'); sapply( unique(dat.brca$Training_Data), function(x) {print(x); t.test(dat.brca[ dat.brca$Training_Data==x,'PF_Auc' ], dat.brca[ dat.brca$Training_Data==x,'BM_Auc'] )$p.value } )

cat("=============================\n")
cat("THYR - All models, split by ancestry\n")
cat("=============================\n")
cat('Recall'); sapply( unique(dat.thyr$Training_Data), function(x) {print(x); t.test(dat.thyr[ dat.thyr$Training_Data==x,'PF_Recall' ], dat.thyr[ dat.thyr$Training_Data==x,'BM_Recall'] )$p.value } )
cat('MCC'); sapply( unique(dat.thyr$Training_Data), function(x) {print(x); t.test(dat.thyr[ dat.thyr$Training_Data==x,'PF_Mcc' ], dat.thyr[ dat.thyr$Training_Data==x,'BM_Mcc'] )$p.value} )
cat('Precision'); sapply( unique(dat.thyr$Training_Data), function(x) {print(x); t.test(dat.thyr[ dat.thyr$Training_Data==x,'PF_Precision' ], dat.thyr[ dat.thyr$Training_Data==x,'BM_Precision'] )$p.value } )
cat('AUC'); sapply( unique(dat.thyr$Training_Data), function(x) {print(x); t.test(dat.thyr[ dat.thyr$Training_Data==x,'PF_Auc' ], dat.thyr[ dat.thyr$Training_Data==x,'BM_Auc'] )$p.value } )

cat("=============================\n")
cat("UCEC - All models, split by ancestry\n")
cat("=============================\n")
cat('Recall'); sapply( unique(dat.ucec$Training_Data), function(x) {print(x); t.test(dat.ucec[ dat.ucec$Training_Data==x,'PF_Recall' ], dat.ucec[ dat.ucec$Training_Data==x,'BM_Recall'] )$p.value } )
cat('MCC'); sapply( unique(dat.ucec$Training_Data), function(x) {print(x); t.test(dat.ucec[ dat.ucec$Training_Data==x,'PF_Mcc' ], dat.ucec[ dat.ucec$Training_Data==x,'BM_Mcc'] )$p.value} )
cat('Precision'); sapply( unique(dat.ucec$Training_Data), function(x) {print(x); t.test(dat.ucec[ dat.ucec$Training_Data==x,'PF_Precision' ], dat.ucec[ dat.ucec$Training_Data==x,'BM_Precision'] )$p.value } )
cat('AUC'); sapply( unique(dat.ucec$Training_Data), function(x) {print(x); t.test(dat.ucec[ dat.ucec$Training_Data==x,'PF_Auc' ], dat.ucec[ dat.ucec$Training_Data==x,'BM_Auc'] )$p.value } )

sink()
dim(dat.thyr[!(is.na(dat.thyr$pf_precision) | is.na(dat.thyr$bm_precision)),])
#checking
dat.thyr[is.na(dat.thyr$pf_precision) | is.na(dat.thyr$bm_precision),]
model_5            eur    batch4
thca_preds <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/PhyloFrame/results/TCGA_Thyroid_Gnomad4_corrected/model_runs/sample_batches_and_predictions.tsv")

eur_model5 <- thca_preds[thca_preds$model_train_data == "eur" & thca_preds$model_num == "model_5" &
                           thca_preds$sample_batch_num  == "batch9" & thca_preds$sample_ancestry_batch == "eur",]

eur_model5$subtype <- as.factor(eur_model5$subtype)
eur_model5$.pred_class <- as.factor(eur_model5$.pred_class)
eur_model5 <- eur_model5[eur_model5$version == "benchmark",]
yardstick::conf_mat(eur_model5, "subtype",".pred_class")
####
dat.brca[is.na(dat.brca$pf_precision) | is.na(dat.brca$bm_precision),]
#model_5            eur    batch4
thca_preds <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/PhyloFrame/results/TCGA_Breast_Gnomad4_corrected/model_runs/sample_batches_and_predictions.tsv")

eur_model5 <- thca_preds[thca_preds$model_train_data == "eur" & thca_preds$model_num == "model_14" &
                           thca_preds$sample_batch_num  == "batch8" & thca_preds$sample_ancestry_batch == "eur",]
eur_model5 <- eur_model5[eur_model5$version == "phyloFrame",]

eur_model5$subtype <- as.factor(eur_model5$subtype)
eur_model5$.pred_class <- as.factor(eur_model5$.pred_class)
yardstick::conf_mat(eur_model5, "subtype",".pred_class")
##


dim(eur_model5)
eur_model5$subtype <- as.factor(eur_model5$subtype)
eur_model5$.pred_class <- as.factor(eur_model5$.pred_class)
yardstick::conf_mat(eur_model5, "subtype",".pred_class")

mcc
dim(dat.ucec[!(is.na(dat.ucec$pf_recall) | is.na(dat.ucec$bm_recall)),])

length(dat.ucec$pf_auc)
dim(dat.brca[!(is.na(dat.brca$pf_precision) | is.na(dat.brca$bm_precision)),])
## investigate
preds <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/PhyloFrame/results/TCGA_Breast_Gnomad4_corrected/model_runs/sample_batches_and_predictions.tsv")
eur_mod14 <- preds[preds$model_train_data == "eur" & preds$model_num == "model_14" & preds$sample_batch_num == "batch8",]
TP = nrow(eur_mod14[eur_mod14$subtype == "Basal" &  eur_mod14$.pred_class == "Basal",])
TN = nrow(eur_mod14[eur_mod14$subtype == "Luminal" &  eur_mod14$.pred_class == "Luminal",])
FP = nrow(eur_mod14[eur_mod14$subtype == "Luminal" &  eur_mod14$.pred_class == "Basal",])
FN = nrow(eur_mod14[eur_mod14$subtype == "Basal" &  eur_mod14$.pred_class == "Luminal",])
FP
MCC = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
MCC = (TP * TN - FP * FN) / sqrt((3) * (TP + FN) * (66) * (TN + FN))

library(tidyverse)
library(stargazer)
library(patchwork)
model <- lm(pf_mcc ~ bm_mcc, data = dat.brca)
cols <- c(eur="#CA4136",mixed='peachpuff4',admix="#FDC652",eas="#496849",afr="#1F619E") # Keep colors same across cancer plots

p1 <- dat.ucec %>%
  ggplot(aes( x = bm_recall, y = pf_recall)) +
  #scale_colour_manual(values=cols) +
  #scale_fill_manual(values = cols) +
  geom_point(aes(color = train_data),alpha = 0.6, size = 5) +
  geom_smooth(color = "grey20", fill = "grey40",method = "lm", se=TRUE, linewidth = 2) +
  geom_abline(linetype = 2, intercept=0)+
  lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols)+
  labs(title = "95% Confidence Interval") +
  theme_bw() +
  #theme(plot.title = element_text(face = "bold",hjust = 0.5))
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))
p1
ggsave(here::here("figures","Metrics_Figure","most_recent","UCEC_RECALL_CI_LM.png"),plot = last_plot(),width = 7, height = 5)
ggsave(here::here("figures","Metrics_Figure","most_recent","UCEC_RECALL_CI_LM.svg"),plot = last_plot(),width = 7, height = 5)


dim(dat.brca)



p1 <- dat.brca %>%
  ggplot(aes( x = bm_auc, y = pf_auc)) +
  #scale_colour_manual(values=cols) +
  #scale_fill_manual(values = cols) +
  geom_point(aes(color = train_data),alpha = 0.6, size = 5) +
  geom_smooth(color = "grey20", fill = "grey40",method = loess, linewidth = 2, span = 0.9) +
  geom_abline(linetype = 2, intercept=0)+
  lims( x=c(0,1), y=c(0,1) )+ scale_colour_manual(values=cols)+
  labs(title = "95% Confidence Interval") +
  theme_bw() +
  #theme(plot.title = element_text(face = "bold",hjust = 0.5))
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"))
p1



supp_table <- data.frame("Cancer" = dat$Cancer, "Model" = paste(dat$train_data,dat$model_num, sep = "_"), "Training_Data_Ancestry" = dat$train_data,
                         "Test_Data" = paste(dat$ancestry_batch,dat$batch_num, sep = "_"), "Test_Data_Ancestry" = dat$ancestry_batch,
                         "PF_ACCURACY"= dat$pf_accuracy, "BM_ACCURACY"= dat$bm_accuracy,
                         "PF_AUC" = dat$pf_auc,"BM_AUC" = dat$bm_auc, "PF_COHENSKAPPA"= dat$pf_kapp, "BM_COHENSKAPPA" = dat$bm_kapp,"PF_FMEAS"= dat$pf_fmeas,
                         "BM_FMEAS" = dat$bm_fmeas,"PF_MCC" = dat$pf_mcc,"BM_MCC" = dat$bm_mcc,"PF_PRECISION" =dat$pf_precision,"BM_PRECISION" = dat$bm_precision,
                         "PF_RECALL" = dat$pf_recall,"BM_RECALL" = dat$bm_recall)
