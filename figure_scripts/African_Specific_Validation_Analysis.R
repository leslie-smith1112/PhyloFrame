### African Validation From Kiley
#dat <- read.table('all_predictions_with_metadata.tsv', sep='\t', header=T)
dat <- all_vali_preds
dat <- dat[ dat$train_data!='mixed',] # We no longer consider the MIXED models - remove them

dat$ancestry_assigned_by_us[ dat$ancestry_assigned_by_us=='None' ] <- 'Unknown'
dat$ancestry_assigned_by_study[ dat$ancestry_assigned_by_study=='Other' ] <- NA

dat$PF_pred_correct <- dat$pf_pred_class==dat$subtype
dat$BM_pred_correct <- dat$bm_pred_class==dat$subtype

# Make metadata
dat.meta <- dat[ , c("sample_id","ancestry_assigned_by_study","subtype","ancestry_assigned_by_us")]
dat.meta <- dat.meta[!duplicated(dat.meta),]

dat2 <- dat[,c('sample_id','PF_pred_correct','BM_pred_correct','ancestry_assigned_by_study','ancestry_assigned_by_us','model_num')]
dat3 <- xtabs(PF_pred_correct ~ sample_id + model_num, dat2)
dat4 <- xtabs(BM_pred_correct ~ sample_id + model_num, dat2)
tmp <- dat2[,c('sample_id','ancestry_assigned_by_study','ancestry_assigned_by_us')]
tmp$all <- paste0( tmp$sample_id, tmp$ancestry_assigned_by_study, tmp$ancestry_assigned_by_us )
tmp2 <- tmp[!duplicated(tmp$all),]
num_models  <- length(unique(paste0(dat$train_data,dat$model_num)))
tmp2$PF_pct <- apply(dat3, 1, sum)/num_models * 100
tmp2$BM_pct <- apply(dat4, 1, sum)/num_models * 100


## Adding back the metadata to the per-sample performance info
rownames(tmp2) <- tmp2$sample_id
rownames(dat.meta) <- dat.meta$sample_id
dat.all <- merge(dat.meta, tmp2[,c('PF_pct','BM_pct')], by='row.names',all=TRUE)
rownames(dat.all) <- dat.all$Row.names
dat.all <- dat.all[,-1] # Remove the row.names ccolumn
dat.all$EvW <- 'a non-African'  # Adding the a in front to simplify my life. Should be releveling the data
dat.all[ dat.all$ancestry_assigned_by_study %in% c('Kenyan','African Ethiopian'), 'EvW' ] <- 'East African'
dat.all[ dat.all$ancestry_assigned_by_study %in% c('African American', 'African Ghanaian'), 'EvW' ] <- 'West African'
# East = Kenya, Ethiopia

sink(here::here("results","files","Validation_African_Comparisons.log"), split = TRUE) # open connection to log
## Run statistical tests before flattening data frame
t.test(dat.all[ dat.all$EvW=='East African', 'PF_pct'], dat.all[dat.all$EvW=='East African', 'BM_pct'], paired = TRUE, alternative = "greater")
t.test(dat.all[ dat.all$EvW=='West African', 'PF_pct'], dat.all[dat.all$EvW=='West African', 'BM_pct'], paired = TRUE, alternative = "greater")
t.test(dat.all[ dat.all$EvW=='a non-African', 'PF_pct'], dat.all[dat.all$EvW=='a non-African', 'BM_pct'], paired = TRUE, alternative = "greater")

# East vs West PF, repeat for BM
t.test(dat.all[ dat.all$EvW=='East African', 'PF_pct'], dat.all[dat.all$EvW=='West African', 'PF_pct'], paired = FALSE, alternative = "greater")
t.test(dat.all[ dat.all$EvW=='East African', 'BM_pct'], dat.all[dat.all$EvW=='West African', 'BM_pct'], paired = FALSE, alternative = "greater") # Just checking :)
t.test(dat.all[ dat.all$EvW=='a non-African', 'PF_pct'], dat.all[dat.all$EvW!='a non-African', 'PF_pct'], paired = FALSE, alternative = "greater")
sink() #close connections

library(tidyr)
dat.all <- pivot_longer(dat.all, c('PF_pct','BM_pct'), names_to='ModelType', values_to='performance')
head(dat.all)
dat.all$ModelType <- factor(dat.all$ModelType, levels = c("PF_pct", "BM_pct"))
# Plot by AFR ancestry region only
library(ggplot2)
ggplot(dat.all) + geom_boxplot(aes(x=EvW, y=performance, fill=ModelType, color = ModelType)) +
  scale_fill_manual(values=c(PF_pct="#1F619E", BM_pct="#9fccfa")) +
  scale_color_manual(values = c(PF_pct="#153e64", BM_pct="#7b9dbf")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(fill='Model Type', y='% Models Correct', x='West or East African')
ggsave('boxplot_validation_EastVsWestAfrican.png', width=10, height=7)


ggplot(as.data.frame(table(dat.all[dat.all$ModelType == "PF_pct",]$ancestry_assigned_by_us, dat.all[dat.all$ModelType == "PF_pct",]$subtype)), aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(position = "stack", stat = "identity")  + scale_fill_manual(values = c("#b993aa","#dfc6a9")) + theme_minimal()

#c("#b993aa","#e7bad5"))

ggplot(data, aes(fill=condition, y=value, x=specie)) +
  geom_bar(position="stack", stat="identity")











