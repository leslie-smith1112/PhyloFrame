
  ####################################  READ IN ANCESTRY INFORMATION ####################################
  devtools::load_all() #load package to get associated data
  sample_dat <- readr::read_tsv(here::here("results","Multi_Study_Validation_corrected","model_runs","final_results.tsv"))
  colnames(sample_dat) <- c("sample_id", "subtype",".pred_class", ".pred_Basal",".pred_Luminal","train_data","model_num","version")
  head(sample_dat)
  sample_dat <- sample_dat[sample_dat$train_data != "mixed",]
  sample_dat <- sample_dat[,-9]
  sample_dat <- sample_dat[sample_dat$sample_id != "sample_id",]
  bm <- sample_dat[sample_dat$version == "benchmark",]
  pf <- sample_dat[sample_dat$version == "phyloFrame",]
  dim(pf)
  dim(bm)
  sample_list <- unique(sample_dat$sample_id)
  ## do once for phyloFrame then once for benchmark (below)
  ## PHYLOFRAME
  pf_dat <- merge(pf, all_vali_meta, by = "sample_id")
  pf_dat <- pf_dat |> dplyr::select(-subtype.y)
  colnames(pf_dat)[colnames(pf_dat) == "subtype.x"] <- "subtype"

  devtools::load_all()
  dat <- all_vali_preds[all_vali_preds$train_data != "mixed",]
  dim(dat)

  dat$ancestry_assigned_by_us[ dat$ancestry_assigned_by_us=='None' ] <- 'Unknown'
  dat$ancestry_assigned_by_study[ dat$ancestry_assigned_by_study=='Other' ] <- NA
  dat$PF_pred_correct <- dat$pf_pred_class==dat$subtype
  dat$BM_pred_correct <- dat$bm_pred_class==dat$subtype
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
  dat.all <- pivot_longer(tmp2, c('PF_pct','BM_pct'), names_to='ModelType', values_to='performance')
  head(dat.all)
  dat.all$axis <- paste0(dat.all$ancestry_assigned_by_us, dat.all$ModelType)
  dat.all$axis <- factor(dat.all$axis, levels = c("AfricanPF_pct","AfricanBM_pct","AmericanPF_pct","AmericanBM_pct","AsianPF_pct",
                                                  "AsianBM_pct","EuropeanPF_pct","EuropeanBM_pct","UnknownPF_pct","UnknownBM_pct"))


  p <- ggplot(dat.all, aes(x=axis, y=performance, fill = axis, color = axis)) + geom_boxplot() +
    scale_fill_manual(values = c("#1F619E", "#9fccfa","#FDC652","#ffe943","#496849","#63b34c","#CA4136","#f5a2a5", "#a986e6","#e1d0fe"
    )) + scale_color_manual(values = c("#153e64", "#7b9dbf","#ca9e43","#cab836","#2b3d2b","#427533","#7c2a23","#b97d7f", "#5b31a5","#26075c"))+ theme_minimal()
  p



#
#
#   pf_dat <- dat[dat$]
#   new_pf <- data.frame(matrix(ncol = 5, nrow = 0))
#   colnames(new_pf ) <- c("Sample_id",  "Continental_Ancestry", "SubContinental_Ancestry", "Model", "Percent_correct")
#   for(i in 1:length(sample_list)){
#     temp.dat <- pf_dat[pf_dat$sample_id == sample_list[i],] #get all predictions of sample
#     total <- nrow(temp.dat)
#     correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
#     percent <- (correct/total) * 100
#     new.list <- c(temp.dat$sample_id[1], temp.dat$continental_ancestry[1],temp.dat$subcontinental_ancestry[1],
#                   paste0(temp.dat$continental_ancestry[1],"PhyloFrame"),percent )
#     new_pf[nrow(new_pf) + 1,] <-new.list
#   }
#   new_pf$type <- "PhyloFrame"
#   #BENCHMARK
#   bm_dat <- merge(bm, all_vali_meta, by = "sample_id")
#   bm_dat <- bm_dat |> dplyr::select(-subtype.y)
#   new_bm <- data.frame(matrix(ncol = 5, nrow = 0))
#   colnames(new_bm ) <- c("Sample_id", "Continental_Ancestry", "SubContinental_Ancestry","Model", "Percent_correct")
#   for(i in 1:length(sample_list)){
#     temp.dat <- bm_dat[bm_dat$sample_id == sample_list[i],] #get all predictions of sample
#     total <- nrow(temp.dat)
#     correct <- nrow(temp.dat[temp.dat$subtype == temp.dat$.pred_class,])
#     percent <- (correct/total) * 100
#     new.list <- c(temp.dat$sample_id[1], temp.dat$continental_ancestry[1],temp.dat$subcontinental_ancestry[1],
#                   paste0(temp.dat$continental_ancestry[1],"Benchmark"),percent )
#     new_bm[nrow(new_bm) + 1,] <-new.list
#   }
#   new_bm$type <- "Benchmark"
#
#   all_samples <- rbind(new_pf, new_bm)
#   all_samples$Percent_correct <- as.numeric(all_samples$Percent_correct)
#   all_samples$Model <- factor(all_samples$Model, levels = c( "AfricanPhyloFrame", "AfricanBenchmark", ""AmericanBenchmark AmericanPhyloFrame,
#                                                              AsianBenchmark    AsianPhyloFrame  EuropeanBenchmark,
#                                                              EuropeanPhyloFrame      NoneBenchmark     NonePhyloFrame ))
#   p <- ggplot(all_samples, aes(x=Model, y=Percent_correct, fill = Model, color = Model)) + geom_boxplot() +
#       scale_fill_manual(values = c("#1F619E", "#9fccfa","#FDC652","#ffe943","#496849","#63b34c","#CA4136","#f5a2a5", "#a986e6","#e1d0fe"
#       )) + scale_color_manual(values = c("#153e64", "#7b9dbf","#ca9e43","#cab836","#2b3d2b","#427533","#7c2a23","#b97d7f", "#5b31a5","#26075c"))+ theme_minimal()
#     p
#
#
#   dat <- mega_validation_no_mixed
#   # for when we want ot subset to a specific model
#   #dat <- mega_validation_meta
#   #dat <- mega_validation_meta[mega_validation_meta$train_data != "mixed",]
#
#   dat$ancestry_assigned_by_us[ dat$ancestry_assigned_by_us=='None' ] <- 'Unknown'
#   dat$ancestry_assigned_by_study[ dat$ancestry_assigned_by_study=='Other' ] <- NA
#
#   dat$PF_pred_correct <- dat$pf_pred_class==dat$subtype
#   dat$BM_pred_correct <- dat$bm_pred_class==dat$subtype
#
#   dat$ancestry_assigned_by_us[ dat$ancestry_assigned_by_us=='None' ] <- 'Unknown'
#   dat$ancestry_assigned_by_study[ dat$ancestry_assigned_by_study=='Other' ] <- NA
#
#   dat$PF_pred_correct <- dat$pf_pred_class==dat$subtype
#   dat$BM_pred_correct <- dat$bm_pred_class==dat$subtype
#
#   dat2 <- dat[,c('sample_id','PF_pred_correct','BM_pred_correct','ancestry_assigned_by_study','ancestry_assigned_by_us','model_num')]
#   dat3 <- xtabs(PF_pred_correct ~ sample_id + model_num, dat2)
#   dat4 <- xtabs(BM_pred_correct ~ sample_id + model_num, dat2)
#   tmp <- dat2[,c('sample_id','ancestry_assigned_by_study','ancestry_assigned_by_us')]
#   tmp$all <- paste0( tmp$sample_id, tmp$ancestry_assigned_by_study, tmp$ancestry_assigned_by_us )
#   tmp2 <- tmp[!duplicated(tmp$all),]
#   tmp2 <- tmp2[tmp2$all != "NANANA",]
#
#   num_models  <- length(unique(paste0(dat$train_data,dat$model_num)))
#
#   tmp2$PF_pct <- apply(dat3, 1, sum)/num_models * 100
#   tmp2$BM_pct <- apply(dat4, 1, sum)/num_models * 100
#
#   afr_sub <- tidyr::gather(tmp2,key ="Ancestry", "Percent", PF_pct, BM_pct)
# library(ggplot2)
#   p <- ggplot(afr_sub, aes(x=ancestry_assigned_by_us , y=Percent, fill = Ancestry, color = Ancestry)) + geom_boxplot() +
#     scale_fill_manual(values = c("#1F619E", "#9fccfa","#FDC652","#ffe943","#496849","#63b34c","#CA4136","#f5a2a5", "#a986e6","#e1d0fe"
#     )) + scale_color_manual(values = c("#153e64", "#7b9dbf","#ca9e43","#cab836","#2b3d2b","#427533","#7c2a23","#b97d7f", "#5b31a5","#26075c"))+ theme_minimal()
#   p
#
#
#   #
# #   ggplot(all_vali_meta, aes(x = continental_ancestry, color = subtype, fill = subtype)) + geom_bar()
# #
# #   p<-ggplot(all_samples, aes(x=Model, y=Percent_correct, fill=Model, color = Model)) +
# #     geom_violin(position=position_dodge(1)) +
# #     scale_fill_manual(values = c("#1F619E", "#9fccfa","#FDC652","#ffe943","#496849","#63b34c","#CA4136","#f5a2a5", "#a986e6","#e1d0fe"
# #     )) + scale_color_manual(values = c("#153e64", "#7b9dbf","#ca9e43","#cab836","#2b3d2b","#427533","#7c2a23","#b97d7f", "#5b31a5","#26075c"))+ theme_minimal() +
# #     geom_violin(trim=TRUE,position=position_dodge(1))+
# #     geom_dotplot(binaxis='y', stackdir='center', dotsize=0.05)
# #     #geom_boxplot(width=0.1, fill="white")
# #   p
# # p
#   library(ggdist)
#
#   all_samples %>%
#     ggplot(aes(x = Ancestry, y = Percent_correct, fill = type)) +
#
#     # add half-violin from {ggdist} package
#     stat_halfeye(
#       # adjust bandwidth
#       adjust = 0.5,
#       # move to the right
#       justification = -0.2,
#       # remove the slub interval
#       .width = 0,
#       point_colour = NA
#     ) +
#     geom_boxplot(
#       width = 0.12,
#       # removing outliers
#       outlier.color = NA,
#       alpha = 0.5
#     ) +
#     stat_dots(
#       # ploting on left side
#       side = "left",
#       # adjusting position
#       justification = 1.1,
#       # adjust grouping (binning) of observations
#       binwidth = 0.25
#     )
#

