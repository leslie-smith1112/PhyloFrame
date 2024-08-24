## COSMIC gene heatmap

## SUPP FIG 3 ###

############ BRCA COSMIC HEATMAP  ##############

set.seed(123)
disease <- "breast"
dir <-  here::here("results","TCGA_Breast_Gnomad4_corrected")
ccgenes <- readr::read_tsv(here::here("Cosmic_Genes.tsv"))
genes <- ccgenes$`Gene Symbol`
############ EUR ############
## define lists for number of models in each ancestry for help reading in files ##
## BREAST model numbers
eur.num <- 1:17
afr.num <- 1:2
eas.num <- 1
all.num <- 20
admixed.num <- 1
mixed.num <- 1:6

eur.model.names <- paste0("eur_model_",eur.num)
afr.model.names <- paste0("afr_model_",afr.num)
eas.model.names <- paste0("eas_model_",eas.num)
admixed.names <- paste0("admixed_model_",admixed.num) # only one model
mixed.names <- paste0("mixed_model_", mixed.num)

## read in signatures from phyloframe and bechmark ancestry models
pf.eur.dir <- paste0(dir,"/model_runs/phyloFrame/eur/")
pf.afr.dir <- paste0(dir,"/model_runs/phyloFrame/afr/")
pf.admixed.dir <- paste0(dir,"/model_runs/phyloFrame/admix/")
pf.eas.dir <- paste0(dir,"/model_runs/phyloFrame/eas/")
pf.mixed.dir <- paste0(dir,"/model_runs/phyloFrame/mixed/")

eur.dir <- paste0(dir,"/model_runs/benchmark/eur/")
afr.dir <- paste0(dir,"/model_runs/benchmark/afr/")
admixed.dir <- paste0(dir,"/model_runs/benchmark/admix/")
eas.dir <- paste0(dir,"/model_runs/benchmark/eas/")
mixed.dir <- paste0(dir,"/model_runs/benchmark/mixed/")

pf.eur.names <- paste0(pf.eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt")
pf.eur.myfiles <- lapply(pf.eur.names, readr::read_tsv)
pf.afr.names <- paste0(pf.afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt")
pf.afr.myfiles <- lapply(pf.afr.names, readr::read_tsv)
pf.eas.names <- paste0(pf.eas.dir,"model_",eas.num,"/model_",eas.num,"_all_sig.txt")
pf.eas.myfiles <- lapply(pf.eas.names, readr::read_tsv)
pf.admixed.names <- paste0(pf.admixed.dir,"model_",admixed.num,"/model_",admixed.num,"_all_sig.txt")
pf.admixed.myfiles <- lapply(pf.admixed.names, readr::read_tsv)
pf.mixed.names <- paste0(pf.mixed.dir,"model_",mixed.num,"/model_",mixed.num,"_all_sig.txt")
pf.mixed.myfiles <- lapply(pf.mixed.names, readr::read_tsv)
#bind all PF signatures together
eur.pf <- do.call(rbind,pf.eur.myfiles)
afr.pf <- do.call(rbind,pf.afr.myfiles)
eas.pf <- do.call(rbind,pf.eas.myfiles)
admixed.pf <- do.call(rbind,pf.admixed.myfiles)
mixed.pf <- do.call(rbind,pf.mixed.myfiles)
all.pf <- rbind(eur.pf, afr.pf, eas.pf, admixed.pf, mixed.pf)

df.eur.names <- paste0(eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt")
eur.myfiles <- lapply(df.eur.names, readr::read_tsv)
df.afr.names <- paste0(afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt")
afr.myfiles <- lapply(df.afr.names, readr::read_tsv)
df.eas.names <- paste0(eas.dir,"model_",eas.num,"/model_",eas.num,"_all_sig.txt")
eas.myfiles <- lapply(df.eas.names, readr::read_tsv)
admixed.names <- paste0(admixed.dir,"model_",admixed.num,"/model_",admixed.num,"_all_sig.txt")
admixed.myfiles <- lapply(admixed.names, readr::read_tsv)
mixed.names <- paste0(mixed.dir,"model_",mixed.num,"/model_",mixed.num,"_all_sig.txt")
mixed.myfiles <- lapply(mixed.names, readr::read_tsv)
#bind all BM signatures together
eur.bm <- do.call(rbind,eur.myfiles)
afr.bm <- do.call(rbind,afr.myfiles)
eas.bm <- do.call(rbind,eas.myfiles)
admixed.bm <- do.call(rbind,admixed.myfiles)
mixed.bm <- do.call(rbind,mixed.myfiles)
all.bm <- rbind(eur.bm, afr.bm, eas.bm, admixed.bm, mixed.bm)

all.list <- c(pf.eur.myfiles,pf.afr.myfiles,pf.eas.myfiles,pf.admixed.myfiles,pf.mixed.myfiles,
              eur.myfiles,afr.myfiles,eas.myfiles,admixed.myfiles,mixed.myfiles)

allsies <- do.call(rbind,all.list)
in.genes <- allsies$Variable[allsies$Variable %in% genes]
in.genes <- unique(in.genes)
length(in.genes)

# ------------------- EUR ------------------- #
# data frame for phyloFrame and benchmark signtature percentages - 20 models total (17 eur, 2 afr, 1 eas)
eur.df <- data.frame(matrix(ncol = length(in.genes) + 1,nrow = 0))
colnames(eur.df) <- c("model",in.genes)
#loop through each signatures genes, if cosmic gene is in there mark TRUE else FALSE
#add row by row
for(i in 1:length(pf.eur.myfiles)){
  current <- pf.eur.myfiles[i][[1]]
  to.add <- in.genes %in% current$Variable
  to.add[to.add == TRUE] <- 1
  to.add[to.add == FALSE] <- 0
  # to.add <- as.numeric(to.add)
  to.add <- c(paste0("pfeurmodel",i), to.add)
  eur.df[nrow(eur.df) + 1,] <- to.add
}

for(i in 1:length(pf.afr.myfiles)){
  current <- pf.afr.myfiles[i][[1]]
  to.add <- in.genes %in% current$Variable
  to.add[to.add == TRUE] <- 2
  to.add[to.add == FALSE] <- 0
  # to.add <- as.numeric(to.add)
  to.add <- c(paste0("pfafrmodel",i), to.add)
  eur.df[nrow(eur.df) + 1,] <- to.add
}
for(i in 1:length(pf.eas.myfiles)){
  current <- pf.eas.myfiles[i][[1]]
  to.add <- in.genes %in% current$Variable
  to.add[to.add == TRUE] <- 3
  to.add[to.add == FALSE] <- 0
  # to.add <- as.numeric(to.add)
  to.add <- c(paste0("pfeasmodel",i), to.add)
  eur.df[nrow(eur.df) + 1,] <- to.add
}
for(i in 1:length(pf.admixed.myfiles)){
  current <- pf.admixed.myfiles[i][[1]]
  to.add <- in.genes %in% current$Variable
  to.add[to.add == TRUE] <- 4
  to.add[to.add == FALSE] <- 0
  # to.add <- as.numeric(to.add)
  to.add <- c(paste0("pfadmixedmodel",i), to.add)
  eur.df[nrow(eur.df) + 1,] <- to.add
}
# for(i in 1:length(pf.mixed.myfiles)){
#   current <- pf.mixed.myfiles[i][[1]]
#   to.add <- in.genes %in% current$Variable
#   to.add[to.add == TRUE] <- 5
#   to.add[to.add == FALSE] <- 0
#   # to.add <- as.numeric(to.add)
#   to.add <- c(paste0("pfmixedmodel",i), to.add)
#   eur.df[nrow(eur.df) + 1,] <- to.add
# }
library(tibble)
temp <- column_to_rownames(eur.df, "model")
temp <- as.matrix(temp)
temp2 <- apply(temp, 1, as.numeric)
rownames(temp2) <- colnames(temp)
pheatmap(temp2, color =  c("#eeeee9", "#CA4136","#1F619E","#496849","#FDC652"),
         cluster_rows = FALSE, cluster_cols = FALSE,fontsize = 8 )

gene_names <- data.frame("genes" = rownames(temp2))
write.table(gene_names, here::here("results", "files", "BRCA_cosmic_heatmap_genes.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

bmeur.df <- data.frame(matrix(ncol = length(in.genes) + 1,nrow = 0))
colnames(bmeur.df) <- c("model",in.genes)
for(i in 1:length(eur.myfiles)){
  current <- eur.myfiles[i][[1]]
  to.add <- in.genes %in% current$Variable
  to.add[to.add == TRUE] <- 1
  to.add[to.add == FALSE] <- 0
  # to.add <- as.numeric(to.add)
  to.add <- c(paste0("eurmodel",i), to.add)
  bmeur.df[nrow(bmeur.df) + 1,] <- to.add
}

for(i in 1:length(afr.myfiles)){
  current <- afr.myfiles[i][[1]]
  to.add <- in.genes %in% current$Variable
  to.add[to.add == TRUE] <- 2
  to.add[to.add == FALSE] <- 0
  # to.add <- as.numeric(to.add)
  to.add <- c(paste0("afrmodel",i), to.add)
  bmeur.df[nrow(bmeur.df) + 1,] <- to.add
}
for(i in 1:length(eas.myfiles)){
  current <- eas.myfiles[i][[1]]
  to.add <- in.genes %in% current$Variable
  to.add[to.add == TRUE] <- 3
  to.add[to.add == FALSE] <- 0
  # to.add <- as.numeric(to.add)
  to.add <- c(paste0("easmodel",i), to.add)
  bmeur.df[nrow(bmeur.df) + 1,] <- to.add
}
for(i in 1:length(admixed.myfiles)){
  current <- admixed.myfiles[i][[1]]
  to.add <- in.genes %in% current$Variable
  to.add[to.add == TRUE] <- 4
  to.add[to.add == FALSE] <- 0
  # to.add <- as.numeric(to.add)
  to.add <- c(paste0("admixedmodel",i), to.add)
  bmeur.df[nrow(bmeur.df) + 1,] <- to.add
}
# for(i in 1:length(mixed.myfiles)){
#   current <- mixed.myfiles[i][[1]]
#   to.add <- in.genes %in% current$Variable
#   to.add[to.add == TRUE] <- 5
#   to.add[to.add == FALSE] <- 0
#   # to.add <- as.numeric(to.add)
#   to.add <- c(paste0("mixedmodel",i), to.add)
#   bmeur.df[nrow(bmeur.df) + 1,] <- to.add
# }

temp.bm <- column_to_rownames(bmeur.df, "model")
temp.bm <- as.matrix(temp.bm)
temp.bm[1:5,1:5]

tempbm2 <- apply(temp.bm, 1, as.numeric)
rownames(tempbm2) <- colnames(temp.bm)
tempbm2[1:5,1:5]
pheatmap(tempbm2, color =  c("#eeeee9", "#CA4136","#1F619E","#496849","#FDC652"),
         cluster_rows = FALSE, cluster_cols = FALSE,fontsize = 8 )
