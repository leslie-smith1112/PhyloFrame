## Signature Signature Correlation Plot ##
## Only Implemented for BRCA ##
library(ComplexHeatmap)
set.seed(123)

## - Results directory naems:
TCGA_Breast_Gnomad4_range2


result_dir <- "TCGA_Breast_Gnomad4_corrected"
disease <- "breast"
############ EUR ############
## define lists for number of models in each ancestry for help reading in files ##

if(disease == "breast"){ #TODO this is the worst possible way to do this - needs
  ## BREAST
  eur.num <- 1:17
  afr.num <- 1:2
  eas.num <- 1
  all.num <- 20
  admix.num <- 1
}else if(disease == "thyroid"){
  ## THYROID:
  eur.num <- 1:23
  afr.num <- 1
  eas.num <- 1:3
  all.num <- 27
  admix.num <- 1
}else if(disease == "uterine"){
  ## UTERINE
  eur.num <- 1:12
  afr.num <- 1:2
  eas.num <- 1
  all.num <- 15
  admix.num <- 1
}else{
  stop("must enter valid disease")
}
# get number of models for each ancestry to help read in signature results
eur.count <- length(eur.num)
afr.count <- length(afr.num)
eas.count <- length(eas.num)
admix.count <- length(admix.num)

eur.model.names <- paste0("eur_model_",eur.num)
afr.model.names <- paste0("afr_model_",afr.num)
eas.model.names <- paste0("eas_model_",eas.num)
admix.names <- paste0("admix_model_",admix.num) # only one model

## read in signatures from phyloframe and bechmark ancestry models
pf.eur.dir <- paste0(here::here("results",result_dir), "/model_runs/phyloFrame/eur/")
pf.afr.dir <- paste0(here::here("results",result_dir), "/model_runs/phyloFrame/afr/")
pf.admix.dir <- paste0(here::here("results",result_dir), "/model_runs/phyloFrame/admix/")
pf.eas.dir <- paste0(here::here("results",result_dir), "/model_runs/phyloFrame/eas/")
pf.admix.dir <- paste0(here::here("results",result_dir), "/model_runs/phyloFrame/admix/")

eur.dir <- paste0(here::here("results",result_dir), "/model_runs/benchmark/eur/")
afr.dir <- paste0(here::here("results",result_dir), "/model_runs/benchmark/afr/")
admix.dir <- paste0(here::here("results",result_dir), "/model_runs/benchmark/admix/")
eas.dir <- paste0(here::here("results",result_dir), "/model_runs/benchmark/eas/")
admix.dir <- paste0(here::here("results",result_dir), "/model_runs/benchmark/admix/")

pf.eur.names <- paste0(pf.eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt")
pf.eur.myfiles <- lapply(pf.eur.names, readr::read_tsv)
pf.afr.names <- paste0(pf.afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt")
pf.afr.myfiles <- lapply(pf.afr.names, readr::read_tsv)
pf.eas.names <- paste0(pf.eas.dir,"model_",eas.num,"/model_",eas.num,"_all_sig.txt")
pf.eas.myfiles <- lapply(pf.eas.names, readr::read_tsv)
pf.admix.names <- paste0(pf.admix.dir,"model_",admix.num,"/model_",admix.num,"_all_sig.txt")
pf.admix.myfiles <- lapply(pf.admix.names, readr::read_tsv)


df.eur.names <- paste0(eur.dir,"model_",eur.num,"/model_",eur.num,"_all_sig.txt")
eur.myfiles <- lapply(df.eur.names, readr::read_tsv)
df.afr.names <- paste0(afr.dir,"model_",afr.num,"/model_",afr.num,"_all_sig.txt")
afr.myfiles <- lapply(df.afr.names, readr::read_tsv)
df.eas.names <- paste0(eas.dir,"model_",eas.num,"/model_",eas.num,"_all_sig.txt")
eas.myfiles <- lapply(df.eas.names, readr::read_tsv)
df.admix.names <- paste0(admix.dir,"model_",admix.num,"/model_",admix.num,"_all_sig.txt")
admix.myfiles <- lapply(df.admix.names, readr::read_tsv)
### SIGNAURE SIGNATURE CORRELATION ##
# ------------------- EUR ------------------- #
# data frame for phyloFrame and benchmark signtature percentages - 20 models total (17 eur, 2 afr, 1 eas)
eur.df <- data.frame(matrix(ncol = all.num,nrow = 0))
pf.eur.df <- data.frame(matrix(ncol = all.num,nrow = 0))
#for every ancestry signature from every model - divide intersection/union with signatures * 100 to get
#% overlap of each signature and build matrix
# do for each ancestry model
#*.myfiles is the list of the ancestry model signatures
for(i in 1:eur.count){
  temp <- data.frame(matrix(ncol = all.num,nrow = 0))
  pf.temp <- data.frame(matrix(ncol = all.num,nrow = 0))
  for(j in 1:eur.count){
    ### check the intersection of model i with model j
    the.int <- length(intersect(eur.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable))
    ### check how many genes total from the 2 models
    gene.union <- union(eur.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j] <- overlap
    ## phyloframe ##
    pf.the.int <- length(intersect(pf.eur.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eur.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j] <- pf.overlap
  }
  for(j in 1:afr.count){
    the.int <- length(intersect(eur.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable))
    gene.union <- union(eur.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j+eur.count] <- overlap
    ## phyloframe ##
    pf.the.int <- length(intersect(pf.eur.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eur.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j+eur.count] <- pf.overlap
  }
  for(j in 1:eas.count){
    the.int <- length(intersect(eur.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable))
    gene.union <- union(eur.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j+ (eur.count + afr.count)] <- overlap
    ## phyloframe ##
    pf.the.int <- length(intersect(pf.eur.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eur.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j + (eur.count + afr.count) ] <- pf.overlap
  }
  # for(j in 1:admix.count){
  #   the.int <- length(intersect(eur.myfiles[i][[1]]$Variable,admix.myfiles[j][[1]]$Variable))
  #   gene.union <- union(eur.myfiles[i][[1]]$Variable,admix.myfiles[j][[1]]$Variable)
  #   denom <- length(gene.union)
  #   overlap <- (the.int/denom) * 100
  #   temp[1,j+ (eur.count + afr.count + eas.count)] <- overlap
  #   ## phyloframe ##
  #   pf.the.int <- length(intersect(pf.eur.myfiles[i][[1]]$Variable,pf.admix.myfiles[j][[1]]$Variable))
  #   pf.gene.union <- union(pf.eur.myfiles[i][[1]]$Variable,pf.admix.myfiles[j][[1]]$Variable)
  #   pf.denom <- length(pf.gene.union)
  #   pf.overlap <- (pf.the.int/pf.denom) * 100
  #   pf.temp[1,j + (eur.count + afr.count + eas.count) ] <- pf.overlap
  # }
  eur.df <- rbind(eur.df, temp)
  pf.eur.df <- rbind(pf.eur.df, pf.temp)
}

rownames(eur.df) <- eur.model.names
rownames(pf.eur.df) <- eur.model.names
rownames(eur.df) <- paste0(rownames(eur.df))
rownames(pf.eur.df) <- paste0(rownames(pf.eur.df))
c.name <- c(eur.model.names, afr.model.names,eas.model.names)
colnames(eur.df) <- c.name
colnames(pf.eur.df) <- c.name

afr.df <- data.frame(matrix(ncol = all.num,nrow = 0))
pf.afr.df <- data.frame(matrix(ncol = all.num,nrow = 0))
#### AFR
for(i in 1:afr.count){
  temp <- data.frame(matrix(ncol = all.num,nrow = 0))
  pf.temp <- data.frame(matrix(ncol = all.num,nrow = 0))
  for(j in 1:eur.count){
    the.int <- length(intersect(afr.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable))
    gene.union <- union(afr.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j] <- overlap
    ## phyloframe
    pf.the.int <- length(intersect(pf.afr.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.afr.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j] <- pf.overlap
  }
  for(j in 1:afr.count){
    the.int <- length(intersect(afr.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable))
    gene.union <- union(afr.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j+ eur.count] <- overlap
    ## phyloFrame
    pf.the.int <- length(intersect(pf.afr.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.afr.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j+eur.count] <- pf.overlap
  }
  for(j in 1:eas.count){
    the.int <- length(intersect(afr.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable))
    gene.union <- union(afr.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j+ (eur.count + afr.count)] <- overlap
    ## phyloframe ##
    pf.the.int <- length(intersect(pf.afr.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.afr.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j + (eur.count + afr.count) ] <- pf.overlap
  }
  # for(j in 1:admix.count){
  #   the.int <- length(intersect(afr.myfiles[i][[1]]$Variable,admix.myfiles[j][[1]]$Variable))
  #   gene.union <- union(afr.myfiles[i][[1]]$Variable,admix.myfiles[j][[1]]$Variable)
  #   denom <- length(gene.union)
  #   overlap <- (the.int/denom) * 100
  #   temp[1,j+ (eur.count + afr.count + eas.count)] <- overlap
  #   ## phyloframe ##
  #   pf.the.int <- length(intersect(pf.afr.myfiles[i][[1]]$Variable,pf.admix.myfiles[j][[1]]$Variable))
  #   pf.gene.union <- union(pf.afr.myfiles[i][[1]]$Variable,pf.admix.myfiles[j][[1]]$Variable)
  #   pf.denom <- length(pf.gene.union)
  #   pf.overlap <- (pf.the.int/pf.denom) * 100
  #   pf.temp[1,j + (eur.count + afr.count + eas.count) ] <- pf.overlap
  # }
  afr.df <- rbind(afr.df, temp)
  pf.afr.df <- rbind(pf.afr.df, pf.temp)
}
rownames(afr.df) <- afr.model.names
rownames(pf.afr.df) <- afr.model.names
rownames(afr.df) <- paste0(rownames(afr.df))
rownames(pf.afr.df) <- paste0(rownames(pf.afr.df))
colnames(afr.df) <- c.name
colnames(pf.afr.df) <- c.name

#### EAS
eas.df <- data.frame(matrix(ncol = all.num,nrow = 0))
pf.eas.df <- data.frame(matrix(ncol = all.num,nrow = 0))

for(i in 1:eas.count){
  temp <- data.frame(matrix(ncol = all.num,nrow = 0))
  pf.temp <- data.frame(matrix(ncol = all.num,nrow = 0))
  for(j in 1:eur.count){
    ### checking model1 of phyloframe against all other models of benchamrk
    the.int <- length(intersect(eas.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable))
    gene.union <- union(eas.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j] <- overlap
    ## phyloFrame ##
    pf.the.int <- length(intersect(pf.eas.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eas.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j] <- pf.overlap
  }
  for(j in 1:afr.count){
    ### checking model1 of phyloframe against all other models of benchamrk
    the.int <- length(intersect(eas.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable))
    gene.union <- union(eas.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j+eur.count] <- overlap
    # phyloframe
    pf.the.int <- length(intersect(pf.eas.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eas.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j+ eur.count] <- pf.overlap
  }
  for(j in 1:eas.count){
    ### checking model1 of phyloframe against all other models of benchamrk
    the.int <- length(intersect(eas.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable))
    gene.union <- union(eas.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable)
    denom <- length(gene.union)
    overlap <- (the.int/denom) * 100
    temp[1,j+ (eur.count + afr.count)] <- overlap
    ## phyloframe ##
    pf.the.int <- length(intersect(pf.eas.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable))
    pf.gene.union <- union(pf.eas.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable)
    pf.denom <- length(pf.gene.union)
    pf.overlap <- (pf.the.int/pf.denom) * 100
    pf.temp[1,j + (eur.count + afr.count) ] <- pf.overlap
  }
  # for(j in 1:admix.count){
  #   the.int <- length(intersect(eas.myfiles[i][[1]]$Variable,admix.myfiles[j][[1]]$Variable))
  #   gene.union <- union(eas.myfiles[i][[1]]$Variable,admix.myfiles[j][[1]]$Variable)
  #   denom <- length(gene.union)
  #   overlap <- (the.int/denom) * 100
  #   temp[1,j+ (eur.count + afr.count + eas.count)] <- overlap
  #   ## phyloframe ##
  #   pf.the.int <- length(intersect(pf.eas.myfiles[i][[1]]$Variable,pf.admix.myfiles[j][[1]]$Variable))
  #   pf.gene.union <- union(pf.eas.myfiles[i][[1]]$Variable,pf.admix.myfiles[j][[1]]$Variable)
  #   pf.denom <- length(pf.gene.union)
  #   pf.overlap <- (pf.the.int/pf.denom) * 100
  #   pf.temp[1,j + (eur.count + afr.count + eas.count) ] <- pf.overlap
  # }
  eas.df <- rbind(eas.df, temp)
  pf.eas.df <- rbind(pf.eas.df, pf.temp)
}
#
rownames(eas.df) <- eas.model.names
rownames(pf.eas.df) <- eas.model.names
colnames(eas.df) <- c.name
colnames(pf.eas.df) <- c.name

#### admix
# admix.df <- data.frame(matrix(ncol = all.num,nrow = 0))
# pf.admix.df <- data.frame(matrix(ncol = all.num,nrow = 0))
#
# for(i in 1:eur.count){
#   temp <- data.frame(matrix(ncol = all.num,nrow = 0))
#   pf.temp <- data.frame(matrix(ncol = all.num,nrow = 0))
#   for(j in 1:eur.count){
#     ### checking model1 of phyloframe against all other models of benchamrk
#     the.int <- length(intersect(admix.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable))
#     gene.union <- union(admix.myfiles[i][[1]]$Variable,eur.myfiles[j][[1]]$Variable)
#     denom <- length(gene.union)
#     overlap <- (the.int/denom) * 100
#     temp[1,j] <- overlap
#     ## phyloFrame ##
#     pf.the.int <- length(intersect(pf.admix.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable))
#     pf.gene.union <- union(pf.admix.myfiles[i][[1]]$Variable,pf.eur.myfiles[j][[1]]$Variable)
#     pf.denom <- length(pf.gene.union)
#     pf.overlap <- (pf.the.int/pf.denom) * 100
#     pf.temp[1,j] <- pf.overlap
#   }
#   for(j in 1:afr.count){
#     ### checking model1 of phyloframe against all other models of benchamrk
#     the.int <- length(intersect(admix.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable))
#     gene.union <- union(admix.myfiles[i][[1]]$Variable,afr.myfiles[j][[1]]$Variable)
#     denom <- length(gene.union)
#     overlap <- (the.int/denom) * 100
#     temp[1,j+eur.count] <- overlap
#     # phyloframe
#     pf.the.int <- length(intersect(pf.admix.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable))
#     pf.gene.union <- union(pf.admix.myfiles[i][[1]]$Variable,pf.afr.myfiles[j][[1]]$Variable)
#     pf.denom <- length(pf.gene.union)
#     pf.overlap <- (pf.the.int/pf.denom) * 100
#     pf.temp[1,j+ eur.count] <- pf.overlap
#   }
#   for(j in 1:eas.count){
#     ### checking model1 of phyloframe against all other models of benchamrk
#     the.int <- length(intersect(admix.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable))
#     gene.union <- union(admix.myfiles[i][[1]]$Variable,eas.myfiles[j][[1]]$Variable)
#     denom <- length(gene.union)
#     overlap <- (the.int/denom) * 100
#     temp[1,j+ (eur.count + afr.count)] <- overlap
#     ## phyloframe ##
#     pf.the.int <- length(intersect(pf.admix.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable))
#     pf.gene.union <- union(pf.admix.myfiles[i][[1]]$Variable,pf.eas.myfiles[j][[1]]$Variable)
#     pf.denom <- length(pf.gene.union)
#     pf.overlap <- (pf.the.int/pf.denom) * 100
#     pf.temp[1,j + (eur.count + afr.count) ] <- pf.overlap
#   }
#   for(j in 1:admix.count){
#     the.int <- length(intersect(admix.myfiles[i][[1]]$Variable,admix.myfiles[j][[1]]$Variable))
#     gene.union <- union(admix.myfiles[i][[1]]$Variable,admix.myfiles[j][[1]]$Variable)
#     denom <- length(gene.union)
#     overlap <- (the.int/denom) * 100
#     temp[1,j+ (eur.count + afr.count + eas.count)] <- overlap
#     ## phyloframe ##
#     pf.the.int <- length(intersect(pf.admix.myfiles[i][[1]]$Variable,pf.admix.myfiles[j][[1]]$Variable))
#     pf.gene.union <- union(pf.admix.myfiles[i][[1]]$Variable,pf.admix.myfiles[j][[1]]$Variable)
#     pf.denom <- length(pf.gene.union)
#     pf.overlap <- (pf.the.int/pf.denom) * 100
#     pf.temp[1,j + (eur.count + afr.count + eas.count) ] <- pf.overlap
#   }
#   admix.df <- rbind(admix.df, temp)
#   pf.admix.df <- rbind(pf.admix.df, pf.temp)
# }



#bind all ancestry dataframes and create row labels
benchmark <-rbind(eur.df, afr.df,eas.df)
rownames(benchmark) <- paste0("benchmark_", rownames(benchmark))
phyloFrame <- rbind(pf.eur.df, pf.afr.df, pf.eas.df)
rownames(phyloFrame) <- paste0("PF_", rownames(phyloFrame))

all <- rbind(phyloFrame, benchmark)
#all <- rownames_to_column(all, "model_num")
all.mat <- as.matrix(all)
df.row <- c(rep("eur",eur.count), rep("afr",afr.count), rep("eas",eas.count),rep("eur",eur.count), rep("afr",afr.count), rep("eas",eas.count))
df <- c(rep("eur",eur.count), rep("afr",afr.count), rep("eas",eas.count))

#column annotation
for.annotation <- data.frame(Ancestry = df)
rownames(for.annotation) <- c.name
for.annotation$Ancestry <- as.factor(for.annotation$Ancestry)
#row annotation
for.row.annotation <- data.frame(Ancestry = df.row)
rownames(for.row.annotation) <- rownames(all)
for.row.annotation$Ancestry <- as.factor(for.row.annotation$Ancestry)

## correlation matrix
res <- cor(all.mat)
round(res, 2)

anc.color <- c("#CA4136","#1F619E","#496849")
names(anc.color) <- c("eur","afr","eas")
the.color <- list(Ancestry = anc.color)
library(RColorBrewer)
display.brewer.all()
# display.brewer.pal(n = 9, name = 'YlOrBr')
# brewer.pal(n = 9, name = "Reds")
library(pheatmap)
heatmap <- pheatmap(
  all.mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_colors = the.color,
  annotation_col = for.annotation,
  annotation_row = for.row.annotation,
  main = "Model Signature Overlap",
  colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"))(25)

)

heatmap




############################### ADDED FOR FIGURE 3 H AND I ###############################
## added to get percentage of top 6 COSMIC in each model PF and BM ##
#genes = "COL2A1","MUC16","S100A7","SLC34A2","AR","FOXA1","ESR1","GATA3","AR","ERBB4"
cosmic <- readr::read_tsv(here::here("Cosmic_Genes.tsv"))

## - find most common genes - ##
signatures <- readr::read_tsv(here::here("results", "TCGA_Breast_Gnomad4_corrected","model_runs","final_model_signatures.txt"))
signatures <- signatures[signatures$train_data != "mixed",]
signatures <- signatures[signatures$train_data != "admix",]

pf_sigs <- signatures[signatures$version == "phyloFrame",]
bm_sigs <- signatures[signatures$version == "benchmark",]

pf_tt <- pf_sigs[pf_sigs$gene %in% cosmic$`Gene Symbol`,]
dat_pf <- table(pf_tt$gene)
sorted_pf <- dat_pf[order(dat_pf, decreasing = TRUE)]
head(sorted_pf)
sorted_pf[1:30]

bm_tt <- bm_sigs[bm_sigs$gene %in% cosmic$`Gene Symbol`,]
dat_bm <- table(bm_tt$gene)
sorted_bm <- dat_bm[order(dat_bm, decreasing = TRUE)]
head(sorted_bm)

## - phyloFrame: CNTNAP2  COL2A1   MUC16   PTPRT    ESR1   FOXA1
## - benchmark: ESR1 FOXA1 GATA3    AR ERBB4 IL6ST

dat <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(dat) <- c("gene", "phyloFrame","Benchmark")
# note there are 27 total breast cancer models #

perc.pf <- (nrow(pf_sigs[pf_sigs$gene == "ESR1",])/20)*100
perc.bm <- (nrow(bm_sigs[bm_sigs$gene == "ESR1",])/20)*100
dat[nrow(dat)+1,] <- c("ESR1",perc.pf, perc.bm)

#PTPRT
perc.pf <- (nrow(pf_sigs[pf_sigs$gene == "PTPRT",])/20)*100
perc.bm <- (nrow(bm_sigs[bm_sigs$gene == "PTPRT",])/20)*100
dat[nrow(dat)+1,] <- c("PTPRT",perc.pf, perc.bm)

#MUC16
perc.pf <- (nrow(pf_sigs[pf_sigs$gene == "MUC16",])/20)*100
perc.bm <- (nrow(bm_sigs[bm_sigs$gene == "MUC16",])/20)*100
dat[nrow(dat)+1,] <- c("MUC16",perc.pf, perc.bm)
#COL2A1
perc.pf <- (nrow(pf_sigs[pf_sigs$gene == "CNTNAP2",])/20)*100
perc.bm <- (nrow(bm_sigs[bm_sigs$gene == "CNTNAP2",])/20)*100
dat[nrow(dat)+1,] <- c("CNTNAP2",perc.pf, perc.bm)
#COL2A1
perc.pf <- (nrow(pf_sigs[pf_sigs$gene == "COL2A1",])/20)*100
perc.bm <- (nrow(bm_sigs[bm_sigs$gene == "COL2A1",])/20)*100
dat[nrow(dat)+1,] <- c("COL2A1",perc.pf, perc.bm)

# CNTNAP2    MUC16  COL2A1     PTPRT     ESR1
# [1,] 96.296296 92.59259 96.2963 85.185185 81.48148
# [2,]  3.703704 14.81481  0.0000  3.703704 85.18519
dat <- tibble::column_to_rownames(dat, "gene")
dat <- as.matrix(dat)
dat2 <- apply(dat, 1, as.numeric)
rownames(dat2) <- c("PhyloFrame", "Benchmark")
library(pheatmap)
heatmap <- pheatmap(
  dat2,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Model Signature Overlap",
  colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"))(25)

  #colorRampPalette(c("#FCFBFD", "#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D"))(25)

  #colorRampPalette(c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704"))(25)
)
## logging the percentages
sink(here::here("results","files","cosmic_model_sig_sig_heatmap.log"))
cat("################## PHYLOFRAME MOST COMMONLY FOUND COSMIC GENES ##################")
dat2
sink()
heatmap

## get specific models
pf_genes <- rownames(dat)
pf_models <- pf_sigs[pf_sigs$gene %in% pf_genes,]
paste0(pf_models$gene ,pf_models$train_data, "_",pf_models$model_num)
pf_models[pf_models$gene == "CNTNAP2",]
pf_models[pf_models$gene == "COL2A1",]
pf_models[pf_models$gene == "MUC16",]
pf_models[pf_models$gene == "PTPRT",]
pf_models[pf_models$gene == "ESR1",]
#do same for benchmark genes
dat <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(dat) <- c("gene", "phyloFrame","Benchmark")
#
#ERBB4
perc.pf <- (nrow(pf_sigs[pf_sigs$gene == "ERBB4",])/20)*100
perc.bm <- (nrow(bm_sigs[bm_sigs$gene == "ERBB4",])/20)*100
dat[nrow(dat)+1,] <- c("ERBB4",perc.pf, perc.bm)
#AR
perc.pf <- (nrow(pf_sigs[pf_sigs$gene == "AR",])/20)*100
perc.bm <- (nrow(bm_sigs[bm_sigs$gene == "AR",])/20)*100
dat[nrow(dat)+1,] <- c("AR",perc.pf, perc.bm)
#GATA3
perc.pf <- (nrow(pf_sigs[pf_sigs$gene == "GATA3",])/20)*100
perc.bm <- (nrow(bm_sigs[bm_sigs$gene == "GATA3",])/20)*100
dat[nrow(dat)+1,] <- c("GATA3",perc.pf, perc.bm)
#FOXA1
perc.pf <- (nrow(pf_sigs[pf_sigs$gene == "FOXA1",])/20)*100
perc.bm <- (nrow(bm_sigs[bm_sigs$gene == "FOXA1",])/20)*100
dat[nrow(dat)+1,] <- c("FOXA1",perc.pf, perc.bm)
#ESR1
perc.pf <- (nrow(pf_sigs[pf_sigs$gene == "ESR1",])/20)*100
perc.bm <- (nrow(bm_sigs[bm_sigs$gene == "ESR1",])/20)*100
dat[nrow(dat)+1,] <- c("ESR1",perc.pf, perc.bm)

head(dat)
dat <- tibble::column_to_rownames(dat, "gene")
dat <- as.matrix(dat)
dat2 <- apply(dat, 1, as.numeric)
rownames(dat2) <- c("PhyloFrame","Benchmark")
library(pheatmap)
breaksList = seq(0, 100, by = 4)
heatmap <- pheatmap(
  dat2,
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Model Signature Overlap",
  colorRampPalette(c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"))(length(breaksList)),
  breaks = breaksList

  #colorRampPalette(c("#FCFBFD", "#EFEDF5", "#DADAEB", "#BCBDDC", "#9E9AC8", "#807DBA", "#6A51A3", "#54278F", "#3F007D"))(25)

  #colorRampPalette(c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801", "#A63603", "#7F2704"))(25)
)

pheatmap(expressionData[1:10, ], # Plots the first 10 genes of the dataset
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList) # Sets the breaks of the color scale as in breaksList


dat
sink(here::here("results","files","cosmic_model_sig_sig_heatmap.log"),append = TRUE)
cat("################## BENCHMARK MOST COMMONLY FOUND COSMIC GENES ##################")
dat2
sink()

bm_sigs[bm_sigs$gene == "ESR1",]
bm_sigs[bm_sigs$gene == "FOXA1",]
bm_sigs[bm_sigs$gene == "GATA3",]
bm_sigs[bm_sigs$gene == "AR",]
bm_sigs[bm_sigs$gene == "ERBB4",]


