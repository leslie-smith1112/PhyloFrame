eaf <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/PhyloFrame/data-raw/mean_enhancedAF_exome.tsv", col_names = TRUE)
# see how many genes have symbol names #
eaf_cut <-  eaf[!(is.na(eaf$gene)),]
eaf_cut <-  eaf_cut[!(eaf_cut$gene == "NA"),]
message("EAF_SYMBOL_DIM:", dim(eaf_cut))
message("EAF_SYMBOL_GENE_COUNT", length(unique(eaf_cut$gene)))
#to look at the dimensions of how many snps we lose.
eaf_non <- eaf[(is.na(eaf$gene)),]
message("EAF_NON_DIM:", dim(eaf_non))
message("EAF_NON_GENE_COUNT", length(unique(eaf_non$gene)))

#eaf_cut
library(magrittr)
library(tidyr)
library(ggplot2)

cutGenome <- eaf_cut %>% dplyr::select(-chrom, -position, -rs_id, -ref_allele, -alt_allele)
data.long <- gather(cutGenome,ancestry, EAF, nfe_seu:oth, factor_key = TRUE)

data.long$ancestry <- factor(data.long$ancestry, levels = c("afr", "amr", "asj","eas", "fin", "nfe", "nfe_est","nfe_nwe","nfe_onf","nfe_seu", "oth"))
e <- ggplot(data = data.long, aes(x=EAF,  fill = ancestry, color= ancestry)) + geom_density(alpha = 0.7) + theme_minimal() + scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 100000))+
  scale_fill_manual(values=c("#1F619E","#F4C867","#B2ABD2","#496849","#CE1256","#CA4136", "#A50F15", "#EF3B2C","#FB6A4A","#FC9272" ,"#BF812D"))  +
  scale_color_manual(values=c("#1F619E","#F4C867","#B2ABD2","#496849","#CE1256","#CA4136", "#A50F15", "#EF3B2C","#FB6A4A","#FC9272" ,"#BF812D")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")

e

ggsave(file = "/blue/kgraim/leslie.smith1/Repositories/PhyloFrame/data-raw/gnomad4.1_annotated/Gnomad2_EAF1.png", plot = e, width = 10, height = 8)

data.long$ancestry <- factor(data.long$ancestry, levels = c("afr", "amr", "asj","eas", "fin", "nfe", "nfe_est","nfe_nwe","nfe_onf","nfe_seu", "oth"))
e <- ggplot(data = data.long, aes(x=EAF,  fill = ancestry, color= ancestry)) + geom_density(alpha = 0.7) + theme_minimal() + scale_x_continuous(limits = c(0.5, 1)) + scale_y_continuous(limits = c(0, 100000))+
  scale_fill_manual(values=c("#1F619E","#F4C867","#B2ABD2","#496849","#CE1256","#CA4136", "#A50F15", "#EF3B2C","#FB6A4A","#FC9272" ,"#BF812D"))  +
  scale_color_manual(values=c("#1F619E","#F4C867","#B2ABD2","#496849","#CE1256","#CA4136", "#A50F15", "#EF3B2C","#FB6A4A","#FC9272" ,"#BF812D")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")

e

ggsave(file = "/blue/kgraim/leslie.smith1/Repositories/PhyloFrame/data-raw/gnomad4.1_annotated/Gnomad2_EAF.51.png", plot = e, width = 10, height = 8)

