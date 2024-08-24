## read in the disease network
## get the enhanced alleles for each ancestry
## barplot: will need the ancestry #EAFAlleles


# eaf <- readr::read_tsv(here::here("data-raw","gnomad4.1_annotated","gnomad.exomes.all.tsv"), col_names = TRUE)
# # see how many genes have symbol names #
# eaf_cut <-  eaf[!(is.na(eaf$SYMBOL)),]
# eaf_cut <-  eaf_cut[!(eaf_cut$SYMBOL == "NA"),]
# message("EAF_SYMBOL_DIM:", dim(eaf_cut))
# message("EAF_SYMBOL_GENE_COUNT", length(unique(eaf_cut$SYMBOL)))
# write.table(eaf_cut, here::here("data-raw", "gnomad4.1_annotated", "gnomad.exomes.symbol.tsv"), col.names = TRUE, row.names = FALSE)
# #to look at the dimensions of how many snps we lose.
# eaf_non <- eaf[(is.na(eaf$SYMBOL)),]
# write.table(eaf_non, here::here("data-raw", "gnomad4.1_annotated", "gnomad.exomes.no.symbol.tsv"), col.names = TRUE, row.names = FALSE)
# message("EAF_NON_DIM:", dim(eaf_non))
# message("EAF_NON_GENE_COUNT", length(unique(eaf_non$SYMBOL)))
# print(head(unique(eaf_non$GENE)))

############### added for GUI execution ###################
exome <- readr::read_tsv(here::here("data-raw","gnomad4.1_annotated","gnomad.exomes.all.tsv"), col_names = TRUE)
# exome <- exome[exome$SYMBOL != "NA",]
# exome <- exome[na.omit(exome$SYMBOL),]

#eaf_cut
library(magrittr)
library(tidyr)
library(ggplot2)
library(ggridges)
cutGenome <- exome |> dplyr::select(-CHROM, -POS, -RSID, -REF, -ALT, -GENE, -SYMBOL, -IMPACT, -CONSEQUENCE, -ALLELE, -DISTANCE)
data.long <- gather(cutGenome,ancestry, EAF, AFR_EAF:SAS_EAF, factor_key = TRUE)


head(data.long) #.00025
data.long$ancestry <- factor(data.long$ancestry, levels = c("AFR_EAF"	,"AMR_EAF",	"ASJ_EAF",	"EAS_EAF",	"FIN_EAF",	"MID_EAF",	"NFE_EAF",	"SAS_EAF"))
message(range(data.long$EAF))
message("Dimensions: ", dim(data.long))

# e <- ggplot(data = data.long, aes(x=EAF, y = ancestry, fill = ancestry)) +
#   geom_density_ridges() +
#   theme_ridges() +
#   theme(legend.position = "none") +
#   #scale_x_continuous(limits = c(-1, 1)) + #scale_y_continuous(limits = c(0, 100000))+
#   scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))
#   #scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#cb8de6","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
# e
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_none.png"), plot = e, width = 15, height = 8)
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_none.svg"), plot = e, width = 15, height = 8)
#
#
# e <- ggplot(data = data.long, aes(x=EAF, y = ancestry, fill = ancestry)) +
#   geom_density_ridges() +
#   theme_ridges() +
#   theme(legend.position = "none") +
#   scale_x_continuous(limits = c(0, 1)) + #scale_y_continuous(limits = c(0, 100000))+
#   scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))
# #scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#cb8de6","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
# e
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_0_1.png"), plot = e, width = 15, height = 8)
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_0_1.svg"), plot = e, width = 15, height = 8)

e <- ggplot(data = data.long, aes(x=EAF, y = ancestry, fill = ancestry)) +
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0.001, 1)) + #scale_y_continuous(limits = c(0, 100000))+
  scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))
#scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#cb8de6","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
e
ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_p001.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_p001.svg"), plot = e, width = 15, height = 8)

#
# e <- ggplot(data = data.long, aes(x=EAF, y = ancestry, fill = ancestry)) +
#   geom_density_ridges() +
#   theme_ridges() +
#   theme(legend.position = "none") +
#   scale_x_continuous(limits = c(0.2, 1)) + #scale_y_continuous(limits = c(0, 5000))+
#   scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))
#   #scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#cb8de5","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
# e
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_zoomed_2.png"), plot = e, width = 15, height = 8)
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_zoomed_2.svg"), plot = e, width = 15, height = 8)

# e <- ggplot(data = data.long, aes(x=EAF,  fill = ancestry, color= ancestry)) + geom_density(alpha = 0.6) + theme_bw() + #xlim(0,1) + ylim(0,10000) + #scale_x_continuous(limits = c(0.1, 1)) + #scale_y_continuous(limits = c(0, 50000))+
#   scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))  +
#   scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
# e
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_annotated_EAF.1.png"), plot = e, width = 15, height = 8)
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_annotated_EAF.1.svg"), plot = e, width = 15, height = 8)
#
#
# ggplot(data = afr_only, aes(x=EAF,  fill = ancestry, color= ancestry)) + geom_density(alpha = 0.6) + theme_bw() + #xlim(0,1) + ylim(0,10000) + #scale_x_continuous(limits = c(0.1, 1)) + #scale_y_continuous(limits = c(0, 50000))+
#   scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))  +
#   scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))





# # ## - COSMIC
cosmic_genes <- readr::read_tsv(here::here("Cosmic_Genes.tsv"))
cosmic_eaf <- exome[exome$SYMBOL %in% cosmic_genes$`Gene Symbol`,]
#cutCosmic <- cosmic_eaf %>% dplyr::select(AFR_EAF, EAS_EAF, NFE_EAF)
#data.long <- gather(cutCosmic,ancestry, EAF, AFR_EAF:NFE_EAF, factor_key = TRUE)

#data.long$ancestry <- factor(data.long$ancestry, levels = c("AFR_EAF", "EAS_EAF",	"NFE_EAF"))
cutGenome <- cosmic_eaf |> dplyr::select(-CHROM, -POS, -RSID, -REF, -ALT, -GENE, -SYMBOL, -IMPACT, -CONSEQUENCE, -ALLELE, -DISTANCE)

data.long <- gather(cutGenome,ancestry, EAF, AFR_EAF:SAS_EAF, factor_key = TRUE)

data.long$ancestry <- factor(data.long$ancestry, levels = c("AFR_EAF"	,"AMR_EAF",	"ASJ_EAF",	"EAS_EAF",	"FIN_EAF",	"MID_EAF",	"NFE_EAF",	"SAS_EAF"))

e <- ggplot(data = data.long, aes(x=EAF, y = ancestry, fill = ancestry)) +
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0.001, 1)) + #scale_y_continuous(limits = c(0, 100000))+
  scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))
#scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#cb8de6","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
e
ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_COSMIC.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_COSMIC.svg"), plot = e, width = 15, height = 8)


e <- ggplot(data = data.long, aes(x=EAF, y = ancestry, fill = ancestry)) +
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0.1, 1)) + #scale_y_continuous(limits = c(0, 100000))+
  scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))
#scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#cb8de6","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
e
ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_COSMIC1.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_COSMIC1.svg"), plot = e, width = 15, height = 8)


#message(range(data.long$EAF))
#cosmic_data.long <- pivot_longer(cutCosmic, c("AFR_EAF"	,	"EAS_EAF",	"NFE_EAF"),names_to='ancestry', values_to='EAF')
#cosmic_data.long$ancestry <- factor(data.long$ancestry, levels = c("AFR_EAF"	,	"EAS_EAF",	"NFE_EAF"))

# e <- ggplot(data = data.long, aes(x=EAF, y = ancestry,  fill = ancestry)) +
#   geom_density_ridges() +
#   theme_ridges() +
#   theme(legend.position = "none") +
#   #scale_x_continuous(limits = c(-1, 1)) + #scale_y_continuous(limits = c(0, 100000)) +
#   scale_fill_manual(values=c("#1F619E","#496849","#CA4136")) # +
#   #scale_color_manual(values=c("#1F619E","#496849","#CA4136")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
# e
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_afr_eas_nfe_none.png"), plot = e, width = 15, height = 8)
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_afr_eas_nfe_none.svg"), plot = e, width = 15, height = 8)
#
#
#
#
# e <- ggplot(data = data.long, aes(x=EAF,  y = ancestry, fill = ancestry)) +
#   geom_density_ridges() +
#   theme_ridges() +
#   theme(legend.position = "none") +
#   scale_x_continuous(limits = c(0, 1)) + #scale_y_continuous(limits = c(0, 5000))+
#   scale_fill_manual(values=c("#1F619E","#496849","#CA4136")) # +
#   #scale_color_manual(values=c("#1F619E","#496849","#CA4136")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
# e
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_zoomed_afr_eas_nfe_0_1.png"), plot = e, width = 15, height = 8)
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_zoomed_afr_eas_nfe_0_1.svg"), plot = e, width = 15, height = 8)
#
#
# e <- ggplot(data = data.long, aes(x=EAF, y = ancestry,  fill = ancestry)) +
#   geom_density_ridges() +
#   theme_ridges() +
#   theme(legend.position = "none") +
#   scale_x_continuous(limits = c(0.1, 1)) + #scale_y_continuous(limits = c(0, 100000)) +
#   scale_fill_manual(values=c("#1F619E","#496849","#CA4136")) # +
# #scale_color_manual(values=c("#1F619E","#496849","#CA4136")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
# e
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_afr_eas_nfe_p1.png"), plot = e, width = 15, height = 8)
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_afr_eas_nfe_p1.svg"), plot = e, width = 15, height = 8)
#
# e <- ggplot(data = data.long, aes(x=EAF, y = ancestry,  fill = ancestry)) +
#   geom_density_ridges() +
#   theme_ridges() +
#   theme(legend.position = "none") +
#   scale_x_continuous(limits = c(0.2, 1)) + #scale_y_continuous(limits = c(0, 100000)) +
#   scale_fill_manual(values=c("#1F619E","#496849","#CA4136")) # +
# #scale_color_manual(values=c("#1F619E","#496849","#CA4136")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
# e
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_afr_eas_nfe_zoomed_2.png"), plot = e, width = 15, height = 8)
# ggsave(file = here::here("figures","EAF_Density","gnomad4.1_eaf_ridge_afr_eas_nfe_zoomed_2.svg"), plot = e, width = 15, height = 8)
#
#
#











# ## load PhyloFrame package
# devtools::load_all()
#
# ## --- functions used the get enhanced alleles for each disease network --- ##
#
# ## function to be called on each ancestry to find enhanced alleles
# find_ancestry_enhanced <- function(eafs, curr_ancestry, alt_anc_1,alt_anc_2, alt_anc_3,
#                                    alt_anc_4, alt_anc_5, alt_anc_6, alt_anc_7){
#   ## for the current ancestry only keep alleles with frequency in defined EAF range
#   valid_af <- eafs[(eafs[,curr_ancestry] >= 0.2) & (eafs[,curr_ancestry] <= 1),]
#   ## only keep alleles enhanced in current ancestry
#   anc_snps <- valid_af[((valid_af[,curr_ancestry] > 0) & (valid_af[,alt_anc_1] < 0) & (valid_af[,alt_anc_2] < 0) & (valid_af[,alt_anc_3] < 0) & (valid_af[,alt_anc_4] < 0) &
#                           (valid_af[,alt_anc_5] < 0) & (valid_af[,alt_anc_6] < 0) & (valid_af[,alt_anc_7] < 0)) ,]
#   ## return unique gene list
#   return(na.omit(unique(anc_snps$SYMBOL)))
# }
#
# ## fucntion to be called on each disease network to get the enhanced alleles for that network
# find_enhanced <- function(network, eaf_exome){
#
#   ## only keep alleles on genes within the network
#   cut_eaf <- eaf_exome[eaf_exome$SYMBOL %in% network$Gene1 | eaf_exome %in% network$Gene2,]
#   afr <- find_ancestry_enhanced(cut_eaf, "AFR_EAF", "AMR_EAF","ASJ_EAF", "EAS_EAF","FIN_EAF", "NFE_EAF", "MID_EAF", "SAS_EAF")
#   amr <- find_ancestry_enhanced(cut_eaf, "AMR_EAF", "AFR_EAF","ASJ_EAF", "EAS_EAF", "FIN_EAF", "NFE_EAF", "MID_EAF", "SAS_EAF")
#   asj <- find_ancestry_enhanced(cut_eaf, "ASJ_EAF", "AFR_EAF","AMR_EAF", "EAS_EAF", "FIN_EAF", "NFE_EAF", "MID_EAF", "SAS_EAF")
#   eas <- find_ancestry_enhanced(cut_eaf, "EAS_EAF", "AFR_EAF","AMR_EAF", "ASJ_EAF", "FIN_EAF", "NFE_EAF", "MID_EAF", "SAS_EAF")
#   fin <- find_ancestry_enhanced(cut_eaf, "FIN_EAF", "AFR_EAF","AMR_EAF", "ASJ_EAF", "EAS_EAF", "NFE_EAF", "MID_EAF", "SAS_EAF")
#   nfe <- find_ancestry_enhanced(cut_eaf, "NFE_EAF", "AFR_EAF","AMR_EAF", "ASJ_EAF", "EAS_EAF", "FIN_EAF", "MID_EAF", "SAS_EAF")
#   mid <- find_ancestry_enhanced(cut_eaf, "MID_EAF", "AFR_EAF","AMR_EAF", "ASJ_EAF", "EAS_EAF", "FIN_EAF","NFE_EAF", "SAS_EAF")
#   sas <- find_ancestry_enhanced(cut_eaf, "SAS_EAF", "AFR_EAF","AMR_EAF", "ASJ_EAF", "EAS_EAF", "FIN_EAF","NFE_EAF", "MID_EAF")
#   ancestry_genes <- data.frame("ancestry" = c("afr","amr","asj","eas","fin","nfe","mid","sas"), "eaf_count" = c(length(afr),length(amr),length(asj), length(eas),length(fin),length(nfe),length(mid),length(sas)))
#   return(ancestry_genes)
# }
#
# ## --- end of function definitions --- ##
# exome <- readr::read_tsv(here::here("data-raw","gnomad4.1_annotated","gnomad.exomes.all.tsv"), col_names = TRUE)
#
# ## load disease networks
# breast_net <- load_breast_network()
# thyroid_net <- load_thyroid_network()
# uterine_net <- load_uterine_network()
#
#
# # exome <- exome[exome$SYMBOL != "NA",]
# # exome <- exome[na.omit(exome$SYMBOL),]
# brca_counts <- find_enhanced(breast_net, exome)
# thca_counts <- find_enhanced(thyroid_net, exome)
# ucec_counts <- find_enhanced(uterine_net, exome)
#
# ## write to file
# write.table(brca_counts, here::here("results","files","brca_network_eaf_counts.tsv"),sep = "\t",col.names = TRUE, row.names = FALSE)
# write.table(thca_counts, here::here("results","files","thca_network_eaf_counts.tsv"),sep = "\t",col.names = TRUE, row.names = FALSE)
# write.table(ucec_counts, here::here("results","files","ucec_network_eaf_counts.tsv"),sep = "\t",col.names = TRUE, row.names = FALSE)
#
# ## read in from file (typically first hald of script done in command line with increased memory)
# brca_counts <- readr::read_tsv(here::here("results","files","brca_network_eaf_counts.tsv"))
# thca_counts <- readr::read_tsv(here::here("results","files","thca_network_eaf_counts.tsv"))
# ucec_counts <- readr::read_tsv(here::here("results","files","ucec_network_eaf_counts.tsv"))
#
# ## brca
# e <- ggplot(data = brca_counts, aes(x=ancestry, y=eaf_count, fill = ancestry, color= ancestry)) + geom_bar(stat = "identity") + theme_bw() + #+ scale_x_continuous(limits = c(0.2, 1)) + #scale_y_continuous(limits = c(0, 50000))+
#   scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))  +
#   scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
# e
# ggsave(file = here::here("figures","Supp_Disease_Network_Eafs","brca.png"), plot = e, width = 10, height = 8)
# ggsave(file = here::here("figures","Supp_Disease_Network_Eafs","brca.svg"), plot = e, width = 10, height = 8)
#
# ## thca
# e <- ggplot(data = thca_counts, aes(x=ancestry, y=eaf_count, fill = ancestry, color= ancestry)) + geom_bar(stat = "identity") + theme_bw() +#+ scale_x_continuous(limits = c(0.2, 1)) + #scale_y_continuous(limits = c(0, 50000))+
# scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))  +
#   scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
# e
# ggsave(file = here::here("figures","Supp_Disease_Network_Eafs","thca.png"), plot = e, width = 10, height = 8)
# ggsave(file = here::here("figures","Supp_Disease_Network_Eafs","thca.svg"), plot = e, width = 10, height = 8)
#
# ## ucec
# e <- ggplot(data = ucec_counts, aes(x=ancestry, y=eaf_count, fill = ancestry, color= ancestry)) + geom_bar(stat = "identity") + theme_bw() +#+ scale_x_continuous(limits = c(0.2, 1)) + #scale_y_continuous(limits = c(0, 50000))+
# scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))  +
#   scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
# e
# ggsave(file = here::here("figures","Supp_Disease_Network_Eafs","ucec.png"), plot = e, width = 10, height = 8)
# ggsave(file = here::here("figures","Supp_Disease_Network_Eafs","ucec.svg"), plot = e, width = 10, height = 8)
#
# # p<-ggplot(data=df, aes(x=dose,
# #   geom_bar(stat="identity")
#
# # trimmed_graph <- network[((network$Connection >= 0.2) & (network$Connection < 0.51)),]
