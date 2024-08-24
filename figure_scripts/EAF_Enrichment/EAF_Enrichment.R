library(magrittr)
library(tidyr)
library(ggplot2)
library(ggridges)
exome <- readr::read_tsv(here::here("data-raw","gnomad4.1_annotated","gnomad.exomes.all.tsv"), col_names = TRUE)

## all gene plots
cutGenome <- exome |> dplyr::select(-CHROM, -POS, -RSID, -REF, -ALT, -GENE, -SYMBOL, -IMPACT, -CONSEQUENCE, -ALLELE, -DISTANCE)
data.long <- gather(cutGenome,ancestry, EAF, AFR_EAF:SAS_EAF, factor_key = TRUE)
data.long$ancestry <- factor(data.long$ancestry, levels = c("AFR_EAF"	,"AMR_EAF",	"ASJ_EAF",	"EAS_EAF",	"FIN_EAF",	"MID_EAF",	"NFE_EAF",	"SAS_EAF"))

#trying a bunch of thresholds here because hipergator is going down -.-

e <- ggplot(data = data.long, aes(x=EAF,  fill = ancestry, color= ancestry)) + geom_density(alpha = 0.6) + theme_bw() + coord_cartesian(xlim = c(-1, 1), ylim = c(0,25000))+
scale_x_continuous(breaks = seq(-1, 1, by = 0.2)) +   scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))  +
  scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b")) + geom_vline(xintercept=0.2)
e

ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_density_ylim25000.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_density_ylim25000.svg"), plot = e, width = 15, height = 8)

e <- ggplot(data = data.long, aes(x=EAF,  fill = ancestry, color= ancestry)) + geom_density(alpha = 0.6) + theme_bw() + coord_cartesian(xlim = c(-1, 1))+
scale_x_continuous(breaks = seq(-1, 1, by = 0.2)) +   scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))  +
  scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b")) + geom_vline(xintercept=0.2)
e

ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_density_hi.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_density_hi.svg"), plot = e, width = 15, height = 8)

e <- ggplot(data = data.long, aes(x=EAF,  fill = ancestry, color= ancestry)) + geom_density(alpha = 0.6) + theme_bw() + coord_cartesian(xlim = c(0.2, 0.6),ylim=c(0,10000))+
scale_x_continuous(breaks = seq(0.2, 0.6, by = 0.2)) +   scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))  +
  scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.2)
e

ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_density_zoom.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_density_zoom.svg"), plot = e, width = 15, height = 8)


e <- ggplot(data = data.long, aes(x=EAF,  fill = ancestry, color= ancestry)) + geom_density(alpha = 0.6) + theme_bw() + coord_cartesian(xlim = c(-1, 1), ylim = c(0,25000))+
scale_x_continuous(breaks = seq(-1, 1, by = 0.2)) +   scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))  +
  scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b")) + geom_vline(xintercept=0.2)
e
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_density_ylim25000.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_density_ylim25000.svg"), plot = e, width = 15, height = 8)


## make ridge plot as well
e <- ggplot(data = data.long, aes(x=EAF, y = ancestry, fill = ancestry)) +
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0.01, 1))+
 #coord_cartesian(xlim = c(0.01, 1)+
  geom_vline(xintercept=0.2) + #scale_y_continuous(limits = c(0, 100000))+
  scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))
#scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#cb8de6","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
e
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_ridgeBB.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_ridgeBB.svg"), plot = e, width = 15, height = 8)

## small scale change
e <- ggplot(data = data.long, aes(x=EAF, y = ancestry, fill = ancestry)) +
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0.1, 1))+
  #coord_cartesian(xlim = c(0.1, 1) +
  geom_vline(xintercept=0.2) + #scale_y_continuous(limits = c(0, 100000))+
  scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))
#scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#cb8de6","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
e
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_ridgep1BB.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_ridgep1BB.svg"), plot = e, width = 15, height = 8)


e <- ggplot(data = data.long, aes(x=EAF, y = ancestry, fill = ancestry)) +
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0.1, 1))+
  #coord_cartesian(xlim = c(0, 1) +
  geom_vline(xintercept=0.2) + #scale_y_continuous(limits = c(0, 100000))+
  scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))
#scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#cb8de6","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
e
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_ridgep0_nocartesian.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_ridgep0_nocartesian.svg"), plot = e, width = 15, height = 8)


e <- ggplot(data = data.long, aes(x=EAF, y = ancestry, fill = ancestry)) +
  geom_density_ridges() +
  theme_ridges() +
  theme(legend.position = "none") +
  scale_x_continuous(limits = c(0.2, 1))+
  #coord_cartesian(xlim = c(0.2, 1) +
  geom_vline(xintercept=0.2) + #scale_y_continuous(limits = c(0, 100000))+
  scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))
#scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#cb8de6","#CA4136",  "#86dc5b")) #+ geom_vline(xintercept=0.001, size=0.7, color="#7ffa04")
e
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_ridgep1_cartesianp2BB.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_eaf_ridgep1_cartesianp2BB.svg"), plot = e, width = 15, height = 8)


# Cosmic Genes
cosmic_genes <- readr::read_tsv(here::here("Cosmic_Genes.tsv"))
cosmic_eaf <- exome[exome$SYMBOL %in% cosmic_genes$`Gene Symbol`,]
cutGenome <- cosmic_eaf |> dplyr::select(-CHROM, -POS, -RSID, -REF, -ALT, -GENE, -SYMBOL, -IMPACT, -CONSEQUENCE, -ALLELE, -DISTANCE)

data.long <- gather(cutGenome,ancestry, EAF, AFR_EAF:SAS_EAF, factor_key = TRUE)
data.long$ancestry <- factor(data.long$ancestry, levels = c("AFR_EAF"	,"AMR_EAF",	"ASJ_EAF",	"EAS_EAF",	"FIN_EAF",	"MID_EAF",	"NFE_EAF",	"SAS_EAF"))

e <- ggplot(data = data.long, aes(x=EAF,  fill = ancestry, color= ancestry)) + geom_density(alpha = 0.6) + theme_bw() + coord_cartesian(xlim = c(-1, 1))+
scale_x_continuous(breaks = seq(-1, 1, by = 0.2)) +   scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))  +
  scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b")) + geom_vline(xintercept=0.2)
e

ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_cosmic_eaf_density.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_cosmic_eaf_density.svg"), plot = e, width = 15, height = 8)


## started from here: run at command line then remove next 3 lines

e <- ggplot(data = data.long, aes(x=EAF,  fill = ancestry, color= ancestry)) + geom_density(alpha = 0.6) + theme_bw() + coord_cartesian(xlim = c(-1, 1), ylim = c(0,100000))+
scale_x_continuous(breaks = seq(-1, 1, by = 0.2)) +   scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))  +
  scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b")) + geom_vline(xintercept=0.2)
e

ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_cosmic_eaf_density_cut.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_cosmic_eaf_density_cut.svg"), plot = e, width = 15, height = 8)

## if you have time:

e <- ggplot(data = data.long, aes(x=EAF,  fill = ancestry, color= ancestry)) + geom_density(alpha = 0.6) + theme_bw() + coord_cartesian(xlim = c(-1, 1), ylim = c(0,25000))+
scale_x_continuous(breaks = seq(-1, 1, by = 0.2)) +   scale_fill_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b"))  +
  scale_color_manual(values=c("#1F619E","#FFA500","#710765","#496849","#f61d6b","#be5c05","#CA4136",  "#86dc5b")) + geom_vline(xintercept=0.2)
e
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_cosmic_eaf_density_ylim25000.png"), plot = e, width = 15, height = 8)
ggsave(file = here::here("figures","EAF_Plots","gnomad4.1_coscmic_eaf_density_ylim25000.svg"), plot = e, width = 15, height = 8)




