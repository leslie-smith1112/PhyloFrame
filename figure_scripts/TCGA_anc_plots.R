## TCGA PLOT
head(ancestry)
dat <- readxl::read_xlsx(here::here("data-raw","tcga_estimated_ancestry.xlsx"),skip = 1)
table(dat$consensus_ancestry)

## create map for ancestries
anc_map <- data.frame("KEY" = c("eas_admix","admix","sas_admix","sas","amr","eur_admix","afr_admix","afr","eas","eur"),
                      "VALUE" = c("Asian","N/A","Asian", "Asian","American","European","African","African","Asian","European"))

dat$assigned <- plyr::mapvalues(dat$consensus_ancestry, anc_map$KEY, anc_map$VALUE)
table(dat$assigned)

diseases <- unique(dat$tumor_type)
ancestries <- unique(dat$assigned)
ancestries <- ancestries[-5] ## - gettingrid of NA
library(ggplot2)
library(svglite)
lapply(diseases, plot_ancestry)
lapply(ancestries, plot_disease)

## function creates 33 plots ( 1 for each of the diseases ) with sample counts for each ancestry
plot_ancestry <- function(disease){
  graph_colors <- data.frame("African" = "#1F619E", "American" = "#f7b018", "Asian" = "#85a93c" , "European" = "#CA4136")
  disease_dat <- dat[dat$tumor_type == disease,]
  anc_count <- as.data.frame(table(disease_dat$assigned))
  colnames(anc_count) <- c("Ancestry","Sample Count")
  p <- ggplot(data=anc_count, aes(x=Ancestry, y=`Sample Count`, color = Ancestry, fill = Ancestry)) + geom_bar(stat="identity") + theme_minimal() + coord_cartesian(ylim=c(0,850)) +
    scale_color_manual(values=c("African" = "#1F619E", American = "#f7b018", Asian = "#85a93c" , European = "#CA4136")) +
    scale_fill_manual(values=c("African" = "#1F619E", American = "#f7b018", Asian = "#85a93c" , European = "#CA4136")) +
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) + labs(title=paste0(disease," Ancestry Sample Count"))

  #ggplot(dat.brca, aes(color=Training_Data,x=BM_Recall, y=PF_Recall)) + geom_point() + geom_abline(intercept=0) + lims(y=c(0,850) )+ scale_colour_manual(values=cols);
  #ggsave('BRCA_recall.png',width=5.1,height=4)
  p
  ggsave(here::here("figures","TCGA_figure",paste0(disease,"_tcga_ancestry.png")),plot = last_plot(),width = 7, height = 5, bg = "white")
  ggsave(here::here("figures","TCGA_figure",paste0(disease,"_tcga_ancestry.svg")),plot = last_plot(),width = 7, height = 5, bg = "white")

}
## function creates 4 plots ( 1 for each of the assigned ancestries ) with sample counts for each disease
plot_disease <- function(ancestry){
  graph_colors <- data.frame("African" = "#1F619E", "American" = "#f7b018", "Asian" = "#85a93c" , "European" = "#CA4136")
  anc_color = graph_colors[,ancestry]
  ancestry_dat <- dat[dat$assigned == ancestry,]
  disease_count <- as.data.frame(table(ancestry_dat$tumor_type))
  colnames(disease_count) <- c("Disease","Sample Count")
  p <- ggplot(data=disease_count, aes(x=Disease, y=`Sample Count`)) + geom_bar(stat="identity",color = anc_color, fill = anc_color) + theme_minimal() + coord_cartesian(ylim=c(0,850)) +
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=14,face="bold")) + labs(title=paste0(ancestry," Disease Sample Count"))
    #scale_color_manual(values=c("African" = "#1F619E", American = "#f7b018", Asian = "#85a93c" , European = "#CA4136")) +
    #scale_fill_manual(values=c("African" = "#1F619E", American = "#f7b018", Asian = "#85a93c" , European = "#CA4136")) +

  p
  ggsave(here::here("figures","TCGA_figure",paste0(ancestry,"_tcga_disease.png")),plot = last_plot(),width = 9, height = 5, bg = "white")
  ggsave(here::here("figures","TCGA_figure",paste0(ancestry,"_tcga_disease.svg")),plot = last_plot(),width = 9, height = 5, bg = "white")

}

# library
library(ggplot2)
library(svglite)

# create a datase

diseases <- unique(dat$tumor_type)
ancestries <- unique(dat$assigned)
#dat
master <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(master) <- c("Disease","Ancestry", "Sample_count")
for(i in 1:length(diseases)){
  curr_disease <- diseases[i]
  curr_dat <- dat[dat$tumor_type == curr_disease,]
  count <- as.data.frame(table(curr_dat$assigned))
  disease_row <- data.frame("Disease" = rep(curr_disease,nrow(count)), "Ancestry" = count$Var1,
             "Sample_count" = count$Freq)
  master <- rbind(master, disease_row)
}
dim(master)
length(unique(master$Disease))

p <- ggplot(master, aes(fill=Ancestry, y=Sample_count, x=Disease)) +
  geom_bar(position="stack", stat="identity") +
  scale_color_manual(values=c(African = "#1F619E", American = "#FFA500", Asian = "#496849" , European = "#CA4136")) +
  scale_fill_manual(values=c(African = "#1F619E", American = "#FFA500", Asian = "#496849" , European = "#CA4136")) +
  theme(axis.text=element_text(size=8, face = "bold"),
        axis.title=element_text(size=14,face="bold"),
        panel.grid.major = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank())

p
ggsave("test.svg", plot=p, width=10, height=8, device = "svg")


library(grid)
library(Cairo)






