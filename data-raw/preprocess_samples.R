
#----------------- Add ancestry data to project -----------------#
ancestry <- readxl::read_xlsx(here::here("data-raw","tcga_estimated_ancestry.xlsx"),skip = 1)
ancestry  <- ancestry[ancestry$consensus_ancestry %in% c("eur", "afr","eas" ,"eas_admix","afr_admix","eur_admix","admix"),]
# paste on -01 to match epxression matrices
ancestry$patient <- paste0(ancestry$patient, "-01")

usethis::use_data(ancestry, overwrite = FALSE)

#----------------- Preprocess TCGA Breast -----------------#
expression <- readr::read_tsv(here::here("data-raw","brca_data_mrna_seq_v2_rsem.txt"), col_names = TRUE)
clinical <- readr::read_tsv(here::here("data-raw","brca_tcga_pan_can_atlas_2018_clinical_data.tsv"), col_names = TRUE)
ancestry_breast <- ancestry[ancestry$tumor_type == "BRCA",]
ancestry_breast <- ancestry_breast |> dplyr::select(patient, consensus_ancestry)

#make sure we dont have any NAs in the gene name
expression <- expression[!is.na(expression$Hugo_Symbol),] #prev in trim.expr.matrix
#get rid of entrez id
expression <- subset(expression, select =  -Entrez_Gene_Id) #prev in trim.expr.matrix
expression <- expression |> #prev in main script
  dplyr::group_by(Hugo_Symbol) |>
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)

# add subtype prediction to matrix
genes <- expression$Hugo_Symbol

#add.subtype.breast <- function(expression, tcga.clinical){
samples <- colnames(expression)
#exclude gene name
samples <- samples[-1]
#get rid of genes - later added as row names
expression <- expression[,-1]
#now we have genes as column names and samples on the rows
tran <- data.table::transpose(expression)
tran <- log2(tran + 1)
tran$sample_id <- samples

#cut clinical data
clin.cut <- data.frame(clinical$`Sample ID`,clinical$Subtype)
colnames(clin.cut) <- c( "sample_id","subtype")
temp <- merge(x = tran, y = clin.cut, by = "sample_id")
#expression_breast <- temp[temp$subtype %in% c("BRCA_Basal","BRCA_LumA","BRCA_LumB"),]
##BELOW IS THE TEMP ADD TO GET THE SAME ORDERING OF SAMPLES AS THE INITAL PHYLOFRAME RUN ##
tcga.binomial <- temp[temp$subtype == "BRCA_Basal",] #160 patients
tcga.binomial1 <- temp[temp$subtype == "BRCA_LumA",] #20 patients
binomial <- rbind(tcga.binomial,tcga.binomial1)
tcga.binomial2 <- temp[temp$subtype == "BRCA_LumB",]
binomial2 <- rbind(binomial,tcga.binomial2)
expression_breast <- na.omit(binomial2)

pattern1 <- "BRCA_LumA"
pattern2 <- "BRCA_LumB"
pattern3 <- "BRCA_Basal"

#### replace subtypes specifics for logistic regression model####
expression_breast$subtype <- stringi::stri_replace_all_fixed(expression_breast$subtype, pattern2,"Luminal")
expression_breast$subtype <- stringi::stri_replace_all_fixed(expression_breast$subtype, pattern1,"Luminal")
expression_breast$subtype <- stringi::stri_replace_all_fixed(expression_breast$subtype, pattern3,"Basal")
expression_breast$subtype <- as.factor(expression_breast$subtype)

colnames(expression_breast) <- c("sample_id", genes, "subtype")
rownames(expression_breast) <- expression_breast$sample_id
expression_breast <- expression_breast[,-1]
#expression_breast <- expression_breast[rownames(expression_breast) %in% ancestry_breast$patient,]
temp_breast <- expression_breast[rownames(expression_breast) %in% ancestry_breast$patient,]

ancestry_breast  <- ancestry_breast[ancestry_breast$patient %in% rownames(expression_breast),]
usethis::use_data(expression_breast, overwrite = FALSE)
usethis::use_data(ancestry_breast, overwrite = FALSE) ## ancestry for the samples we are using in our study

#----------------- Preprocess TCGA Thyroid -----------------#
expression <- readr::read_tsv(here::here("data-raw","thca_data_mrna_seq_v2_rsem.txt"), col_names = TRUE)
clinical <- readr::read_tsv(here::here("data-raw","thca_tcga_pan_can_atlas_2018_clinical_data.tsv"))
ancestry_thyroid <- ancestry[ancestry$tumor_type == "THCA",]
ancestry_thyroid <- ancestry_thyroid |> dplyr::select(patient, consensus_ancestry)
#narrow samples down to those actually in the expression matrix

#make sure we dont have any NAs in the gene name
expression <- expression[!is.na(expression$Hugo_Symbol),]
#get rid of entrez id
expression <- subset(expression, select =  -Entrez_Gene_Id)

expression <- expression |>
  dplyr::group_by(Hugo_Symbol) |>
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)

genes <- expression$Hugo_Symbol
samples <- colnames(expression)
samples <- samples[-1] #get rid of Hugo Symbol at the beginning of the list

expression <- expression[,-1] #get rid of genes - later added as col names
tran <- data.table::transpose(expression)
tran <- log2(tran + 1)
tran$sample_id <- samples

#cut clinical data
clin.cut <- data.frame( clinical$`Sample ID`, clinical$`American Joint Committee on Cancer Metastasis Stage Code`)
colnames(clin.cut) <- c( "sample_id","subtype")
#narrow down to subtypes we want
clin.cut <- clin.cut[clin.cut$subtype == "M0" | clin.cut$subtype == "MX",]

expression_thyroid <- merge(x = tran, y = clin.cut, by = "sample_id")
expression_thyroid$subtype <- as.factor(expression_thyroid$subtype)

colnames(expression_thyroid) <- c("sample_id", genes, "subtype")
rownames(expression_thyroid) <- expression_thyroid$sample_id
expression_thyroid <- expression_thyroid[,-1]
expression_thyroid <- expression_thyroid[rownames(expression_thyroid) %in% ancestry_thyroid$patient,]
ancestry_thyroid  <- ancestry_thyroid[ancestry_thyroid$patient %in% rownames(expression_thyroid),]
usethis::use_data(expression_thyroid, overwrite = FALSE)
usethis::use_data(ancestry_thyroid, overwrite = FALSE) ## ancestry only for those in our study


#----------------- Preprocess TCGA Uterine -----------------#
expression <- readr::read_tsv(here::here("data-raw","ucec_data_mrna_seq_v2_rsem.txt"), col_names = TRUE)
clinical <- readr::read_tsv(here::here("data-raw","ucec_tcga_pan_can_atlas_2018_clinical_data.tsv"), col_names = TRUE)
ancestry_uterine <- ancestry[ancestry$tumor_type == "UCEC",]
ancestry_uterine <- ancestry_uterine |> dplyr::select(patient, consensus_ancestry)

#make sure we dont have any NAs in the gene name
expression <- expression[!is.na(expression$Hugo_Symbol),]
#get rid of entrez id
expression <- subset(expression, select =  -Entrez_Gene_Id)

expression <- expression |>
  dplyr::group_by(Hugo_Symbol) |>
  dplyr::summarise_if(is.numeric, mean, na.rm = TRUE)

genes <- expression$Hugo_Symbol
samples <- colnames(expression)
samples <- samples[-1] #get rid of Hugo Symbol at the beginning of the list

expression <- expression[,-1] #get rid of genes - later added as col names
tran <- data.table::transpose(expression)
tran <- log2(tran + 1)
tran$sample_id <- samples

#cut clinical data
clin.cut <- data.frame(clinical$`Sample ID`, clinical$`Tumor Type`)
colnames(clin.cut) <- c( "sample_id","subtype")

#narrow down to subtypes we want
clin.cut <- clin.cut[clin.cut$subtype == "Endometrioid Endometrial Adenocarcinoma" | clin.cut$subtype == "Serous Endometrial Adenocarcinoma",]
pattern1 <- "Endometrioid Endometrial Adenocarcinoma"
pattern2 <- "Serous Endometrial Adenocarcinoma"

clin.cut$subtype <- stringi::stri_replace_all_fixed(clin.cut$subtype, pattern1,"Endometrioid")
clin.cut$subtype <- stringi::stri_replace_all_fixed(clin.cut$subtype, pattern2,"Serous")

temp <- merge(x = tran, y = clin.cut, by = "sample_id")
expression_uterine <- merge(x = tran, y = clin.cut, by = "sample_id")
expression_uterine$subtype <- as.factor(expression_uterine$subtype)

colnames(expression_uterine) <- c("sample_id", genes, "subtype")
rownames(expression_uterine) <- expression_uterine$sample_id
expression_uterine <- expression_uterine[,-1]
expression_uterine <- expression_uterine[rownames(expression_uterine) %in% ancestry_uterine$patient,]
ancestry_uterine  <- ancestry_uterine[ancestry_uterine$patient %in% rownames(expression_uterine),]
usethis::use_data(expression_uterine, overwrite = FALSE)
usethis::use_data(ancestry_uterine, overwrite = FALSE)



