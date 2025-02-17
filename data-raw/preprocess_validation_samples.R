# preprocessing for validation studies
#----------------- Preprocess Multi Validation Study -----------------#

##------------- GET METADATA FOR STUDIES -------------##
## - Script used to harmonize the metadata from studies used in validation set - ##
## initially harmonizes newly added studies - then combines with previously used M. Davis subsaharan study

#old file: metadata <- readr::read_csv(here::here("data-raw","validation","Breast_Validation_Metadata.csv"))
#file from the map created with collaborators for mapping of subtypes and ancestry:
metadata <- readxl::read_xlsx(here::here("data-raw","validation","Metadata_to_upload.xlsx"))
colnames(metadata)
metadata[1:5,1:10]

## - studies to exclude:GSE113184
exclude <- c("GSE223181","GSE113184","GSE48390","GSE15852", "GSE142102","PROTEOGENOMIC_CPTACCELL2020")
metadata <- metadata[!(metadata$STUDY_ID %in% exclude),]
dim(metadata)

#concatinate the subtype columns for mapping
metadata$concat_sub <- paste(metadata$SUBTYPE_CONSENSUS,metadata$IHC_SUBTYPE,metadata$ER_STATUS, metadata$PR,
                             metadata$HER2,metadata$PAM50_SUBTYPE,sep = "_")

## - read in subtype and ancestry map file, made in bash with Makefile - ##
mappings <- readxl::read_xlsx(here::here("data-raw","validation","subtype_map.xlsx"))
dim(mappings)
mappings[1:5,1:5]

## - map the subtype information and the ancestry information
metadata$subtype_assigned_by_us <- plyr::mapvalues(metadata$concat_sub, mappings$SUBTYPE_KEY, mappings$SUBTYPE_VALUE)
metadata$ancestry_assigned_by_us <- plyr::mapvalues(metadata$`RACE/ETHNICITY`, mappings$`UNIQUE VALUES TOGETHER`, mappings$Continental)

## - adding the other ancestry info for fun
metadata$ancestry_assigned_by_us_tcga <- plyr::mapvalues(metadata$`RACE/ETHNICITY`, mappings$`UNIQUE VALUES TOGETHER`, mappings$TCGA_Level)
metadata$ancestry_assigned_by_us_admix_info <- plyr::mapvalues(metadata$`RACE/ETHNICITY`, mappings$`UNIQUE VALUES TOGETHER`, mappings$`Likely Admixed`)

metadata[1:5,c("concat_sub","subtype_assigned_by_us", "RACE/ETHNICITY","ancestry_assigned_by_us")]
## - before cutting invalid subtypes:  1536 samples, after: 1030
valid_subtypes <- metadata[metadata$subtype_assigned_by_us %in% c("Luminal","Basal"),]
valid_subtypes <- valid_subtypes[valid_subtypes$`TUMOR/NORMAL` == "tumor",]
dim(valid_subtypes)

table(valid_subtypes$ancestry_assigned_by_us, valid_subtypes$subtype_assigned_by_us)
#Continental:
#           Basal Luminal
# African    105      46
# American    50     166
# Asian       59     225
# European   112     250
# None        14       3

table(valid_subtypes$ancestry_assigned_by_us_tcga, valid_subtypes$subtype_assigned_by_us)
#TCGA:
#       Basal Luminal
# AFR    105      46
# AMR     50     166
# EAS     59     225
# EUR    112     250
# None    14       3

## GSE225846 has S_ appended in the metadata but not the expression data for some reason
valid_subtypes[valid_subtypes$STUDY_ID == "GSE225846",]$SAMPLE_ID <- sub('..','',valid_subtypes[valid_subtypes$STUDY_ID == "GSE225846",]$SAMPLE_ID)

#write.table(valid_subtypes, here::here("data-raw","validation","Breast_Validation_Metadata_annotated.tsv"), sep = "\t",col.names = TRUE, row.names = FALSE)

#TODO REMOVE
## - meta for the GEO studies:
# validation_raw_expr[1:5,1:5]
# validation_met_trimmed[1:5,1:5]
# dim(validation_met_trimmed)

## - add subsaharan metadata
meta <- readr::read_tsv(here::here("data-raw","validation","GSE211167","TNBC_African_validation_set.tsv"))
colnames(meta)
meta$Label <- substr(meta$ID, 17, 24)
head(meta$Label)
## meta - here we are adding the subsaharan data set to the new validatio set
meta$race_cut <- plyr::mapvalues(meta$race, c("self-reported race: African American", "self-reported race: African/Ethiopian",
                                              "self-reported race: African/Ghanaian"), c("African American", "African Ethiopian",
                                                                                         "African Ghanaian" ))
meta$ancestry_assigned_by_us <- "African"
meta$subtype_assigned_by_us <- "Basal"
# creating dataframe to for sub saharan dataset to match the information in dataframe with other validation sets
sub_add <- data.frame("SAMPLE_ID" = meta$Label, "PATIENT_ID" = meta$geo_accession, "TUMOR/NORMAL" = rep("tumor", nrow(meta)), "STUDY_ID " = meta$series,
                      "RACE" = meta$race_cut, "ETHNICITY" = rep(NA, nrow(meta)), `RACE/ETHNICITY` = meta$race_cut, "AGE" = rep(NA, nrow(meta)),
                      "SUBTYPE_CONSENSUS" = rep("Basal", nrow(meta)), "IHC_SUBTYPE" = rep(NA, nrow(meta)), "ER_STATUS" = rep(NA, nrow(meta)),
                      "PR" = rep(NA, nrow(meta)),"HER2" = rep(NA, nrow(meta)),"PAM5_SUBTYPE" = rep("Basal", nrow(meta)),
                      "CANCER_TYPE_DETAILED" = rep(NA, nrow(meta)),"concat_sub" = rep(NA, nrow(meta)), "subtype_assigned_by_us" = rep("Basal", nrow(meta)),
                      "ancestry_assigned_by_us" = rep("African", nrow(meta)),"ancestry_assigned_by_us_tcga" = rep("AFR", nrow(meta)))
sub_add$ancestry_assigned_by_us_admix_info <- plyr::mapvalues(sub_add$RACE.ETHNICITY,c("African American", "African Ethiopian",
                                                                                       "African Ghanaian"), c("Yes, EUR", "No", "No"))
#fixed things so that things like African-American is changed to African American
valid_subtypes$race_ethnicity_fixed <- plyr::mapvalues(valid_subtypes$`RACE/ETHNICITY`, c("African-American", "European-American",
                                                                                                          "Mexico Hispanic", "ASIAN"), c("African American", "European American",
                                                                                                                                         "Mexican", "Asian" ))
sub_add$race_ethnicity_fixed <- sub_add$RACE.ETHNICITY # <- this got changed from RACE/ETHNICITY in the function call
colnames(sub_add) <- colnames(valid_subtypes)
vali_prep <- valid_subtypes |> dplyr::select(SAMPLE_ID, subtype_assigned_by_us,race_ethnicity_fixed,ancestry_assigned_by_us, STUDY_ID )
colnames(vali_prep) <- c("sample_id", "subtype", "subcontinental_ancestry","continental_ancestry", "study")

sub_prep <- sub_add |> dplyr::select(SAMPLE_ID, subtype_assigned_by_us,race_ethnicity_fixed,ancestry_assigned_by_us, STUDY_ID )
colnames(sub_prep) <- c("sample_id", "subtype", "subcontinental_ancestry","continental_ancestry", "study")

# validation_metadata will be the rda object saved at end of script
validation_metadata <- rbind(vali_prep, sub_prep)
dim(validation_metadata)

##------------- GET EXPRESSION FOR STUDIES -------------##

studies <- c("Breast_Cancer_SMC2018","GSE101927","GSE142258","GSE46581",
             "GSE142731", "GSE59590","GSE2109","GSE62502", "GSE75678", "GSE225846", "GSE86374" , "GSE33095", "GSE37751")

length(studies)
suffix <- "_expr.txt"
file_list <- here::here("data-raw","validation",studies,paste0(studies,suffix))
## - read in all expression matrices
dat <- lapply(file_list, readr::read_tsv)

## - checking all matrices were loaded
length(studies) == length(dat)
library(dplyr)
## - combine all dataframes by SYMBOL id
expr_dat <- dat |> purrr::reduce(inner_join)
dim(expr_dat) #12052  3745 -- after removal of proteomic study: 12157 x 3623
expr_dat[1:5,1:5]

## - make sure there are no duplicated - ##
expr_dat[duplicated(expr_dat$SYMBOL),]

## move symbols to row names
expr_dat <- tibble::column_to_rownames(expr_dat, "SYMBOL")
## - only keep samples that are also in the validation metadata
expr_cut <- expr_dat[,colnames(expr_dat) %in% validation_metadata$sample_id]
dim(expr_cut) #12052  1342 -- after revmoal of proteogenomic study 12157 x 1239

## - same as above but opposite - keeping samples in expressino set

head(validation_metadata)
dim(validation_metadata)
validation_metadata[1:5,1:5]

table( validation_metadata$continental_ancestry,validation_metadata$subtype)
## after cuts:
#
# Basal Luminal
# African    331      46
# American    50     166
# Asian       55     212
# European   112     250
# None        14       3

## prematurley saving validation data and matrix just in case
#validation_met_trimmed <- validation_cut # this is used in the subsaharan study
#usethis::use_data(validation_met_trimmed, overwrite = TRUE)


val_sample <- validation_metadata |> dplyr::select(sample_id, subtype)
colnames(val_sample) <- c("sample_id", "subtype")
expr_temp <- expr_cut
expr_temp[1:5,1:5]
expr_temp <- t(expr_temp)
expr_temp <- tibble::rownames_to_column(as.data.frame(expr_temp), "SAMPLE_ID")
validation_raw_expr <- merge(expr_temp, val_sample, by.x = "SAMPLE_ID", by.y = "sample_id")

####Prep Multi-Validation Expression Data
expr_dat <- readr::read_tsv(here::here("data-raw","validation","GSE211167","GSE211167_tpm_matrix.txt"))
#-- replace gnees with alias in the validation set, then only keep common genes between the datasets
names(expr_dat)[names(expr_dat) == "HID1"] <- "C17orf28"
names(expr_dat)[names(expr_dat) == "CYP2B7P"] <- "CYP2B7P1"
names(expr_dat)[names(expr_dat) == "MISP"] <- "C19orf21"
names(expr_dat)[names(expr_dat) == "BRINP3"] <- "FAM5C"
names(expr_dat)[names(expr_dat) == "CT83"] <- "CXorf61"
names(expr_dat)[names(expr_dat) == "BPIFB1"] <- "C20orf114"
names(expr_dat)[names(expr_dat) == "FDCSP"] <- "C4orf7"
names(expr_dat)[names(expr_dat) == "SOWAHA"] <- "ANKRD43"
names(expr_dat)[names(expr_dat) == "BRINP3"] <- "FAM5C"
names(expr_dat)[names(expr_dat) == "MS4A8"] <- "MS4A8B"
names(expr_dat)[names(expr_dat) == "NEURL1"] <- "NEURL"
names(expr_dat)[names(expr_dat) == "NSG2"] <- "HMP19"
names(expr_dat)[names(expr_dat) == "GPC1"] <- "GPC1-AS1"
expr_dat$subtype <- "Basal"
#expr_dat$subtype <- as.factor(expr_dat$subtype)
expr_dat <- tibble::column_to_rownames(expr_dat,"Label")
validation_raw_expr <- tibble::column_to_rownames(validation_raw_expr, "SAMPLE_ID")
common_genes <- intersect(colnames(expr_dat), colnames(validation_raw_expr)) #12126
length(common_genes)

GSE_dat <- validation_raw_expr[,colnames(validation_raw_expr) %in% common_genes]
GSE_dat <- GSE_dat[,common_genes]
sub_dat <- expr_dat[,colnames(expr_dat) %in% common_genes]
sub_dat <- sub_dat[,common_genes]
all.equal(colnames(sub_dat), colnames(GSE_dat))
all.equal(colnames(GSE_dat), common_genes)
validation_expression <- dplyr::bind_rows(GSE_dat, sub_dat)
dim(validation_expression)
validation_expression[1:5,1:5]
validation_expression$subtype <- as.factor(validation_expression$subtype)
usethis::use_data(validation_expression, overwrite = FALSE)

## - get TCGA training matrix with the common genes
tcga_validation_training <- expression_breast[,colnames(expression_breast) %in% common_genes]
dim(tcga_validation_training)
usethis::use_data(tcga_validation_training, overwrite = FALSE)

# TODO  MOVE THESE ALL DOWN TO WHEN WE SAVE EXPRESSION AMTRIX
validation_metadata <- validation_metadata[validation_metadata$sample_id %in% rownames(validation_expression),]
# final validation metadata object:
usethis::use_data(validation_metadata  , overwrite = FALSE)

valid_subtypes <- valid_subtypes[valid_subtypes$SAMPLE_ID%in% rownames(validation_expression),]
all.equal(colnames(valid_subtypes), colnames(sub_add))
validation_metadata_extended <- rbind(valid_subtypes, sub_add)
#we dont ue GSE142102
#validation_metadata_extended  <- validation_metadata_extended[validation_metadata_extended$STUDY_ID != "GSE142102",]
usethis::use_data(validation_metadata_extended  , overwrite = FALSE)



