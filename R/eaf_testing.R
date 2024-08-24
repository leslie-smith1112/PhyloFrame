## Testing EAF ##
## comparing the old EAF code to the result of the new EAF code to make sure nothing has changed.
# network <- readr::read_tsv(here::here("data-raw","mammary_epithelium_symbol.top"))
# network <- network[,-4]
# eaf <- readr::read_tsv(here::here("data-raw", "mean_enhancedAF_exome.tsv"))
#
# #select only the ancestries in
# colnames(eaf)
# library(magrittr)
# cut_eaf <- eaf %>% dplyr::select(chrom, position, rs_id, ref_allele,
#                                  alt_allele, afr, sas, amr, eas, nfe,
#                                  asj, fin, gene, consequence, type, distance)
# sig <- readr::read_tsv("/home/leslie.smith1/blue_kgraim/leslie.smith1/Repositories/PhyloFrame/results/TCGA_changed_samples_order/model_runs/phyloFrame/base/eur/model_1/model_1_all_sig.txt")
# head(sig)
# gene_sig <- sig$Variable
#
# expression_breast
# new_function <- phyloFrame_gene_selection_tester(network,gene_sig,expression_breast,cut_eaf)
# new_function
# old_function <- order_frequencies_tester(network, gene_sig, expression_breast, cut_eaf)
# old_function
#
# ## changing the gnomad 4 selectiom ##

EAF_selection_tester <- function(eafs, train_expression, curr_ancestry, alt_anc_1,alt_anc_2, alt_anc_3,
                          alt_anc_4, alt_anc_5, alt_anc_6){
  # genes_keep <- 30
  # upper_af <- 1
  # lower_af <- 0.001
  # upper_bound <- 10000
  # lower_bound <- 0
  curr_ancestry
  anc_snps <- eafs[((eafs[,curr_ancestry] > 0) & (eafs[,alt_anc_1] < 0) & (eafs[,alt_anc_2] < 0) & (eafs[,alt_anc_3] < 0) & (eafs[,alt_anc_4] < 0) &
                      (eafs[,alt_anc_5] < 0) & (eafs[,alt_anc_6] < 0)) ,]
  #print(dim(anc_snps))
  ##possibly remove this
  # anc_genes <- table(anc_snps$gene)
  # anc_mut_order <- anc_genes[order(anc_genes, decreasing = TRUE)]
  # top_mut <- anc_mut_order[anc_mut_order < 10000 & anc_mut_order > 0]
  # valid_muts <- anc_snps[anc_snps$gene %in% names(top_mut),] #mutaiton bounds very large - no excluded genes
  # only keeps snps within enhanced allele frequency range
  valid_af <- anc_snps[(anc_snps[,curr_ancestry] > 0.001) & (anc_snps[,curr_ancestry] < 1),]
  #print(dim(valid_af))
  valid_snps <- na.omit(unique(valid_af$gene))
  valid_snps <- c(valid_snps, "subtype")
  anc_expression <- train_expression[,colnames(train_expression) %in% valid_snps]
  #print(dim(anc_expression))
  valid_genes <- get_variable_genes(anc_expression, 10)
  valid_genes <- unique(valid_genes[!(valid_genes == "subtype")])
  return(valid_genes)
}


phyloFrame_gene_selection_tester <- function(network, base_genes, train_expression, eafs){
  network_genes <- network_walk(network, base_genes)
  #CHROM	POS	RSID	REF	ALT	AFR_EAF	amr	asj	eas	fin	MID_EAF	nfe	SAS_EAF	GENE	SYMBOL	IMPACT	CONSEQUENCE	ALLELE	DISTANCE
  cut_eaf <- eafs[eafs$gene %in% network_genes,]
  #cut_eaf <- eafs[eafs$gene %in% network_genes,]
  afr <- EAF_selection_tester(cut_eaf, train_expression,"afr", "amr","asj", "eas","fin", "nfe", "sas")
  names(afr) <- rep("AFR", length(afr))
  message("Completed AFR")
  amr <- EAF_selection_tester(cut_eaf, train_expression,"amr", "afr","asj", "eas", "fin", "nfe", "sas")
  names(amr) <- rep("AMR", length(amr))
  asj <- EAF_selection_tester(cut_eaf, train_expression,"asj", "afr","amr", "eas", "fin", "nfe", "sas")
  names(asj) <- rep("ASJ", length(asj))
  eas <- EAF_selection_tester(cut_eaf, train_expression,"eas", "afr","amr", "asj", "fin", "nfe", "sas")
  names(eas) <- rep("EAS", length(eas))
  message("Completed EAS")
  fin <- EAF_selection_tester(cut_eaf, train_expression,"fin", "afr","amr", "asj", "eas", "nfe", "sas")
  names(fin) <- rep("FIN", length(fin))
  nfe <- EAF_selection_tester(cut_eaf, train_expression,"nfe", "afr","amr", "asj", "eas", "fin", "sas")
  names(nfe) <- rep("NFE", length(nfe))
  sas <- EAF_selection_tester(cut_eaf, train_expression,"sas", "afr","amr", "asj", "eas", "fin","nfe")
  names(sas) <- rep("SAS", length(sas))
  message("Completed SAS")
  ancestry_genes <- c(afr,amr,asj,eas,fin,nfe,sas)
  ## - make sure selected genes are in our expression matrix - ##
  valid_ancestry_genes <- ancestry_genes[(ancestry_genes %in% colnames(train_expression)) & (!ancestry_genes %in% base_genes)]
  message("Ancestry-specific variants identified.")
  return(valid_ancestry_genes)

}


### old filtering ##
order_frequencies_tester <- function(network, base_genes, expression, eafs){
  #### keep only include genes that are relevant in network ####
  # include genes = genes from netowrk walk
  network_genes <- network_walk(network, base_genes)
  #CHROM	POS	RSID	REF	ALT	AFR_EAF	AMR_EAF	ASJ_EAF	EAS_EAF	FIN_EAF	MID_EAF	NFE_EAF	SAS_EAF	GENE	SYMBOL	IMPACT	CONSEQUENCE	ALLELE	DISTANCE
  #cut_eaf <- eafs[eafs$SYMBOL %in% network_genes,]
  no_sex <- eafs[eafs$gene %in% network_genes,]
  keep.num <- 10
  upper.bound <- 1
  lower.bound <- 0.001
  mutationburden.upper <- 10000 # these parameters were defined for testing - threshold was made large enough no genes get excluded from this.
  mutationburden.lower <- 0

  ########################################################################
  ### For each ancestry, get range where right peak of enhanced frequencies
  ### is highest and keep only genes with snps in that range. Add genes to
  ### list holding all ancestry selected genes.
  #######################################################################

  ########################################################################
  ### AFR ###
  #######################################################################
  #--- get enhanced genes for each ancestry ---#
  afr_snps <- no_sex[no_sex$afr > 0 & no_sex$sas < 0 &  no_sex$amr < 0 & no_sex$eas < 0  & no_sex$nfe < 0 & no_sex$fin < 0 & no_sex$asj < 0 , ]
  afr.t <- table(afr_snps$gene)
  afr.t <- afr.t[order(afr.t, decreasing = TRUE)]
  afr.k <- afr.t[afr.t < mutationburden.upper & afr.t > mutationburden.lower]
  afr.nsnps <- afr_snps[afr_snps$gene %in% names(afr.k),] # mutaiton bounds very large - no excluded genes
  afr.snps.k <- afr.nsnps[(afr.nsnps$afr > lower.bound & afr.nsnps$afr < upper.bound),] # only keeps snps within enhanced allele frequency range

  #afr.snps.k <- afr_snps[(afr_snps$afr > lower.bound & afr_snps$afr < upper.bound),]
  afr.keep <- na.omit(afr.snps.k$gene)
  afr.keep <- c(afr.keep, "subtype")
  temp.expr <- expression[,colnames(expression) %in% afr.keep]
  afr.keep <- get_variable_genes(temp.expr, keep.num)
  dat <- afr.keep

  ########################################################################
  ### EAS ###
  #######################################################################
  # for ancestries with subancestries, subancestries are excluded from calculation as many mutations are commonly enahnced between them
  eas_snps <- no_sex[no_sex$eas > 0 & no_sex$sas < 0  & no_sex$amr < 0 & no_sex$afr < 0  & no_sex$nfe < 0 & no_sex$fin < 0 & no_sex$asj < 0 , ]
  eas.t <- table(eas_snps$gene)
  eas.t <- eas.t[order(eas.t, decreasing = TRUE)]
  eas.k <- eas.t[eas.t < mutationburden.upper & eas.t > mutationburden.lower]
  eas.nsnps <- eas_snps[eas_snps$gene %in% names(eas.k),]
  eas.snps.k <- eas.nsnps[(eas.nsnps$eas > lower.bound & eas.nsnps$eas < upper.bound),]
  #eas.snps.k <- eas_snps[(eas_snps$eas > lower.bound & eas_snps$eas < upper.bound),]
  eas.keep <- na.omit(eas.snps.k$gene)
  eas.keep <- c(eas.keep, "subtype")
  temp.expr <- expression[,colnames(expression) %in% eas.keep]
  eas.keep <- get_variable_genes(temp.expr, keep.num)

  dat <- c(dat,eas.keep)

  # ########################################################################
  # ### NFE ###
  # #######################################################################
  # for ancestries with subancestries, subancestries are excluded from calculation as many mutations are commonly enahnced between them
  nfe_snps <- no_sex[no_sex$nfe > 0 & no_sex$sas < 0 & no_sex$amr < 0 & no_sex$afr < 0 & no_sex$eas < 0 & no_sex$asj < 0 &  no_sex$fin < 0, ]
  nfe.t <- table(nfe_snps$gene)
  nfe.t <- nfe.t[order(nfe.t, decreasing = TRUE)]
  nfe.k <- nfe.t[nfe.t < mutationburden.upper & nfe.t > mutationburden.lower]
  nfe.nsnps <- nfe_snps[nfe_snps$gene %in% names(nfe.k),]
  nfe.snps.k <- nfe.nsnps[(nfe.nsnps$nfe > lower.bound & nfe.nsnps$nfe < upper.bound),]
  #nfe.snps.k <- nfe_snps[(nfe_snps$nfe > lower.bound & nfe_snps$nfe < upper.bound),]
  nfe.keep <- na.omit(nfe.snps.k$gene)
  nfe.keep <- c(nfe.keep, "subtype")
  temp.expr <- expression[,colnames(expression) %in% nfe.keep]
  nfe.keep <- get_variable_genes(temp.expr, keep.num)

  dat <- c(dat,nfe.keep)

  # ########################################################################
  # ### AMR ###
  # #######################################################################
  amr_snps <- no_sex[no_sex$amr > 0 & no_sex$sas < 0 &  no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$eas < 0 & no_sex$fin < 0 & no_sex$asj < 0 , ]

  amr.t <- table(amr_snps$gene)
  amr.t <- amr.t[order(amr.t, decreasing = TRUE)]
  amr.k <- amr.t[amr.t < mutationburden.upper & amr.t > mutationburden.lower]
  amr.nsnps <- amr_snps[amr_snps$gene %in% names(amr.k),]
  amr.snps.k <- amr.nsnps[(amr.nsnps$amr > lower.bound & amr.nsnps$amr < upper.bound),]
  #amr.snps.k <- amr_snps[(amr_snps$amr > lower.bound & amr_snps$amr < upper.bound),]
  amr.keep <- na.omit(amr.snps.k$gene)
  amr.keep <- c(amr.keep, "subtype")
  temp.expr <- expression[,colnames(expression) %in% amr.keep]
  amr.keep <- get_variable_genes(temp.expr, keep.num)

  dat <- c(dat,amr.keep)


  # ########################################################################
  # ### FIN ###
  # #######################################################################
  fin_snps <- no_sex[no_sex$fin > 0 & no_sex$sas < 0 & no_sex$afr < 0 &  no_sex$nfe < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0, ]
  fin.t <- table(fin_snps$gene)
  fin.t <- fin.t[order(fin.t, decreasing = TRUE)]
  fin.k <- fin.t[fin.t < mutationburden.upper & fin.t > mutationburden.lower]
  fin.nsnps <- fin_snps[fin_snps$gene %in% names(fin.k),]
  fin.snps.k <- fin.nsnps[(fin.nsnps$fin > lower.bound & fin.nsnps$fin < upper.bound),]
  #fin.snps.k <- fin_snps[(fin_snps$fin > lower.bound & fin_snps$fin < upper.bound),]
  fin.keep <- na.omit(fin.snps.k$gene)# top most enhanced genes
  fin.keep <- c(fin.keep, "subtype")
  temp.expr <- expression[,colnames(expression) %in% fin.keep]
  fin.keep <- get_variable_genes(temp.expr, keep.num)

  dat <- c(dat,fin.keep)


  # ########################################################################
  # ### SAS ###
  # #######################################################################
  sas_snps <- no_sex[no_sex$sas > 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$fin < 0 & no_sex$eas < 0 & no_sex$amr < 0 & no_sex$asj < 0, ]
  sas.t <- table(sas_snps$gene)
  sas.t <- sas.t[order(sas.t, decreasing = TRUE)]
  sas.k <- sas.t[sas.t < mutationburden.upper & sas.t > mutationburden.lower]
  sas.nsnps <- sas_snps[sas_snps$gene %in% names(sas.k),]
  sas.snps.k <- sas.nsnps[(sas.nsnps$sas > lower.bound & sas.nsnps$sas < upper.bound),]
  #sas.snps.k <- sas_snps[(sas_snps$sas > lower.bound & sas_snps$sas < upper.bound),]
  sas.keep <- na.omit(sas.snps.k$gene)
  sas.keep <- c(sas.keep, "subtype")
  temp.expr <- expression[,colnames(expression) %in% sas.keep]
  sas.keep <- get_variable_genes(temp.expr, keep.num)

  dat <- c(dat,sas.keep)


  # ########################################################################
  # ### asj ###
  # #######################################################################
  asj_snps <- no_sex[no_sex$asj > 0 & no_sex$sas < 0 & no_sex$nfe < 0 & no_sex$afr < 0 & no_sex$fin < 0 & no_sex$eas < 0 & no_sex$amr < 0 , ]
  asj.t <- table(asj_snps$gene)
  asj.t <- asj.t[order(asj.t, decreasing = TRUE)]
  asj.k <- asj.t[asj.t < mutationburden.upper & asj.t > mutationburden.lower]
  asj.nsnps <- asj_snps[asj_snps$gene %in% names(asj.k),]
  asj.snps.k <- asj.nsnps[(asj.nsnps$asj > lower.bound & asj.nsnps$asj < upper.bound),]
  #asj.snps.k <- asj_snps[(asj_snps$asj > lower.bound & asj_snps$asj < upper.bound),]
  asj.keep <- na.omit(asj.snps.k$gene)
  asj.keep <- c(asj.keep, "subtype")
  temp.expr <- expression[,colnames(expression) %in% asj.keep]
  asj.keep <- get_variable_genes(temp.expr, keep.num)

  dat <- c(dat,asj.keep)
  dat <- dat[dat %in% "subtype" == FALSE] # get rid of subtype
  valid.keep <- unique(dat)


  return(valid.keep)
}











