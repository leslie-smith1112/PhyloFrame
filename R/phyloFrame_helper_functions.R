#PhyloFrame helper functions

#' Get Variable Genes
#'
#' @param training_expr The matrix you would like to get the most variable genes from. Expects samples as rows and genes as columns.
#' @param num_genes How many variable genes you want.
#'
#' @return Returns the top n most variable genes in the current training expression matrix.
#' @export
#'
#' @examples most_variable_genes <- get_variable_genes(expression_matrix, 10000)
get_variable_genes <- function(training_expr, num_genes){
  if(!is.null(dim(training_expr))){
    #get rid of subtype so variance operation doesnt give warning
    training_expr <- subset(training_expr, select= -subtype)
      #base::subset(training_expr, select= -subtype)
    # <- training_expr[,-ncol(training_expr)]
    expr_variance <- apply(training_expr,2,var)
    variance <- expr_variance[order(expr_variance, decreasing = TRUE)]
    top_genes <- variance[1:num_genes]
    top_genes <- names(top_genes)
    top_genes <- c(top_genes, "subtype")
    return(top_genes)
  }else{
    stop()
  }
}

#' Network walk
#'
#' @param network The tissue specific functional interaction network from Humanbase.
#' @param base_genes The genes from the initial lasso run.
#'
#' @return Genes within 2 neighbors of the base_genes.
#' @export
#'
#' @examples genes <- network_walk(breast_network, lasso_genes)
network_walk <- function(network, base_genes){
  trimmed_graph <- network[((network$Connection >= 0.2) & (network$Connection < 0.51)),]
  ## -- some of the dz genes are not in network -> trim them out -- ##
  valid_genes <- base_genes[(base_genes %in% trimmed_graph$Gene1 | base_genes %in% trimmed_graph$Gene2)]
  graph_dat <- igraph::graph.data.frame(trimmed_graph, directed =  FALSE)
  ## couldn't just pass a list of nodes into the graph for some reason, so converted them into vertices.
  ## get graph vertices to get start nodes from signature
  vertices <- igraph::V(graph_dat)$name
  vertices <- vertices[vertices %in% valid_genes]
  neighborhood <- igraph::ego(graph_dat, order = 2, nodes = vertices, mode = "all")
  valid_nodes <- unlist(neighborhood)
  nodes <- unique(c(names(valid_nodes)),base_genes)# this is the gene list you need from the neighborhood.
  message("Completed finding gene signature neighborhood, moving to Exome allele frequencies.")
  return(nodes)
}



#' EAF selection
#'
#' @param network Tissue specific functional interaction network.
#' @param base_genes The genes from the initial lasso run.
#' @param train_expression Expression matrix being used in training.
#' @param eafs Enhanced allele frequency for each snp.
#'
#' @return
#' @export
#'
#' @examples ancestry_genes <- phyloFrame_gene_selection(tissue_network, lasso_genes, expression_mat, eaf_dat)
phyloFrame_gene_selection <- function(network, base_genes, train_expression, eafs){
  network_genes <- network_walk(network, base_genes)
  #CHROM	POS	RSID	REF	ALT	AFR_EAF	AMR_EAF	ASJ_EAF	EAS_EAF	FIN_EAF	MID_EAF	NFE_EAF	SAS_EAF	GENE	SYMBOL	IMPACT	CONSEQUENCE	ALLELE	DISTANCE
  cut_eaf <- eafs[eafs$SYMBOL %in% network_genes,]
  #cut_eaf <- cut_eaf[cut_eaf$SYMBOL %in% base_genes,]
  #cut_eaf <- eafs[eafs$gene %in% network_genes,]
  afr <- EAF_selection(cut_eaf, train_expression,"AFR_EAF", "AMR_EAF","ASJ_EAF", "EAS_EAF","FIN_EAF", "NFE_EAF", "MID_EAF", "SAS_EAF")
  names(afr) <- rep("AFR", length(afr))
  message("Completed AFR")
  amr <- EAF_selection(cut_eaf, train_expression,"AMR_EAF", "AFR_EAF","ASJ_EAF", "EAS_EAF", "FIN_EAF", "NFE_EAF", "MID_EAF", "SAS_EAF")
  names(amr) <- rep("AMR", length(amr))
  asj <- EAF_selection(cut_eaf, train_expression,"ASJ_EAF", "AFR_EAF","AMR_EAF", "EAS_EAF", "FIN_EAF", "NFE_EAF", "MID_EAF", "SAS_EAF")
  names(asj) <- rep("ASJ", length(asj))
  eas <- EAF_selection(cut_eaf, train_expression,"EAS_EAF", "AFR_EAF","AMR_EAF", "ASJ_EAF", "FIN_EAF", "NFE_EAF", "MID_EAF", "SAS_EAF")
  names(eas) <- rep("EAS", length(eas))
  message("Completed EAS")
  fin <- EAF_selection(cut_eaf, train_expression,"FIN_EAF", "AFR_EAF","AMR_EAF", "ASJ_EAF", "EAS_EAF", "NFE_EAF", "MID_EAF", "SAS_EAF")
  names(fin) <- rep("FIN", length(fin))
  nfe <- EAF_selection(cut_eaf, train_expression,"NFE_EAF", "AFR_EAF","AMR_EAF", "ASJ_EAF", "EAS_EAF", "FIN_EAF", "MID_EAF", "SAS_EAF")
  names(nfe) <- rep("NFE", length(nfe))
  mid <- EAF_selection(cut_eaf, train_expression,"MID_EAF", "AFR_EAF","AMR_EAF", "ASJ_EAF", "EAS_EAF", "FIN_EAF","NFE_EAF", "SAS_EAF")
  names(mid) <- rep("MID", length(mid))
  sas <- EAF_selection(cut_eaf, train_expression,"SAS_EAF", "AFR_EAF","AMR_EAF", "ASJ_EAF", "EAS_EAF", "FIN_EAF","NFE_EAF", "MID_EAF")
  names(sas) <- rep("SAS", length(sas))
  message("Completed SAS")


  ancestry_genes <- c(afr,amr,asj,eas,fin,nfe,mid,sas)
  ## - make sure selected genes are in our expression matrix - ##
  valid_ancestry_genes <- ancestry_genes[(ancestry_genes %in% colnames(train_expression)) & (!ancestry_genes %in% base_genes)]
  message("Ancestry-specific variants identified.")
  return(valid_ancestry_genes)

}
#for testing
#eafs <- readr::read_tsv(here::here("data-raw","testing_eaf.tsv"))

#' EAF filtering
#'
#' @param eafs The data.frame with snp eafs.
#' @param train_expression Expression matrix we are currently training on.
#' @param curr_ancestry The ancestry we are currently finding enhanced snps for.
#' @param alt_anc_1 Alternate ancestry. These are ancestries where we want the snps eaf to be negative.
#' @param alt_anc_2 Alternate ancestry. These are ancestries where we want the snps eaf to be negative.
#' @param alt_anc_3 Alternate ancestry. These are ancestries where we want the snps eaf to be negative.
#' @param alt_anc_4 Alternate ancestry. These are ancestries where we want the snps eaf to be negative.
#' @param alt_anc_5 Alternate ancestry. These are ancestries where we want the snps eaf to be negative.
#' @param alt_anc_6 Alternate ancestry. These are ancestries where we want the snps eaf to be negative.
#' @param alt_anc_7 Alternate ancestry. These are ancestries where we want the snps eaf to be negative.

#'
#' @return list of valid genes for the ancestry.
#' @export
#'
#' @examples EAF_selection(valid_genes, train_expression,"NFE", "AFR","AMR", "ASJ", "EAS", "FIN", "OTH", "SAS" )
EAF_selection <- function(eafs, train_expression, curr_ancestry, alt_anc_1,alt_anc_2, alt_anc_3,
                          alt_anc_4, alt_anc_5, alt_anc_6, alt_anc_7){
  # genes_keep <- 30
  # upper_af <- 1
  # lower_af <- 0.001
  # upper_bound <- 10000
  # lower_bound <- 0
  curr_ancestry
  valid_af <- eafs[(eafs[,curr_ancestry] >= 0.2) & (eafs[,curr_ancestry] <= 1),]

  anc_snps <- valid_af[((valid_af[,curr_ancestry] > 0) & (valid_af[,alt_anc_1] < 0) & (valid_af[,alt_anc_2] < 0) & (valid_af[,alt_anc_3] < 0) & (valid_af[,alt_anc_4] < 0) &
                     (valid_af[,alt_anc_5] < 0) & (valid_af[,alt_anc_6] < 0) & (valid_af[,alt_anc_7] < 0)) ,]
  valid_snps <- na.omit(unique(anc_snps$SYMBOL))
  valid_snps <- c(valid_snps, "subtype")
  ## start of edit
  # anc_snps1 <- unique(valid_af[(valid_af[,curr_ancestry] > 0),]$SYMBOL)
  # anc_snps2 <- unique(valid_af[(valid_af[,alt_anc_1] > 0),]$SYMBOL)
  # anc_snps3 <- unique(valid_af[(valid_af[,alt_anc_2] > 0),]$SYMBOL)
  # anc_snps4 <- unique(valid_af[(valid_af[,alt_anc_3] > 0),]$SYMBOL)
  # anc_snps5 <- unique(valid_af[(valid_af[,alt_anc_4] > 0),]$SYMBOL)
  # anc_snps6 <- unique(valid_af[(valid_af[,alt_anc_5] > 0),]$SYMBOL)
  # anc_snps7 <- unique(valid_af[(valid_af[,alt_anc_6] > 0),]$SYMBOL)
  # anc_snps8 <- unique(valid_af[(valid_af[,alt_anc_7] > 0),]$SYMBOL)
  #
  # all_genes <- c(anc_snps1, anc_snps2, anc_snps3, anc_snps4, anc_snps5, anc_snps6, anc_snps7, anc_snps8)
  # all_genes_count <- table(all_genes)
  # valid_genes <- names(all_genes_count[all_genes_count > 3])
  # valid_genes <- c(valid_genes, "subtype")
  #anc_expression <- train_expression[,colnames(train_expression) %in% valid_genes]
### end of edit

  anc_expression <- train_expression[,colnames(train_expression) %in% valid_snps]

  #print(dim(anc_expression))
  valid_genes <- get_variable_genes(anc_expression, 30)
  valid_genes <- unique(valid_genes[!(valid_genes == "subtype")])
  return(valid_genes)
}


#' PhyloFrame Driver
#'
#' @param curr_ancestry The name of the directory you want the results in. This will be within the "results" directory
#' @param curr_anc_expr The ancestry we are currently modeling.
#' @param cancer_type the current cancer we are working on
#' @param eafs data.frame with the enhanced allele frequencies for each snp
#' @param network tissue specific functional interaction network
#' @param i current model number, mostly for naming files.
#' @param benchmark_out_dir output directory for benchmark models
#' @param phyloFrame_out_dir output directory for phyloFrame models
#' @param base_genes_out_dir output directory for the baseline genes initally selected in first lasso run.
#'
#'
#' @return phyloframe and benchmark models.
#' @export
#'
#' @examples
phyloFrame <- function(curr_ancestry, curr_anc_expr, cancer_type, eafs, network, i, results_dir, single_batch){
  if(single_batch == FALSE){
    benchmark_out_dir <- here::here(results_dir, "benchmark",curr_ancestry, paste0("model_", i))
    phyloFrame_out_dir <- here::here(results_dir, "phyloFrame",curr_ancestry, paste0("model_", i))
    base_genes_out_dir <- here::here(results_dir, "phyloFrame", "base", curr_ancestry, paste0("model_", i))
    ## - read in current sample batch - ##
    current_batch <- readr::read_tsv(here::here("data-raw",paste0(cancer_type,"_samples"),curr_ancestry,paste0("samples",i,".tsv")))
    training_set <- curr_anc_expr[rownames(curr_anc_expr) %in% current_batch$x,]
  }else{
    benchmark_out_dir <- here::here(results_dir, "benchmark")
    phyloFrame_out_dir <- here::here(results_dir, "phyloFrame")
    base_genes_out_dir <- here::here(results_dir, "phyloFrame", "base")
    training_set <- curr_anc_expr
    current_batch <- NA
  }

  ## - make sure that the test set is not empty - for some ancestries (like admixed) we often dont
  ## - have enough samples to test on
  base_genes <- get_variable_genes(training_set, 10000) #change this to batch specific?
  if(is.null(base_genes) == FALSE){
    training_set <- training_set[,colnames(training_set) %in% base_genes]
  }

  lasso_model <- run_elastic_net(training_set, base_genes_out_dir, i, 0.8)
  #read in the genes from lasso model that will be that basis of our gene search
  lasso_genes <- read.table(paste0(base_genes_out_dir,"/model_",i,"_all_sig.txt"))
  lasso_genes <- lasso_genes[-1,]
  colnames(lasso_genes) <- c("Variable","Importance","Sign")
  lasso_genes <- lasso_genes$Variable
  phyloFrame_anc_genes <- phyloFrame_gene_selection(network, lasso_genes, training_set, eafs)
  ########################################################changed function########################################################
  #phyloFrame_anc_genes <- order_frequencies(network, lasso_genes, training_set, eafs)
  phyloFrame_genes <- c(phyloFrame_anc_genes, lasso_genes,"subtype")
  #write ancestry genes to file for FYI
  #anc <- names(phyloFrame_anc_genes)
  pf_write <- data.frame("GENE" = phyloFrame_anc_genes, "ANCESTRY" = names(phyloFrame_anc_genes))
  write.table(pf_write, paste0(phyloFrame_out_dir,"/ancestry_genes.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)

  #sample benchmark genes to match signature of phyloFrame

  benchmark_gene_length <- unique(phyloFrame_anc_genes) #get number of genes to add

  cosmic <- readr::read_tsv(here::here("Cosmic_Genes.tsv"))
  ## added for testing
  #base_genes_bm <- get_variable_genes(training_set, 10000)
  options <- base_genes[!(base_genes %in% lasso_genes)] # genes to choose from
  #options <- base_genes_bm[!(base_genes_bm %in% lasso_genes)]
  options <- options[options %in% cosmic$`Gene Symbol`]
  ## end edit


  benchmark_sampled_genes <- sample(options,length(benchmark_gene_length))
  benchmark_genes <- c(benchmark_sampled_genes, lasso_genes, "subtype")

  #call models
  benchmark_training <- training_set[,colnames(training_set) %in% benchmark_genes]
  benchmark_model <- run_elastic_net(benchmark_training, benchmark_out_dir, i, 0)
  phyloFrame_training <- training_set[,colnames(training_set) %in% phyloFrame_genes]
  phyloFrame_model <- run_elastic_net(phyloFrame_training, phyloFrame_out_dir, i, 0)
  models <- list(phyloFrame = phyloFrame_model, benchmark = benchmark_model, samples = current_batch, benchmark_dir = benchmark_out_dir, phyloFrame_dir =  phyloFrame_out_dir)
  return(models)
}

