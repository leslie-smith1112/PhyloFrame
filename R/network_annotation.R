

#' Humanbase Network Annotation
#'
#' @param network
#'
#' @return
#' @export
#'
#' @examples
annotate_network <- function(network_name){
  network_file <- here::here("data-raw", network_name)
  network <- readr::read_tsv(network_file, col_names = FALSE)
  annot.df <- data.frame("Symbols" = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = as.character(network$X2), column = "SYMBOL", keytype = "ENTREZID"), network)
  annot.df <- annot.df %>% dplyr::select(-X2)
  colnames(annot.df) <- c("Gene2", "X1", "X3")
  annot.df <- data.frame("Symbols" = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = as.character(annot.df$X1), column = "SYMBOL", keytype = "ENTREZID"), annot.df)
  annot.df <- annot.df %>% dplyr::select(-X1)
  colnames(annot.df) <- c("Gene1","Gene2","Connection")
  rownames(annot.df) <- NULL

  out_file <- here::here("data-raw", paste0(network_name, "_symbol.tsv"))

  write.table(annot.df, file = out_file,
              sep = '\t', col.names = TRUE, row.names = FALSE)
}




