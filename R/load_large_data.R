#' Load Enhanced Allele Frequency
#'
#' @return A dataframe holding the enhanced allele frequencies for variants.
#' @export
#'
#' @examples
#' eaf  <- load_EAF()
load_EAF <- function(){
  eaf <- readr::read_tsv(here::here("data-raw","gnomad4.1_annotated","gnomad.exomes.all.tsv"), col_names = TRUE)
  #eaf <- readr::read_tsv(here::here("data-raw", "gnomad4.1_annotated", "gnomad.exomes.all.head"), col_names = TRUE)
  #eaf <- readr::read_tsv(here::here("data-raw", "mean_enhancedAF_exome.tsv"))

  return(eaf)

}

## load functions for all networks:

#' Load Mammary Epithelium Functional Network
#'
#' @return A data frame holding the functional interaction network for mammary epithelium tissue from HumanBase annotated with gene symbols.
#' @export
#'
#' @examples
#' disease_network  <- load_breast_network()
load_breast_network <- function(){
  #network <- readr::read_tsv(here::here("data-raw","mammary_epithelium_symbol.top"))
  network <- readr::read_tsv(here::here("data-raw","mammary_epithelium_symbol.tsv"))
  network <- network %>% dplyr::select(-4)
  return(network)
}

#' Load Thyroid Gland Functional Network
#'
#' @return A data frame holding the functional interaction network for thyroid gland from HumanBase annotated with gene symbols.
#' @export
#'
#' @examples
#' disease_network  <- load_thyroid_network()
load_thyroid_network <- function(){
  network <- readr::read_tsv(here::here("data-raw", "thyroid_gland_symbol.tsv"))
  network <- network %>% dplyr::select(-4)
  return(network)
}

#' Load Uterine Endometrium Functional Network
#'
#' @return A data frame holding the functional interaction network for uterine endometrium from HumanBase annotated with gene symbols.
#' @export
#'
#' @examples
#' disease_network  <- load_uterine_network()
load_uterine_network <- function(){
  network <- readr::read_tsv(here::here("data-raw","uterine_endometrium_symbol.tsv"))
  network <- network %>% dplyr::select(-4)
  return(network)
}

# load_adrenal_gland_network <- function(){
#   network <- readr::read_tsv(here::here("data-raw","network_downloads_annotated","adrenal_gland_symbol.tsv"))
#   network <- network %>% dplyr::select(-4)
#   return(network)
# }
#
# load_blood_network <- function(){
#   network <- readr::read_tsv(here::here("data-raw","network_downloads_annotated","blood_symbol.tsv"))
#   network <- network %>% dplyr::select(-4)
#   return(network)
# }
#
# load_bone_marrow_network <- function(){
#   network <- readr::read_tsv(here::here("data-raw","network_downloads_annotated","bone_marrow_symbol.tsv"))
#   network <- network %>% dplyr::select(-4)
#   return(network)
# }
#
# load_brain_network <- function(){
#   network <- readr::read_tsv(here::here("data-raw","network_downloads_annotated","brain_symbol.tsv"))
#   network <- network %>% dplyr::select(-4)
#   return(network)
# }
#
# load_colon_network <- function(){
#   network <- readr::read_tsv(here::here("data-raw","network_downloads_annotated","colon_symbol.tsv"))
#   network <- network %>% dplyr::select(-4)
#   return(network)
# }
#
# load_esophagus_network <- function(){
#   network <- readr::read_tsv(here::here("data-raw","network_downloads_annotated","esophagus_symbol.tsv"))
#   network <- network %>% dplyr::select(-4)
#   return(network)
# }
#
# load_eye_network <- function(){
#   network <- readr::read_tsv(here::here("data-raw","network_downloads_annotated","eye_symbol.tsv"))
#   network <- network %>% dplyr::select(-4)
#   return(network)
# }
#
# load_gastrointestinal_tract_network <- function(){
#   network <- readr::read_tsv(here::here("data-raw","network_downloads_annotated","gastrointestinal_tract_symbol.tsv"))
#   network <- network %>% dplyr::select(-4)
#   return(network)
# }
#
# load_global_tract_network <- function(){
#   network <- readr::read_tsv(here::here("data-raw","network_downloads_annotated","global_symbol.tsv"))
#   network <- network %>% dplyr::select(-4)
#   return(network)
# }

#  "kidney", "liver", "lung", "lymph_node",
#            "pancreas", "peripheral_nervous_system", "prostate_gland", "skin", "stomach", "testis", "urinary_bladder", "uterine_cervix", "uterus"]



