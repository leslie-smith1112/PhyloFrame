#' @title TCGA Breast Cancer
#'
#' @description A TCGA breast cancer dataset for binary classification task of Basal vs Luminal subtypes.
#'
#' @format A data frame with 842 rows, representing patient samples, and 20512 columns representing genes and subtype.
#' \describe{
#'   \item{samples}{Samples from TCGA BRCA dataset, only basal and luminal samples included. Represent data frame rows.}
#'   \item{genes}{Genes from TCGA BRCA dataset. Represent data frame columns.}
#'   \item{subtype}{Defines the samples subtype. Represent the last column in the data frame.}
#' }
#' @source <https://github.com/leslie-smith1112>
"expression_breast"


#' @title TCGA Thyroid Cancer
#'
#' @description A TCGA thyroid cancer dataset for binary classification task of M0 vs MX (metastasis).
#'
#' @format A data frame with 436 rows, representing patient samples, and 20512 columns representing genes and metastasis.
#' \describe{
#'   \item{samples}{Samples from TCGA THCA dataset, only M0 and MX samples included. Represent data frame rows.}
#'   \item{genes}{Genes from TCGA THCA dataset. Represent data frame columns.}
#'   \item{subtype}{Defines the samples subtype. Represents the last column in the data frame.}
#' }
#' @source <https://github.com/leslie-smith1112>
"expression_thyroid"


#' @title TCGA Uterine Cancer
#'
#' @description A TCGA uterine cancer dataset for binary classification task of Endometrioid vs Serous subtypes.
#'
#' @format A data frame with 491 rows, representing patient samples, and 20512 columns representing genes and subtype.
#' \describe{
#'   \item{samples}{Samples from TCGA UCEC dataset, only endometrioid and serous samples included. Represent data frame rows.}
#'   \item{genes}{Genes from TCGA UCEC dataset. Represent data frame columns.}
#'   \item{subtype}{Defines the samples subtype. Represents the last column in the data frame.}
#' }
#' @source <https://github.com/leslie-smith1112>
"expression_uterine"


#' @title TCGA Sample Ancestry
#'
#' @description The TCGA sample ancestry assignments for all TCGA cancers from Carrot-Zhang et al paper: Comprehensive Analysis of Genetic Ancestry and Its Molecular Correlates in Cancer.
#'
#' @format A data frame with 10586 rows, representing patient samples, and 13 columns representing the ancestry calls for each patient.
#' \describe{
#'   \item{samples}{Samples from all TCGA datasets.}
#'   \item{ancestry}{Ancestry is assigned through consensus ancestry calling as described by Carrot-Zhang et al in Comprehensive Analysis of Genetic Ancestry and Its Molecular Correlates in Cancer.}
#' }
#' @source <https://github.com/leslie-smith1112>
"ancestry"


#' @title Sample Ancestry for BRCA samples included in analysis
#'
#' @description The TCGA sample ancestry assignments from Carrot-Zhang et al for only samples included in our TCGA BRCA datasets.
#'
#' @format A data frame with 842 rows, representing patient samples, and  columns representing the ancestry calls.
#' \describe{
#'   \item{samples}{Samples from only TCGA BRCA samples used in our study.}
#'   \item{ancestry}{Ancestry is assigned through consensus ancestry calling as described by Carrot-Zhang et al in Comprehensive Analysis of Genetic Ancestry and Its Molecular Correlates in Cancer.}
#' }
#' @source <https://github.com/leslie-smith1112>
"ancestry_breast"

#' @title Sample Ancestry for THCA samples included in analysis
#'
#' @description The TCGA sample ancestry assignments from Carrot-Zhang et al for only samples included in our TCGA THCA datasets.
#'
#' @format A data frame with 436 rows, representing patient samples, and 2 columns representing the ancestry calls.
#' \describe{
#'   \item{samples}{Samples from only TCGA THCA samples used in our study.}
#'   \item{ancestry}{Ancestry is assigned through consensus ancestry calling as described by Carrot-Zhang et al in Comprehensive Analysis of Genetic Ancestry and Its Molecular Correlates in Cancer.}
#' }
#' @source <https://github.com/leslie-smith1112>
"ancestry_thyroid"

#' @title Sample Ancestry for UCEC samples included in analysis
#'
#' @description The TCGA sample ancestry assignments from Carrot-Zhang et al for only samples included in our TCGA UCEC datasets.
#'
#' @format A data frame with 491 rows, representing patient samples, and 2 columns representing the ancestry calls.
#' \describe{
#'   \item{samples}{Samples from only TCGA UCEC samples used in our study.}
#'   \item{ancestry}{Ancestry is assigned through consensus ancestry calling as described by Carrot-Zhang et al in Comprehensive Analysis of Genetic Ancestry and Its Molecular Correlates in Cancer.}
#' }
#' @source <https://github.com/leslie-smith1112>
"ancestry_uterine"


#' @title
#'
#' @description
#'
#' @format
#' \describe{
#'   \item{samples}{Samples from only TCGA UCEC samples used in our study.}
#'   \item{ancestry}{Ancestry is assigned through consensus ancestry calling as described by Carrot-Zhang et al in Comprehensive Analysis of Genetic Ancestry and Its Molecular Correlates in Cancer.}
#' }
#' @source <https://github.com/leslie-smith1112>
"ancestry_uterine"







