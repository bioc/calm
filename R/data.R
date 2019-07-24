#' Psoriasis RNA-seq dataset
#'
#' A dataset containing the test statistics to analyze an RNA-seq study of
#' psoriasis.
#'
#' @format A dataset with the following vectors:
#' \describe{
#'   \item{zval}{16490 z-values of genes with matching microarray data}
#'   \item{len_gene}{16490 gene coding region length for zval}
#'   \item{tval_mic}{16490 matching microarray t-statistics}
#' }
#' @source Liang (2019), Empirical Bayes analysis of RNA sequencing
#' experiments with auxiliary information, to appear in Annals of Applied
#' Statistics;
#' @examples
#' data(pso)
#' dim(pso)
#' # total number of genes without matching microarray data
#' sum(is.na(pso$tval_mic))
"pso"

