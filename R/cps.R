#' Genome Wide Association Study with SLOPE
#'
#' Genome-wide association study performed with SLOPE.
#' In first step, SNPs are clumped according to their correlations and
#' distances. Then, SLOPE is performed on data where each clump has
#' one representative.
#'
#' @docType package
#' @name cps
#' @details Version: 0.29.0
#' @importFrom data.table fread
#' @importFrom tcltk tk_choose.files
#' @importFrom SLOPE SLOPE
#' @importFrom SLOPE create_lambda
#' @import ggplot2
#' @import bigmemory
#' @import biganalytics
#' @author Piotr Sobczyk
#'
#' Maintainer: Piotr Sobczyk \email{Piotr.Sobczyk@@pwr.edu.pl}
#' @examples
#' \donttest{
#' clumpedSLOPE(y, SNPs)
#' }
#' \donttest{
#' runExample()
#' }
NULL
