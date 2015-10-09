#' Genome Wide Association Study with SLOPE
#'
#' Genome-wide association study performed with SLOPE.
#' In first step, SNPs are clumped according to their correlations and
#' distances. Then, SLOPE is performed on data where each clump has
#' one representative.
#'
#' @docType package
#' @name cps
#' @details Version: 0.22
#' @importFrom data.table fread
#' @importFrom adegenet read.PLINK
#' @importFrom tcltk tk_choose.files
#' @import foreach
#' @import bigmemory
#' @author Piotr Sobczyk
#'
#' Maintainer: Piotr Sobczyk \email{Piotr.Sobczyk@@pwr.edu.pl}
#'
NULL
