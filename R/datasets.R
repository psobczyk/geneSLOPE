#' Matrix containing SNPs
#'
#' Subset of PLINK example data, available at PLINK
#' website
#'
#' @docType data
#' @usage SNPs
#' @format Large numeric matrix, 75510 elements
#' @keywords datasets
"SNPs"

#' Phenotype
#'
#' Artificial phenotype. Created as a linear combination
#' of small number of columns of \code{\link{SNPs}} matrix plus some noise
#'
#' @docType data
#' @usage aaaggregation
#' @format a numeric vector of length 90. Each element refers
#' to phenotype for one observations and is related to one
#' row in \code{\link{SNPs}} matrix
#' @keywords datasets
"y"
