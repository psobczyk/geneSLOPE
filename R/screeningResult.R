#' screeningResult class
#'
#' A result of procedure for snp clumping produced by \code{readBigSNPs}
#'
#' @details Always a named list of ten elements
#' \enumerate{
#' \item \code{X} matrix of snps that passed screening
#' \item \code{y} phenotype
#' \item \code{X_info} SNP info from .map file
#' \item \code{numberOfSnps} Number of SNPs in input for screening procedure
#' \item \code{selectedSnpsNumbers} Contains information on which rows of \code{X_info}
#' matrix refer to snps that passed screening
#' \item \code{pValMax} p-value used in screening procedure
#' }
#' @seealso \code{\link{clumpingResult}} \code{\link{readBigSNPs}}
#' @name screeningResult
NULL
