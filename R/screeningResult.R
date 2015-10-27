#' screeningResult class
#'
#' A result of procedure for snp clumping produced by \code{\link{readSNPs}}
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
#' @seealso \code{\link{clumpingResult}} \code{\link{readSNPs}}
#' @name screeningResult
NULL

#' Print screeningResult class object
#'
#' @param x screeningResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method print screeningResult
print.screeningResult <- function(x, ...){
  cat("Object of class screeningResult\n")
  cat("$X: Matrix of SNPs\n")
  cat("\t", nrow(x$X), " observations\n")
  cat("\t", ncol(x$X), " snps\n")
  cat("$y: a vector containing phenotype \n")
  cat("$X_info: : Matrix of SNPs info\n")
  cat("\t", ncol(x$X_info), " Info categories\n")
  cat("\t", nrow(x$X_info), " snps\n")
  cat("$numberOfSnps: Number of SNPs in data file\n")
  cat("$selectedSnpsNumbers: Column numbers of selected SNPs in original X matrix\n")
  cat("$pValMax: p-value threshold value used for screening\n")
}


#' Summary screeningResult class object
#'
#' @param object screeningResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method summary screeningResult
summary.screeningResult <- function(object, ...){
  cat("Object of class screeningResult\n")
  cat(nrow(object$X), " observations\n")
  cat(object$numberOfSnps, " SNPs were screened\n")
  cat(ncol(object$X), " snps had p-value under threshold in marginal test\n")
}
