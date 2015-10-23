#' clumpingResult class
#'
#' A result of procedure for snp clumping produced by \code{clumpProcedure} and
#' \code{clumpProcedure2}.
#'
#' @details Always a named list of ten elements
#' \enumerate{
#' \item \code{X} matrix of snps. One snp representative per each clump
#' \item \code{y} phenotype
#' \item \code{SNPnumber} vector of column number in SNP matrix \code{X_all} of
#' clumps representative
#' \item \code{SNPclumps} list of vectors of column numbers in SNP matrix \code{X_all}
#' for each clump member
#' \item \code{X_info} SNP info from .map file. Obtained from screening procedure.
#' \item \code{selectedSnpsNumbers} Contains information on which rows of \code{X_info}
#' matrix refer to selected clump representatives
#' \item \code{X_all} X matrix that was input for clumping procedure
#' \item \code{numberOfSnps} Number of SNPs in input for screening procedure
#' \item \code{selectedSnpsNumbersScreening} Contains information on which rows of
#' \code{X_info} matrix refer to SNPs that are input in clumping procedure
#' \item \code{pValMax} p-value used in screening procedure
#' }
#' @seealso \code{\link{screeningResult}} \code{\link{clumpProcedure2}}
#' @name clumpingResult
NULL