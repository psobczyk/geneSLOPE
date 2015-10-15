#' GWAS with SLOPE
#'
#' Performs GWAS with SLOPE on data already read into R
#'
#' @export
#' @inheritParams clumpProcedure
#' @param fdr, False Discovery Rate for SLOPE
#' @param lambda lambda for SLOPE. See
#' @return data.frame with two columns \itemize{
#' \item snpName names of selected snps
#' \item snpNumber column numbers of selected snps in input matrix X
#' }
#' @examples
#' \donttest{
#' clumpedSLOPE(y, SNPs)
#' }
clumpedSLOPE <- function(y, X, rho = 0.3, fdr = 0.2, lambda="gaussian",
                         verbose = TRUE){
  if(rho>=1 | rho <= 0){
    stop("Rho has to be within range (0,1)")
  }
  if(fdr>=1 | fdr <= 0){
    stop("FDR has to be within range (0,1)")
  }
  if(length(y) != nrow(X)){
    stop("Length of y must match
         number of rows in X")
  }
  X <- apply(X, 2, replace_na_with_mean)
  message("Missing values were replaced by column mean")
  nonZeroSd <- apply(X, 2, sd)!=0
  X <- X[,nonZeroSd]
  if(sum(!nonZeroSd)>0)
    message(paste(sum(!nonZeroSd), "variables with zero variance were removed"))

  clumpedSNPs <- clumpProcedure(y, X, rho, verbose)
  X2 <- clumpedSNPs$SNPs

  slopeResult <- SLOPE(X = X2, y = y, fdr = fdr,
                       lambda = lambda[1:ncol(X2)])

  selectedSNPs <- unlist(clumpedSNPs$SNPnumber)[slopeResult$selected]
  selectedSNPs <- sort(selectedSNPs)
  cSLOPEResult <- data.frame(snpName=colnames(SNPs)[selectedSNPs],
                             snpNumber=selectedSNPs)
  cSLOPEResult
}

