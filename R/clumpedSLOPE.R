#' Clumped SLOPE
#'
#' @export
#' @inheritParams clumpProcedure
#' @param fdr, False Discovery Rate for SLOPE
clumpedSLOPE <- function(y, X, rho = 0.3, fdr = 0.2){
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
  # filter all snps with no information
  clumpedSNPs <- clumpProcedure(y, X, 0.3)
  X2 <- clumpedSNPs$SNPs

  slopeResult <- SLOPE(X = X2, y = y, fdr = 0.05)

  selectedSNPs <- unlist(clumpedSNPs$SNPnumber)[slopeResult$selected]

  selectedSNPs
}

