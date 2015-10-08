#' Clumping procedure
#'
#' clumping of SNPs as preprocesing for SLOPE in GWAS
#'
#' @export
#' @param y numeric vector, phenotype
#' @param X numeric matrix, SNPs
#' @param rho numeric, minimal correlation between two SNPs to be
#' classified to the same clump
clumpProcedure <- function(y, X, rho = 0.3, pValMax = 0.1){
  if(rho>=1 | rho <= 0){
    stop("Rho has to be within range (0,1)")
  }
  if(length(y) != nrow(X)){
    stop("Length of y must match
         number of rows in X")
  }
  pVals <- foreach(i=1:ncol(X), .combine=c) %do% {
    summary(lm(y~X[,i]))$coefficients[2,4]
  }

  a <- order(pVals, decreasing = FALSE)
  notClumped <- rep(TRUE, length(a))
  clumps <- list()
  representatives <- list()

  i <- 1
  while(any(notClumped & (pVals[a[i]]<pValMax) )){
    idx = a[i]
    if(notClumped[idx]){
      clump <- apply(X[,notClumped, drop=FALSE], 2, cor, X[,idx])>rho
      clumps[[i]] <- which(notClumped)[clump]
      representatives[[i]] <- idx
      notClumped[ which(notClumped)[clump] ] <- FALSE
    }
    i = i+1
  }
  nullClumps <- sapply(representatives, is.null)
  representatives <- representatives[!nullClumps]
  clumps <- clumps[!nullClumps]
  result <- list(
    SNPs = X[,unlist(representatives)],
    SNPnumber = representatives,
    SNPclumps = clumps
  )
  return(result)
}


