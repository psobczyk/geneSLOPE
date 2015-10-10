#' Clumping procedure
#'
#' clumping of SNPs as preprocesing for SLOPE in GWAS
#'
#' @param y numeric vector, phenotype
#' @param X numeric matrix, SNPs
#' @param rho numeric, minimal correlation between two SNPs to be
#' @param verbose logical, if TRUE (default) progress bar is shown
#' classified to the same clump
clumpProcedure <- function(y, X, rho = 0.3, verbose = TRUE){

  if(verbose){
    total = ncol(X)
    # create progress bar
    pb <- txtProgressBar(min = 0, max = total, style = 3)
  }

  if(rho>=1 | rho <= 0){
    stop("Rho has to be within range (0,1)")
  }
  if(length(y) != nrow(X)){
    stop("Length of y must match
         number of rows in X")
  }
  suma = sum((y-mean(y))^2)
  n = length(y) - 2
  pVals <- apply(X, 2, function(x) pValComp(x,y,n,suma))

  a <- order(pVals, decreasing = FALSE)
  notClumped <- rep(TRUE, length(a))
  clumps <- list()
  representatives <- list()

  i <- 1
  while(any(notClumped)){
    idx = a[i]
    if(notClumped[idx]){
      clump <- abs(apply(X[,notClumped, drop=FALSE], 2, cor, X[,idx]))>rho
      clumps[[i]] <- which(notClumped)[clump]
      representatives[[i]] <- idx
      notClumped[ which(notClumped)[clump] ] <- FALSE
    }
    i = i+1
    if(verbose)
      setTxtProgressBar(pb, i)
  }
  if(verbose)
    close(pb)
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


