#' Clumping procedure
#'
#' Clumping of SNPs previously read and screened by \code{\link{readSNPs}}
#'
#' @export
#' @param screenResult object of class screenResult
#' @param rho numeric, minimal correlation between two SNPs to be
#' @param verbose logical, if TRUE (default) progress bar is shown
#' classified to the same clump
#'
#' @return object of class \code{\link{clumpingResult}}
#'
clumpProcedure <- function(screenResult, rho = 0.3, verbose = TRUE){

  if(verbose){
    message("Clumping procedure has started. Depending on
            size of your data this may take several minutes.")
    total = sqrt(ncol(screenResult$X))
    # create progress bar
    pb <- txtProgressBar(min = 0, max = total, style = 3)
  }

  if(rho>=1 | rho <= 0){
    stop("Rho has to be within range (0,1)")
  }
  if(length(screenResult$y) != nrow(screenResult$X)){
    stop("Length of y must match
         number of rows in X")
  }
  suma = sum((screenResult$y-mean(screenResult$y))^2)
  n = length(screenResult$y) - 2
  pVals <- apply(screenResult$X, 2, function(x) pValComp(x,screenResult$y,n,suma))

  a <- order(pVals, decreasing = FALSE)
  notClumped <- rep(TRUE, length(a))
  clumps <- list()
  representatives <- list()

  i <- 1
  while(any(notClumped)){
    idx = a[i]
    if(notClumped[idx]){
      clump <- abs(apply(screenResult$X[,notClumped, drop=FALSE], 2, cor, screenResult$X[,idx]))>rho
      clumps[[i]] <- which(notClumped)[clump]
      representatives[[i]] <- idx
      notClumped[ which(notClumped)[clump] ] <- FALSE
    }
    i = i+1
    if(verbose)
      setTxtProgressBar(pb, sqrt(i))
  }
  if(verbose)
    close(pb)
  nullClumps <- sapply(representatives, is.null)
  representatives <- representatives[!nullClumps]
  clumps <- clumps[!nullClumps]
  result <- structure(
    list( X = screenResult$X[,unlist(representatives)],
          y = screenResult$y,
          SNPnumber = representatives,
          SNPclumps = clumps,
          X_info = screenResult$X_info,
          selectedSnpsNumbers = screenResult$selectedSnpsNumbers[unlist(representatives)],
          X_all = screenResult$X,
          numberOfSnps = screenResult$numberOfSnps,
          selectedSnpsNumbersScreening = screenResult$selectedSnpsNumbers,
          pVals = screenResult$pVals,
          pValMax = screenResult$pValMax),
    class="clumpingResult")
  return(result)
}
