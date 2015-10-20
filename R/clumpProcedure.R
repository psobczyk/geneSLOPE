#' Clumping procedure
#'
#' Clumping of SNPs. Can be used as a preprocesing tool for using SLOPE in GWAS
#'
#' @export
#' @param y numeric vector, phenotype
#' @param X numeric matrix, SNPs
#' @param rho numeric, minimal correlation between two SNPs to be
#' @param verbose logical, if TRUE (default) progress bar is shown
#' classified to the same clump
clumpProcedure <- function(y, X, rho = 0.3, verbose = TRUE){

  if(verbose){
    message("Clumping procedure has started. Depending on
            size of your data this may take several minutes.")
    total = sqrt(ncol(X))
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
      setTxtProgressBar(pb, sqrt(i))
  }
  if(verbose)
    close(pb)
  nullClumps <- sapply(representatives, is.null)
  representatives <- representatives[!nullClumps]
  clumps <- clumps[!nullClumps]
  result <- structure(
    list( X = X[,unlist(representatives)],
          SNPnumber = representatives,
          SNPclumps = clumps),
    class="clumpingResult")
  return(result)
}

#' Print clumpingResult class object
#'
#' @param x clumpingResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
print.clumpingResult <- function(x, ...){
  cat("Object of class clumpingResult\n")
  cat("$X: Matrix\n")
  cat("\t", nrow(x$X), " observations\n")
  cat("\t", ncol(x$X), " snps\n")
  cat("$X_info: information about SNPs")
  cat("$SNPnumber: list with snp representatives for clumps \n")
  cat("\t[", paste(head(x$SNPnumber), collapse=","), "..., ]\n")
  cat("$SNPclumps: list of vector, number of snps in clumps\n")
}

#' Summary clumpingResult class object
#'
#' @param x clumpingResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
summary.clumpingResult <- function(x, ...){
  cat("Object of class clumpingResult\n")
  cat(length(x$SNPclumps), " clumps with ", ncol(x$X_all), " SNPs\n")
  cat("Average clumps size ", mean(unlist(lapply(x$SNPclumps, length))), "\n")
  cat("Smallest clump size ", min(unlist(lapply(x$SNPclumps, length))), "\n")
  cat("Biggest clump size ", max(unlist(lapply(x$SNPclumps, length))), "\n")
}


#' Clumping procedure
#'
#' Clumping of SNPs previously read and screened by \code{\link{readSNPs}}
#'
#' @export
#' @param screenResult object of class screenResult
#' @param rho numeric, minimal correlation between two SNPs to be
#' @param verbose logical, if TRUE (default) progress bar is shown
#' classified to the same clump
clumpProcedure2 <- function(screenResult, rho = 0.3, verbose = TRUE){

  if(verbose){
    message("Clumping procedure has started. Depending on
            size of your data this may take several minutes.")
    total = sqrt(ncol(X))
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
  pVals <- apply(X, 2, function(x) pValComp(x,screenResult$y,n,suma))

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
    list( X = X[,unlist(representatives)],
          X_all = X,
          y = screenResult$y,
          X_info = screenResult$X_info,
          SNPnumber = representatives,
          SNPclumps = clumps,
          selectedSnpsNumbers = screenResult$selectedSnpsNumbers[unlist(representatives)],
          numberOfSnps = screenResult$numberOfSnps,
          selectedSnpsNumbersScreening = screenResult$selectedSnpsNumbers,
          pValMax = screenResult$pValMax),
    class="clumpingResult")
  return(result)
}
