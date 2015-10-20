#' Reading SNPs from one PLINK .tped and .tfam files
#'
#' Reading PLINK files that were previously suitably tranformed
#' in plink - see details
#'
#' @export
#' @param file name of .tped file
#' @param y phenotype
#' @param pValMax p-value threshold value used for screening
#' @param chunk.size size of chunk. The bigger the chunk the faster function works but
#' computer might run out of RAM
#' @param verbose should info be present
#' @return list \itemize{
#' \item X matrix of filtered SNPs
#' \item y phenotype
#' \item X_info SNPs info
#' \item numberOfSnps number of snps read from .tped file
#' \item pValMax p-value threshold value used for screening
#' }
#'
#' @details \strong{Exporting data from PLINK}
#' To import data to R, it needs to be exported from
#' PLINK using the option "--recode12 --transposed"
#' The PLINK command should therefore look like
#' \code{plink --file data --recode12 --transposed --out temp}.
#' For more information, please refer to:
#' \url{http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml}
#'
readSNPs <- function(file, y, pValMax = 0.05, chunk.size = 5e4, verbose = TRUE){
  done = FALSE
  chunk = 1
  numberOfSnps <- 0
  selectedSnpsNumbers <- NULL
  X <- NULL
  X_info <- NULL
  y <- unname(y)
  suma = sum((y-mean(y))^2)
  n = length(y) - 2
  message("Reading and screening data. Depending on data size this might take
          several minutes")
  while(!done)
  {
    temp = data.table::fread(file,skip=(chunk-1)*chunk.size+1,nrow=chunk.size,
                             header = FALSE)
    temp <- as.data.frame(temp)
    SNPinfo <- t(temp[, (1:4)])
    SNPinfo <- SNPinfo[,rep(1:ncol(SNPinfo), each=2)]
    temp <- temp[, -(1:4)]
    temp <- apply(temp, 1, function(row) {
      a <- recodeAD(row)
      c(a[1,], a[2,])
    })
    temp2 <- matrix(nrow = nrow(temp)/2, ncol = 2*ncol(temp))
    temp2[,seq(1, 2*ncol(temp), 2)] <- temp[1:(nrow(temp)/2),]
    temp2[,seq(2, 2*ncol(temp), 2)] <- temp[(nrow(temp)/2+1):nrow(temp),]
    temp <- temp2
    rm("temp2")
    temp <- apply(temp, 2, replace_na_with_mean)
    pVals <- apply(temp, 2,
                   function(x) pValComp(x,y,n,suma))
    smallPVals <- pVals<pValMax
    X <- cbind(X, temp[,smallPVals])
    numberOfSnps <- numberOfSnps + ncol(temp)
    selectedSnpsNumbers <- c(selectedSnpsNumbers,
                             2*(chunk-1)*chunk.size+which(smallPVals))
    SNPinfo <- SNPinfo[,smallPVals]
    X_info <- rbind(X_info, t(SNPinfo))
    if(verbose) message(paste(chunk*chunk.size, "SNPs processed"))
    chunk = chunk + 1
    if(ncol(temp)/2<chunk.size) done = TRUE
  }
  #returning screening result
  result <- structure( list(
    X = X,
    y = y,
    X_info = X_info,
    numberOfSnps = numberOfSnps,
    selectedSnpsNumbers = selectedSnpsNumbers,
    pValMax = pValMax),
    class="screeningResult")
  return(result)
}

#' Print screeningResult class object
#'
#' @param x screeningResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
print.screeningResult <- function(x, ...){
  cat("Object of class screeningResult\n")
  cat("$X: Matrix of SNPs\n")
  cat("\t", nrow(x$X), " observations\n")
  cat("\t", ncol(x$X), " snps\n")
  cat("$y: a vector containing phenotype \n")
  cat("$X_info: : Matrix of SNPs info\n")
  cat("\t", nrow(x$X_info), " Info categories\n")
  cat("\t", ncol(x$X_info), " snps\n")
  cat("$numberOfSnps: Number of SNPs in .tped file\n")
  cat("$selectedSnpsNumbers: Column numbers of selected SNPs in original X matrix\n")
  cat("$pValMax: p-value threshold value used for screening\n")
}


#' Summary screeningResult class object
#'
#' @param x screeningResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
summary.screeningResult <- function(x, ...){
  cat("Object of class screeningResult\n")
  cat(nrow(x$X), " observations\n")
  cat(x$numberOfSnps, " SNPs were screened\n")
  cat(ncol(x$X), " snps had p-value under threshold in marginal test\n")
}
