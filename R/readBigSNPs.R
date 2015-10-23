#' Reading SNPs from one PLINK .raw and .map files
#'
#' Reading PLINK files that were previously suitably tranformed
#' in plink - see details
#'
#' @export
#' @param rawFile name of .raw file
#' @param mapFile name of .map file
#' @param y phenotype
#' @param pValMax p-value threshold value used for screening
#' @param chunk.size size of chunk. The bigger the chunk the faster function works but
#' computer might run out of RAM
#' @param verbose if TRUE (default) information about progress is printed
#' @return object of class \code{\link{screeningResult}}
#'
#' @details \strong{Exporting data from PLINK}
#' To import data to R, it needs to be exported from
#' PLINK using the option "--recodeAD"
#' The PLINK command should therefore look like
#' \code{plink --file input --recodeAD --out output}.
#' For more information, please refer to:
#' \url{http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml}
#'
readBigSNPs <- function(rawFile, mapFile="", y, pValMax=0.05, chunk_size=1e3,
                        verbose = TRUE){
  if(file.exists(mapFile)){
    x_info <- read.table(mapFile)
  } else {
    x_info <- NULL
    message(".map file not found")
  }

  x <- read.big.matrix(filename = rawFile, sep=" ", header = TRUE,
                       type='double', shared=FALSE)
  numberOfSnps <- ncol(x) - 6
  if(verbose) message("Data sucesfully read")

  # see if dimensions match
  if(length(y) != nrow(x))
    stop("Length of phenotype y must match number of rows in matrix of SNPs")


  if(!is.null(x_info) & (nrow(x_info) != numberOfSnps))
    stop("Number of SNPs from .map and .raw data do not match")


  means <- colmean(x, na.rm=TRUE)
  temp <- NULL
  for(i in 1:nrow(x)){
    temp <- which(is.na(x[i,]))
    x[i,temp] <- means[temp]
  }
  if(verbose) message("Missing values imputation completed")

  suma = sum((y-mean(y))^2)
  n = length(y) - 2
  selectedSNPs <- NULL
  chunk=1
  p <- NULL
  x2 <- matrix(0, nrow=nrow(x), ncol=chunk_size)
  total <- floor(numberOfSnps/chunk_size)
  if(verbose){
    message("SNPs screening:")
    pb <- txtProgressBar(min = 0, max = total, style = 3)
  }
  for(i in 1:total){
    x2 <- x[,(7+(chunk-1)*chunk_size):(chunk*chunk_size+6)]
    p <- c(p, apply(x2, 2, function(snp){
      cps:::pValComp(snp, y, n, suma)
    }))
    rm(x2)
    if(verbose) setTxtProgressBar(pb, i)
    chunk = chunk + 1
  }
  if(verbose) close(pb)
  #last chunk
  if(ncol(x)>(7+(chunk-1)*chunk_size)){
    x2 <- x[,(7+(chunk-1)*chunk_size):ncol(x)]
    p <- c(p, apply(x2, 2, function(snp){
      cps:::pValComp(snp, y, n, suma)
    }))
    rm(x2)
  }

  x <- x[,6+which(p<pValMax)]

  result <- structure( list(
    X = x,
    y = y,
    X_info = x_info,#[p<pValMax,],
    numberOfSnps = numberOfSnps,
    selectedSnpsNumbers = which(p<pValMax),
    pValMax = pValMax),
    class="screeningResult")
  return(result)
}
