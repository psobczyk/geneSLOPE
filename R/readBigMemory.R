#' read data with bigmemory package
#'
#' @export
#'
readBigSNPs <- function(filename, y, pValMax=0.05, ncores=1){
  # library(bigmemory)
  # library(foreach)
#   library(doMC)
#   library(parallel)
#   registerDoMC(4)
#   filename = "../plink.tped"
#   y <- cps::y
#   pValMax=0.05
  doMC::registerDoMC(ncores)
  x <- read.big.matrix(filename = filename, sep=" ", shared = TRUE)
  y <- unname(y)
  suma = sum((y-mean(y))^2)
  n = length(y) - 2
  numberOfSnps = 2*nrow(x)
  selectedSNPs <- NULL
  x2 <- foreach(ind=1:nrow(x), .combine=rbind) %dopar% {
      temp = cps:::recodeAD(x[ind,-(1:4)])
      temp = apply(temp, 1, cps:::replace_na_with_mean)
      p <- apply(temp, 2, cps:::pValComp, y, n, suma)<pValMax
      if(p[1] & p[2]){
        selectedSNPs <- c(selectedSNPs, ind, ind)
        rbind(temp[,1], temp[,2])
      } else if (p[1] & !p[2]){
        selectedSNPs <- c(selectedSNPs, ind)
        temp[,1]
      } else if (!p[1] & p[2]){
        selectedSNPs <- c(selectedSNPs, ind)
        temp[,2]
      } else{
        NULL
      }
  }
  rm(x)
  result <- structure( list(
    X = x2,
    y = y,
    X_info = NULL,
    numberOfSnps = numberOfSnps,
    selectedSnpsNumbers = selectedSNPs,
    pValMax = pValMax),
    class="screeningResult")
  return(result)
}
