#' read data with bigmemory package
#'
#' @export
#'
readBigSNPs <- function(rawFile, mapFile, y, pValMax=0.05, chunk_size=1e3){
#   rawFile = "../wgas1.raw"
#   mapFile = "../wgas1.map"
#   y <- cps::y
#   pValMax=0.05
  y <- unname(y)
  suma = sum((y-mean(y))^2)
  n = length(y) - 2
  selectedSNPs <- NULL

  x_info <- read.table(mapFile)
  x <- read.big.matrix(filename = rawFile, sep=" ", header = TRUE,
                       type='double', shared=FALSE)
  numberOfSnps <- ncol(x)

  means <- colmean(x, na.rm=TRUE)
  temp <- NULL
  for(i in 1:nrow(x)){
    temp <- which(is.na(x[i,]))
    x[i,temp] <- means[temp]
  }

  d2 <- ncol(x) - 6
  chunk=1
  p <- NULL
  x2 <- matrix(0, nrow=nrow(x), ncol=chunk_size)
  for(i in 1:floor(d2/chunk_size)){
    x2 <- x[,(7+(chunk-1)*chunk_size):(chunk*chunk_size+6)]
    p <- c(p, apply(x2, 2, function(snp){
      cps:::pValComp(snp, y, n, suma)
    }))
    rm(x2)
    print(chunk)
    chunk = chunk + 1
  }

  x2 <- x[,(7+(chunk-1)*chunk_size):ncol(x)]
  p <- c(p, apply(x2, 2, function(snp){
    cps:::pValComp(snp, y, n, suma)
  }))
  rm(x2)

#   x <- as.matrix(x)
#   x2 <- Matrix(x, sparse = TRUE)
#
#   p <- apply(x2[,-c(1:6)], 2, function(snp){
#     if(any(is.na(temp)))
#        return(1)
#     cps:::pValComp(temp, y, n, suma)
#   })

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
