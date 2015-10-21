#' read data with bigmemory package
#'
#' @export
#'
readBigSNPs <- function(rawFile, mapFile, y, pValMax=0.05){
#   rawFile = "../wgas1.raw"
#   mapFile = "../wgas1.map"
#   y <- cps::y
#   pValMax=0.05
  y <- unname(y)
  suma = sum((y-mean(y))^2)
  n = length(y) - 2
  selectedSNPs <- NULL

  x_info <- read.table(mapFile)
  x_info <- x_info[rep(1:nrow(x_info), each=2),]
  x <- read.big.matrix(filename = rawFile, sep=" ", header = TRUE,
                       type='double', shared = TRUE)
  numberOfSnps <- ncol(x)
  x <- apply(x, 2, cps:::replace_na_with_mean)
  p <- apply(x, 2, function(snp){
    if(any(is.na(snp)))
       return(1)
    cps:::pValComp(snp, y, n, suma)
  })
  x <- x[,p<pValMax]

  result <- structure( list(
    X = x,
    y = y,
    X_info = x_info[p<pValMax,],
    numberOfSnps = numberOfSnps,
    selectedSnpsNumbers = which(p<pValMax),
    pValMax = pValMax),
    class="screeningResult")
  return(result)
}
