#' Reading phenotype
#'
#' Reading phenotype info. It is assumed that first column
#' is family id (FID), second is individual id (IID), \
#' third is Paternal individual ID (PAT),
#' fourth is  Maternal individual ID (MAT), fifth is SEX
#' and six and last is PHENOTYPE
#' If file has only for columns, then it is assumed that PAT and MAT are missing.
#'
#' @export
#' @param filename name of file with phenotype
#'
#' @return list \itemize{
#' \item y phenotype
#' \item y_info other observation info
#' }
#'
readPhenotype <- function(filename, sep=",", header=FALSE, stringAsFactors=FALSE){
  phe.data <- data.table::fread(input = filename, header = header, sep = sep,
                    stringsAsFactors = stringAsFactors)
  phe.data <- as.data.frame(phe.data)
  if(ncol(phe.data)==6){
    y = phe.data[,6]
    y_info = phe.data[,c(1,2,5)]
  } else if(ncol(phe.data)==4){
    y = phe.data[,4]
    y_info = phe.data[,c(1,2,3)]
  } else if(ncol(phe.data)==1){
    y = phe.data[,1]
    y_info = NULL
  } else{
    stop("Incorrect data")
  }

  #returning phenotype data
  result <- structure( list(
    y = y,
    y_info = y_info),
    class="phenotypeData")
  return(result)
}
