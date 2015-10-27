#' Reading phenotype
#'
#' Reading phenotype info. It is assumed, that data comes in .fam file
#' First column is family id (FID), second is individual id (IID),
#' third is Paternal individual ID (PAT),
#' fourth is  Maternal individual ID (MAT), fifth is SEX
#' and six and last is PHENOTYPE
#' If file has only four columns, then it is assumed that PAT and MAT are missing.
#' If there is only one column, then it is assumed that only phenotype is provided.
#'
#' @export
#' @param filename name of file with phenotype
#' @param sep field seperator in file
#' @param header does first row of file contain variables names
#' @param stringAsFactors should character vectors be converted to factors?
#'
#' @return object of class phenotypeData
#'
readPhenotype <- function(filename, sep=" ", header=FALSE, stringAsFactors=FALSE){
  phe.data <- read.table(file = filename, header = header, sep = sep,
                         stringsAsFactors = stringAsFactors)
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
    yInfo = y_info),
    class="phenotypeData")
  return(result)
}
