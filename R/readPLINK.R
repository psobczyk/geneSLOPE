#' Reading PLINK .raw files
#'
#' Reading PLINK files that were previously suitably tranformed
#' in plink - see details
#'
#' @param file name of .raw file
#'
#' @return list \itemize{
#' \item snps matrix of SNPs
#' \item phenotype vector with phenotype values
#' \item positions vector with information about SNP positions in chromosome
#' }
#'
#' @details \strong{Exporting data from PLINK}
#' As stated in documentation for function \code{\link[adegenet]{read.PLINK}},
#' to import data to R, it needs to be exported from
#' PLINK using the option "--recodeA" (and NOT "--recodeAD").
#' The PLINK command should therefore look like:
#' \code{plink --file data --recodeA}.
#' For more information on this topic, please look at this webpage:
#' \url{http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml}
#'
readPLINK <- function(file){
  splittedName <- strsplit(file, "\\.")[[1]]
  if(tail(splittedName,1)!="raw")
    error(paste0("Error. ", file, " - .raw file is required."))
  file = splittedName[1]
  rawFile <- paste0(file, ".raw")
  if(!file.exists(rawFile))
    error("No .raw file. Cannot proceed")
  mapFile <- paste0(file, ".map")
  if(file.exists(mapFile)){
    snpData <- read.PLINK(file = rawFile, map.file = mapFile)
  } else {
    snpData <- read.PLINK(rawFile)
  }

  # we extract information about snp positions, phenotype and snps
  positions <- snpData$other$position
  phenotype <- snpData$other$phenotype
  snps <- sapply(snpData$gen, function(x) as.integer(x))
  result <- list(
    snps = snps,
    phenotype = phenotype,
    positions = positions
  )
  return(result)
}


#' Reading SNPs from one PLINK .raw files
#'
#' Reading PLINK files that were previously suitably tranformed
#' in plink - see details
#'
#' @param file name of .raw file
#'
#' @return list \itemize{
#' \item snps matrix of SNPs
#' \item observationInfo matrix containing non-genetic information including
#'  sex and phenotype
#' \item snpInfo matrix with information about SNPs names and positions in chromosomes
#' }
#'
#' @details \strong{Exporting data from PLINK}
#' To import data to R, it needs to be exported from
#' PLINK using the option "--recodeA" (and NOT "--recodeAD").
#' The PLINK command should therefore look like
#' \code{plink --file data --recodeA}.
#' For more information, please refer to:
#' \url{http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml}
#'
readPLINK2 <- function(file){
  # splittedName <- strsplit(file, "\\.")[[1]]
#   if(tail(splittedName,1)!="raw")
#     stop(paste0("Error. ", file, " - .raw file is required."))
#   file = splittedName[1]
#   rawFile <- paste0(file, ".raw")
#   if(!file.exists(rawFile))
#     stop("No .raw file. Cannot proceed")
  data_mat <- fread(file, stringsAsFactors = FALSE, header = TRUE)
  info <- data_mat[,1:6,with=FALSE]
  data_mat <- data_mat[, !(1:6), with=FALSE]
  splittedName <- strsplit(file, "\\.")[[1]]
  mapFile <- paste0(splittedName[1], ".map")
  if(file.exists(mapFile)){
    map_info <- fread(mapFile, stringsAsFactors = FALSE, header = FALSE)
    names(map_info) <- c("ChromosomeCode", "VariantIdentifier", "Position",
                         "BasePairCoordinate")
  } else {
    map_info = NULL
    message("No .map file found")
  }

  result <- list(
    snps = data_mat,
    observationInfo = info,
    snpInfo = map_info
  )
  return(result)
}
