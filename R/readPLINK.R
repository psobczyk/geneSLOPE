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
