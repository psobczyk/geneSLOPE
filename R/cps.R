#' Genome Wide Association Study with SLOPE
#'
#' Genome-wide association study performed with SLOPE.
#' In first step, SNPs are clumped according to their correlations and
#' distances. Then, SLOPE is performed on data where each clump has
#' one representative.
#'
#' @docType package
#' @name cps
#' @details Version: 0.29.7
#' @importFrom data.table fread
#' @importFrom SLOPE SLOPE
#' @importFrom SLOPE create_lambda
#' @import ggplot2
#' @import bigmemory
#' @import biganalytics
#' @author Piotr Sobczyk
#'
#' Maintainer: Piotr Sobczyk \email{Piotr.Sobczyk@@pwr.edu.pl}
#' @examples
#' \donttest{
#' famFile <- system.file("extdata", "plinkPhenotypeExample.fam", package = "cps")
#' mapFile <- system.file("extdata", "plinkMapExample.map", package = "cps")
#' snpsFile <- system.file("extdata", "plinkDataExample.raw", package = "cps")
#' phe <- readPhenotype(filename = famFile, sep=";")
#' screening <- readBigSNPs(snpsFile, mapFile, phe$y, pValMax = 0.05, chunk_size = 1e2)
#' clumping <- clumpProcedure2(screening, rho = 0.3, verbose = FALSE)
#' slope.result <- genSLOPE(clumping, fdr=0.1)
#' }
#' \dontrun{
#' runExample()
#' }
NULL
