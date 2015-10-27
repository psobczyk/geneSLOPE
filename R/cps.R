#' Genome Wide Association Study with SLOPE
#'
#' Genome-wide association study performed with SLOPE.
#' In first step, SNPs are clumped according to their correlations and
#' distances. Then, SLOPE is performed on data where each clump has
#' one representative.
#'
#' @docType package
#' @name cps
#' @details Version: 0.30.2
#' @importFrom SLOPE SLOPE
#' @importFrom SLOPE create_lambda
#' @import ggplot2
#' @importFrom bigmemory read.big.matrix
#' @author Piotr Sobczyk
#'
#' Maintainer: Piotr Sobczyk \email{Piotr.Sobczyk@@pwr.edu.pl}
#' @examples
#' \donttest{
#' famFile <- system.file("extdata", "plinkPhenotypeExample.fam", package = "cps")
#' mapFile <- system.file("extdata", "plinkMapExample.map", package = "cps")
#' snpsFile <- system.file("extdata", "plinkDataExample.raw", package = "cps")
#' phe <- readPhenotype(filename = famFile, sep=";")
#' screening <- readSNPs(snpsFile, mapFile, phe, pValMax = 0.05, chunkSize = 1e2)
#' clumping <- clumpProcedure(screening, rho = 0.3, verbose = TRUE)
#' slope.result <- genSLOPE(clumping, fdr=0.1)
#' }
#' \dontrun{
#' runExample()
#' }
NULL
