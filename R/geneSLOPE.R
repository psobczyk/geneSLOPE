#' Genome Wide Association Study with SLOPE
#'
#' Package geneSLOPE performes Genome-wide association study with SLOPE.
#' This study is split into three steps.
#' \itemize{
#' \item In the first step data is read using \pkg{\link{bigmemory}} package and immediatly
#' screened using marginal tests for each SNP
#' \item SNPs are clumped using correlations
#' \item SLOPE is performed on data where each clump has
#' one representative (therefore we ensure that variables in linear model
#' are not strognly correlated)
#' }
#'
#'
#' @docType package
#' @name geneSLOPE
#' @details Version: 0.32.5
#' @importFrom SLOPE SLOPE
#' @importFrom SLOPE create_lambda
#' @import ggplot2
#' @importFrom bigmemory read.big.matrix
#' @author Malgorzata Bogdan, Damian Brzyski, Christine Peterson, Chiara Sabatti,
#' Piotr Sobczyk
#'
#' Maintainer: Piotr Sobczyk \email{Piotr.Sobczyk@@pwr.edu.pl}
#'
#' @references \emph{SLOPE -- Adaptive Variable Selection via Convex Optimization},
#' Malgorzata Bogdan, Ewout van den Berg, Chiara Sabatti,
#' Weijie Su and Emmanuel Candes
#'
#' @examples
#' \donttest{
#' famFile <- system.file("extdata", "plinkPhenotypeExample.fam", package = "cps")
#' mapFile <- system.file("extdata", "plinkMapExample.map", package = "cps")
#' snpsFile <- system.file("extdata", "plinkDataExample.raw", package = "cps")
#' phe <- readPhenotype(filename = famFile, sep=";")
#' screening.result <- readSNPs(snpsFile, mapFile, phe, pValMax = 0.05, chunkSize = 1e2)
#' clumping.result <- clumpingProcedure(screening.result, rho = 0.3, verbose = TRUE)
#' slope.result <- genSLOPE(clumping.result, fdr=0.1)
#' }
#' \dontrun{
#' gui_geneSLOPE()
#' }
NULL