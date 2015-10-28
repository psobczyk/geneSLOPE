#' GWAS with SLOPE
#'
#' Performs GWAS with SLOPE on given snp matrix and phenotype.
#' At first clumping procedure is performed. Highly correlated
#' (that is stronger than parameter \emph{rho}) snps are clustered.
#' Then SLOPE is used on snp matrix which contains
#' one representative for each clump.
#'
#' @export
#' @param clumpingResult clumpProcedure output
#' @param fdr, False Discovery Rate for SLOPE
#' @param lambda lambda for SLOPE. See \code{\link[SLOPE]{create_lambda}}
#' @param verbose if TRUE progress bar is printed
#' @return object of class \code{\link{genSlopeResult}}
#'
#' @examples
#' \dontrun{
#' slope.result <- genSLOPE(clumping, fdr=0.1)
#' }
genSLOPE <- function(clumpingResult, fdr = 0.1, lambda="gaussian", verbose = TRUE){
  if(fdr>=1 | fdr <= 0){
    stop("FDR has to be within range (0,1)")
  }
  if(length(clumpingResult$y) != nrow(clumpingResult$X)){
    stop("Length of y must match
         number of rows in X")
  }

  lambda <- SLOPE::create_lambda(length(clumpingResult$y),
                                 clumpingResult$numberOfSnps, fdr, "gaussian")
  lambda <- lambda[1:ncol(clumpingResult$X)]
  slopeResult <- SLOPE::SLOPE(X = clumpingResult$X, y = clumpingResult$y,
                              fdr = fdr, lambda = lambda)

  selectedSNPs <- unlist(clumpingResult$SNPnumber)[slopeResult$selected]
  selectedSNPs <- sort(selectedSNPs)

  X_selected <- clumpingResult$X_all[,selectedSNPs]
  if(length(selectedSNPs)==0)
    X_selected <- rep(1, length(clumpingResult$y))
  # refitting linear model
  lm.fit.summary <- summary(lm(clumpingResult$y~X_selected))


  result <- structure(
    list( X = X_selected,
          effects = lm.fit.summary$coefficients[-1,1],
          R2 = lm.fit.summary$r.squared,
          selectedSNPs = selectedSNPs,
          selectedClumps = clumpingResult$SNPclumps[slopeResult$selected],
          lambda = lambda,
          y = clumpingResult$y,
          clumpRepresentatives = clumpingResult$SNPnumber,
          clumps = clumpingResult$SNPclumps,
          X_info = clumpingResult$X_info,
          X_clumps = clumpingResult$X,
          X_all = clumpingResult$X_all,
          selectedSnpsNumbers = clumpingResult$selectedSnpsNumbersScreening[selectedSNPs],
          selectedSnpsClumpingNumbers = clumpingResult$selectedSnpsNumbersScreening,
          numberOfSnps = clumpingResult$numberOfSnps,
          pValMax = clumpingResult$pValMax),
    class="genSlopeResult")
  return(result)
}

