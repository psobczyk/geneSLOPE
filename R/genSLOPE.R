#' GWAS with SLOPE
#'
#' Performs GWAS with SLOPE on given snp matrix and phenotype.
#' At first clumping procedure is performed. Highly correlated
#' (that is stronger than parameter \emph{rho}) snps are clustered.
#' Then SLOPE is used on snp matrix which contains
#' one representative for each clump.
#'
#' @export
#' @param clumpProcedure
#' @param fdr, False Discovery Rate for SLOPE
#' @param lambda lambda for SLOPE. See \code{\link[SLOPE]{create_lambda}}
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
          numberOfSnps = clumpingResult$numberOfSnps,
          pValMax = clumpingResult$pValMax),
    class="genSlopeResult")
  return(result)
}

#' Print genSlopeResult class object
#'
#' @param x genSlopeResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
print.genSlopeResult <- function(x, ...){
  cat("Object of class genSlopeResult\n")
  cat("$X: Matrix\n")
  cat("\t", nrow(x$X), " observations\n")
  cat("\t", ncol(x$X), " snps\n")
  cat("$effects: effect size for selected snps\n")
  cat("$R2: R squared for fitted lm model\n")
  cat("$selectedSNPs: SNPs selected by SLOPE\n")
  cat("$selectedClumps: clumps that contain SNPs selected by SLOPE\n")
  cat("$lambda: lambdas used by SLOPE\n")
  cat("$y: Phenotype\n")
  cat("$X_clump: Matrix after clumping\n")
  cat("$X_all: Matrix before clumping\n")
  cat("$X_info: Information about snps\n")
  cat("$clumpRepresentatives: clumps representatives\n")
  cat("$clumps: list of clumps\n")
  cat("$selectedSnpsNumbers: snps selected by SLOPE\n")
  cat("$numberOfSnps: Number of SNPs before screening\n")
  cat("$pValMax: p-value threshold\n")
}

#' Summary genSlopeResult class object
#'
#' @param x genSlopeResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
summary.genSlopeResult <- function(x, ...){
  cat("Object of class genSlopeResult\n")
  cat(length(x$selectedSNPs), " snps selected\n")
  cat("R2 of fitted model ", x$R2)
}


#' Plot genSlopeResult class object
#'
#' @param x genSlopeResult class object
#' @param chromosome optional parameter, only selected chromosome will be plotted
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
plot.genSlopeResult <- function(x, chromosomeNumber=NULL, ...){
  plot.data <- data.frame(cbind(chromosome=as.numeric(x$X_info[unlist(x$selectedClumps),1]),
                                snp=as.numeric(x$X_info[unlist(x$selectedClumps),3]),
                                val=2.0))
  granice <- aggregate(plot.data$snp, list(plot.data$chromosome), max)
  granice$x <- c(0,head(cumsum(granice$x),-1))
  for(i in unique(plot.data$chromosome)){
    plot.data$snp[plot.data$chromosome==i] <- granice$x[which(granice$Group.1==i)] +
      plot.data$snp[plot.data$chromosome==i]
  }
  representatives = which(unlist(x$selectedClumps) %in% unlist(x$selectedSNPs))
  plot.data$val[representatives] <- 4
  if(!is.null(chromosomeNumber))
    plot.data <- subset(plot.data, chromosome==chromosomeNumber)
  plot.data <- plot.data[order(plot.data$chromosome, plot.data$snp),]
  representatives <- which(plot.data$val==4)
  ggplot(plot.data) + geom_point(aes(x=snp, y=val, colour = "red", size = 6),
                                 plot.data[representatives,]) +
    geom_segment(aes(x=snp, xend=snp, y=0, yend=val, alpha=val/4)) +
    ylab("") + scale_y_continuous(breaks=NULL) +
    xlab("Genome") + scale_x_continuous(breaks=granice$x, labels=granice$Group.1) +
    scale_alpha_continuous(guide=FALSE) +
    scale_color_discrete(guide=FALSE) +
    scale_size_area(guide=FALSE) +
    theme(panel.background=element_blank())
}
