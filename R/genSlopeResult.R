#' genSlopeResult class
#'
#' A result of procedure for snp clumping produced by \code{\link{genSLOPE}}
#'
#' @details Always a named list of ten elements
#' \enumerate{
#' \item \code{X} matrix of snps. One snp representative per each clump
#' \item \code{effects} coefficients in linear model with SNPs selected by SLOPE
#' \item \code{R2} R-squared in linear model with SNPs selected by SLOPE
#' \item \code{selectedSNPs} SNPs selected by SLOPE
#' \item \code{y} selectedClumps clumps that contain SNPs selected by SLOPE
#' \item \code{lambda} lambda values used by SLOPE procedure
#' \item \code{y} phenotype
#' \item \code{clumpRepresentatives} vector of column number in SNP matrix \code{X_all} of
#' clumps representative
#' \item \code{clumps} list of vectors of column numbers in SNP matrix \code{X_all}
#' for each clump member
#' \item \code{X_info} SNP info from .map file. Obtained from screening procedure.
#' \item \code{selectedSnpsNumbers} Contains information on which rows of \code{X_info}
#' matrix refer to selected clump representatives
#' \item \code{X_clumps} SNP matrix, output of clumping procedure
#' \item \code{X_all} SNP matrix,  input for clumping procedure
#' \item \code{numberOfSnps} Number of SNPs in input for screening procedure
#' \item \code{selectedSnpsNumbersScreening} Contains information on which rows of
#' \code{X_info} matrix refer to SNPs that are input in clumping procedure
#' \item \code{pValMax} p-value used in screening procedure
#' }
#' @seealso \code{\link{screeningResult}} \code{\link{clumpingResult}}
#' @name genSlopeResult
NULL


#' Print genSlopeResult class object
#'
#' @param x genSlopeResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @return Nothing.
#' @export
#'
#' @method print genSlopeResult
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
#' @param object genSlopeResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method summary genSlopeResult
summary.genSlopeResult <- function(object, ...){
  cat("Object of class genSlopeResult\n")
  cat(length(object$selectedSNPs), " snps selected\n")
  cat("R2 of fitted model ", object$R2)
}


#' Plot genSlopeResult class object
#'
#' @param x genSlopeResult class object
#' @param chromosomeNumber optional parameter, only selected chromosome will be plotted
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
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
