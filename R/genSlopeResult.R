#' genSlopeResult class
#'
#' A result of applying SLOPE to matrix of SNPs obtained by
#' clumping produced. Result of function \code{\link{genSLOPE}}
#'
#' @details Always a named list of ten elements
#' \enumerate{
#' \item \code{X} numeric matrix, consists of one snp representative for each clump
#' selected by SLOPE
#' \item \code{effects} numeric vector, coefficients in linear model build on
#' snps selected by SLOPE
#' \item \code{R2} numeric, value of R-squared in linear model build on
#' snps selected by SLOPE
#' \item \code{selectedSNPs} which columns in matrix \code{X_all}
#' are related to snps selected by SLOPE
#' \item \code{y} selectedClumps list of numeric vectors, which columns in SNP matrix
#' \code{X_all} are related to clump members selected by SLOPE
#' \item \code{lambda} numeric vector, lambda values used by SLOPE procedure
#' \item \code{y} numeric vector, phenotype
#' \item \code{clumpRepresentatives} numeric vector, which columns in SNP matrix \code{X_all}
#' are related to clumps representatives
#' \item \code{clumps} list of numeric vectors, which columns in SNP matrix
#' \code{X_all} are related to clump members
#' \item \code{X_info} data.frame, mapping information about SNPs from .map file.
#' Copied from the result of clumping procedure
#' \item \code{selectedSnpsNumbers} numeric vector, which rows of \code{X_info}
#' data.frame are related to snps that were selected by SLOPE
#' \item \code{X_clumps} numeric matrix, consists of one snp representative for each clump
#' \item \code{X_all} numeric matrix, all the snps that passed screening procedure
#' \item \code{selectedSnpsNumbersScreening} numeric vector, which rows of \code{X_info}
#' data.frame are related to snps that passed screening
#' \item \code{numberOfSnps} numeric, total number of SNPs before screening procedure
#' \item \code{pValMax} numeric, p-value used in screening procedure
#' }
#' @seealso \code{\link{screeningResult}} \code{\link{clumpingResult}}
#' \code{\link{genSLOPE}} \code{\link{SLOPE::SLOPE}}
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
  cat("$X: numeric matrix\n")
  cat("\t", nrow(x$X), " rows\n")
  cat("\t", ncol(x$X), " columns\n")
  cat("$effects: numeric vector of length ", length(x$effects), "\n")
  cat("$R2: ", x$R2, "\n")
  cat("$selectedSNPs: numeric vector of length",
      length(x$selectedSNPs), "\n")
  cat("$selectedClumps: list of vectors of length",
      length(x$selectedClumps), "\n")
  cat("$lambda: numeric vector of length",
      length(x$lambda), "\n")
  cat("$y: numeric vector\n")
  cat("$X_clump: Matrix after clumping\n")
  cat("\t", nrow(x$X_clump), " rows\n")
  cat("\t", ncol(x$X_clump), " columns\n")
  cat("$X_all: Matrix before clumping\n")
  cat("\t", nrow(x$X_all), " rows\n")
  cat("\t", ncol(x$X_all), " columns\n")
  cat("$X_info: Information about snps\n")
  cat("\t", nrow(x$X_info), " rows\n")
  cat("\t", ncol(x$X_info), " columns\n")
  cat("$clumpRepresentatives: numeric vector of length",
      length(x$clumpRepresentatives), "\n")
  cat("$clumps: list of numeric vectors of length",
      length(x$clumpRepresentatives), "\n")
  cat("$selectedSnpsNumbers: numeric vector of length",
      length(x$clumps), "\n")
  cat("$selectedSnpsClumpingNumbers: numeric vector of length",
      length(x$selectedSnpsClumpingNumbers), "\n")
  cat("$numberOfSnps: number of SNPs before screening:", x$numberOfSnps, "\n")
  cat("$pValMax: p-value threshold: ", x$pValMax, "\n")
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
  cat("Effect size for selected snps\n")
  cat(object$effects, "\n")
  cat("R square of fitted model ", object$R2)
}


#' Plot genSlopeResult class object
#'
#' @param x genSlopeResult class object
#' @param chromosomeNumber optional parameter, only selected chromosome will be plotted
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
plot.genSlopeResult <- function(x, chromosomeNumber=NULL, ...){
  if(!is.null(x$X_info)){
    plot.data <- NULL
    for(i in 1L:length(x$selectedClumps)){
      plot.data <- rbind(plot.data,
                         cbind(as.numeric(x$X_info[x$selectedSnpsClumpingNumbers[x$selectedClumps[[i]]],1]),
                               as.numeric(x$X_info[x$selectedSnpsClumpingNumbers[x$selectedClumps[[i]]],3]),
                               i, abs(x$effects[i])/2))
    }
    rownames(plot.data) <- NULL
    plot.data <- data.frame(plot.data)
    colnames(plot.data) <- c("chromosome", "snp", "clump", "val")
    granice <- aggregate(plot.data$snp, list(plot.data$chromosome), max)
    granice$x <- c(0,head(cumsum(granice$x),-1))
    for(i in unique(plot.data$chromosome)){
      plot.data$snp[plot.data$chromosome==i] <- granice$x[which(granice$Group.1==i)] +
        plot.data$snp[plot.data$chromosome==i]
    }
    representatives = which(unlist(x$selectedClumps) %in% unlist(x$selectedSNPs))
    plot.data$val[representatives] <- abs(x$effects)
    if(!is.null(chromosomeNumber))
      plot.data <- subset(plot.data, chromosome==chromosomeNumber)
    plot.data$clump <- as.factor(plot.data$clump)
    ggplot(plot.data) + geom_point(aes(x=snp, y=val, colour = clump, size = 6),
                                   plot.data[representatives,]) +
      geom_segment(aes(x=snp, xend=snp, y=0, yend=val, alpha=val/4, color=clump)) +
      ylab("Effect size") + scale_y_continuous() +
      xlab("Genome") + scale_x_continuous(breaks=granice$x, labels=granice$Group.1) +
      scale_alpha_continuous(guide=FALSE) +
      scale_color_discrete("Clump") +
      scale_size_area(guide=FALSE) +
      theme(panel.background=element_blank(),
            panel.grid.major.y=element_line(colour = "grey80"),
            panel.grid.minor.y=element_line(colour = "grey90"),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank())
    } else {
      plot.data <- NULL
      for(i in 1L:length(x$selectedClumps)){
        plot.data <- rbind(plot.data,
                           cbind(x$selectedSnpsClumpingNumbers[unlist(x$selectedClumps[[i]])],
                                 i, abs(x$effects[i])/2))
      }
      plot.data <- data.frame(plot.data)
      colnames(plot.data) <- c("snp", "clump", "val")
      rownames(plot.data) <- NULL
      ggplot(plot.data) + geom_point(aes(x=snp, y=val, colour = "red", size = 6),
                                     plot.data[representatives,]) +
        geom_segment(aes(x=snp, xend=snp, y=0, yend=val, alpha=val/4)) +
        ylab("") + scale_y_continuous() +
        xlab("SNP number") +
        scale_alpha_continuous(guide=FALSE) +
        scale_color_discrete(guide=FALSE) +
        scale_size_area(guide=FALSE) +
        theme(panel.background=element_blank(),
              panel.grid.major.y=element_line(colour = "grey80"),
              panel.grid.minor.y=element_line(colour = "grey90"),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid.major.x=element_blank(),
              panel.grid.minor.x=element_blank())
    }
}
