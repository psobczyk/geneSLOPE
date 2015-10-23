#' clumpingResult class
#'
#' A result of procedure for snp clumping produced by \code{clumpProcedure}
#'
#' @details Always a named list of ten elements
#' \enumerate{
#' \item \code{X} matrix of snps. One snp representative per each clump
#' \item \code{y} phenotype
#' \item \code{SNPnumber} vector of column number in SNP matrix \code{X_all} of
#' clumps representative
#' \item \code{SNPclumps} list of vectors of column numbers in SNP matrix \code{X_all}
#' for each clump member
#' \item \code{X_info} SNP info from .map file. Obtained from screening procedure.
#' \item \code{selectedSnpsNumbers} Contains information on which rows of \code{X_info}
#' matrix refer to selected clump representatives
#' \item \code{X_all} X matrix that was input for clumping procedure
#' \item \code{numberOfSnps} Number of SNPs in input for screening procedure
#' \item \code{selectedSnpsNumbersScreening} Contains information on which rows of
#' \code{X_info} matrix refer to SNPs that are input in clumping procedure
#' \item \code{pValMax} p-value used in screening procedure
#' }
#' @seealso \code{\link{screeningResult}} \code{\link{clumpProcedure}}
#' @name clumpingResult
NULL

#' Print clumpingResult class object
#'
#' @param x clumpingResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method print clumpingResult
print.clumpingResult <- function(x, ...){
  cat("Object of class clumpingResult\n")
  cat("$X: Matrix\n")
  cat("\t", nrow(x$X), " observations\n")
  cat("\t", ncol(x$X), " snps\n")
  cat("$y: vector with phenotype")
  cat("$SNPnumber: list with snp representatives for clumps \n")
  cat("\t[", paste(head(x$SNPnumber), collapse=","), "..., ]\n")
  cat("$SNPclumps: list of vector, number of snps in clumps\n")
  cat("$X_info: information about SNPs\n")
  cat("$selectedSnpsNumbers: ")
}

#' Summary clumpingResult class object
#'
#' @param object clumpingResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method summary clumpingResult
summary.clumpingResult <- function(object, ...){
  cat("Object of class clumpingResult\n")
  cat(length(object$SNPclumps), " clumps with ", ncol(object$X_all), " SNPs\n")
  cat("Average clumps size ", mean(unlist(lapply(object$SNPclumps, length))), "\n")
  cat("Smallest clump size ", min(unlist(lapply(object$SNPclumps, length))), "\n")
  cat("Biggest clump size ", max(unlist(lapply(object$SNPclumps, length))), "\n")
}


#' Plot clumpingResult class object
#'
#' @param x clumpingResult class object
#' @param chromosomeNumber optional parameter, only selected chromosome will be plotted
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
plot.clumpingResult <- function(x, chromosomeNumber=NULL, ...){
  if(!is.null(x$X_info)){
    plot.data <- data.frame(cbind(chromosome=as.numeric(x$X_info[unlist(x$SNPclumps),1]),
                                  snp=as.numeric(x$X_info[unlist(x$SNPclumps),3]),
                                  val=2.0))
    granice <- aggregate(plot.data$snp, list(plot.data$chromosome), max)
    granice$x <- c(0,head(cumsum(granice$x),-1))
    for(i in unique(plot.data$chromosome)){
      plot.data$snp[plot.data$chromosome==i] <- granice$x[which(granice$Group.1==i)] +
        plot.data$snp[plot.data$chromosome==i]
    }
    representatives = which(unlist(x$SNPclumps) %in% unlist(x$SNPnumber))
    plot.data$val[representatives] <- 4
    if(!is.null(chromosomeNumber))
      plot.data <- subset(plot.data, plot.data$chromosome==chromosomeNumber)
    plot.data <- plot.data[order(plot.data$chromosome, plot.data$snp),]
    representatives <- which(plot.data$val==4)
    ggplot(plot.data) + geom_point(aes(x=snp, y=val, colour = "red", size = 6),
                                   data=plot.data[representatives,]) +
      geom_segment(aes(x=snp, xend=snp, y=0, yend=val, alpha=val/4)) +
      ylab("") + scale_y_continuous(breaks=NULL) +
      xlab("Genome") + scale_x_continuous(breaks=granice$x, labels=granice$Group.1) +
      scale_alpha_continuous(guide=FALSE) +
      scale_color_discrete(guide=FALSE) +
      scale_size_area(guide=FALSE) +
      theme(panel.background=element_blank())
  } else {
    plot.data <- data.frame(cbind(snp=x$selectedSnpsNumbersScreening,
                                  val=2.0))
    representatives = which(x$selectedSnpsNumbersScreening %in% x$selectedSnpsNumbers)
    plot.data$val[representatives] <- 4
    representatives <- which(plot.data$val==4)
    ggplot(plot.data) + geom_point(aes(x=snp, y=val, colour = "red", size = 6),
                                   plot.data[representatives,]) +
      geom_segment(aes(x=snp, xend=snp, y=0, yend=val, alpha=val/4)) +
      ylab("") + scale_y_continuous(breaks=NULL) +
      xlab("SNP number") +
      scale_alpha_continuous(guide=FALSE) +
      scale_color_discrete(guide=FALSE) +
      scale_size_area(guide=FALSE) +
      theme(panel.background=element_blank())
  }

}
