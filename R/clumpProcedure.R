#' Clumping procedure
#'
#' Clumping of SNPs previously read and screened by \code{\link{readBigSNPs}}
#'
#' @export
#' @param screenResult object of class screenResult
#' @param rho numeric, minimal correlation between two SNPs to be
#' @param verbose logical, if TRUE (default) progress bar is shown
#' classified to the same clump
#'
#' @return object of class \code{\link{clumpingResult}}
#'
clumpProcedure <- function(screenResult, rho = 0.3, verbose = TRUE){

  if(verbose){
    message("Clumping procedure has started. Depending on
            size of your data this may take several minutes.")
    total = sqrt(ncol(screenResult$X))
    # create progress bar
    pb <- txtProgressBar(min = 0, max = total, style = 3)
  }

  if(rho>=1 | rho <= 0){
    stop("Rho has to be within range (0,1)")
  }
  if(length(screenResult$y) != nrow(screenResult$X)){
    stop("Length of y must match
         number of rows in X")
  }
  suma = sum((screenResult$y-mean(screenResult$y))^2)
  n = length(screenResult$y) - 2
  pVals <- apply(screenResult$X, 2, function(x) pValComp(x,screenResult$y,n,suma))

  a <- order(pVals, decreasing = FALSE)
  notClumped <- rep(TRUE, length(a))
  clumps <- list()
  representatives <- list()

  i <- 1
  while(any(notClumped)){
    idx = a[i]
    if(notClumped[idx]){
      clump <- abs(apply(screenResult$X[,notClumped, drop=FALSE], 2, cor, screenResult$X[,idx]))>rho
      clumps[[i]] <- which(notClumped)[clump]
      representatives[[i]] <- idx
      notClumped[ which(notClumped)[clump] ] <- FALSE
    }
    i = i+1
    if(verbose)
      setTxtProgressBar(pb, sqrt(i))
  }
  if(verbose)
    close(pb)
  nullClumps <- sapply(representatives, is.null)
  representatives <- representatives[!nullClumps]
  clumps <- clumps[!nullClumps]
  result <- structure(
    list( X = screenResult$X[,unlist(representatives)],
          y = screenResult$y,
          SNPnumber = representatives,
          SNPclumps = clumps,
          X_info = screenResult$X_info,
          selectedSnpsNumbers = screenResult$selectedSnpsNumbers[unlist(representatives)],
          X_all = screenResult$X,
          numberOfSnps = screenResult$numberOfSnps,
          selectedSnpsNumbersScreening = screenResult$selectedSnpsNumbers,
          pValMax = screenResult$pValMax),
    class="clumpingResult")
  return(result)
}


#' Print clumpingResult class object
#'
#' @param x clumpingResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
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
#' @param x clumpingResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#' @keywords internal
summary.clumpingResult <- function(x, ...){
  cat("Object of class clumpingResult\n")
  cat(length(x$SNPclumps), " clumps with ", ncol(x$X_all), " SNPs\n")
  cat("Average clumps size ", mean(unlist(lapply(x$SNPclumps, length))), "\n")
  cat("Smallest clump size ", min(unlist(lapply(x$SNPclumps, length))), "\n")
  cat("Biggest clump size ", max(unlist(lapply(x$SNPclumps, length))), "\n")
}


#' Plot clumpingResult class object
#'
#' @param x clumpingResult class object
#' @param chromosome optional parameter, only selected chromosome will be plotted
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
