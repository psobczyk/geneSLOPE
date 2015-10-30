#' identify_clump
#'
#' Enable interactive selection of snps in plot.
#' Return clump number.
#'
#' @title identify_clump
#' @param x appropiate class object
#' @param ... other arguments
#'
#' @rdname identify_clump
#' @export identify_clump
identify_clump <- function(x, ...){
  UseMethod("identify_clump")
}


#' Identify clump number in clumpingResult class plot
#'
#' @param x clumpingResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#'
#' @export
#'
identify_clump.clumpingResult <- function(x, ...) {
  plot.data <- NULL
  for(i in 1L:length(x$SNPclumps)){
    plot.data <- rbind(plot.data,
                       cbind(as.numeric(x$X_info[x$selectedSnpsNumbersScreening[x$SNPclumps[[i]]],1]),
                             as.numeric(x$X_info[x$selectedSnpsNumbersScreening[x$SNPclumps[[i]]],3]),
                             i, -log(x$pVals[x$selectedSnpsNumbersScreening[x$SNPclumps[[i]]]])))
  }
  rownames(plot.data) <- NULL
  plot.data <- data.frame(plot.data)
  colnames(plot.data) <- c("chromosome", "snp", "clump", "val")
  plot.data <- cbind(plot.data,
                     representatives = unlist(x$SNPclumps) %in% unlist(x$SNPnumber))
  granice <- aggregate(x$X_info[,3], list(x$X_info[,1]), max)
  granice_max <- cumsum(granice$x)
  granice$x <- c(0,head(cumsum(granice$x),-1))
  for(i in unique(plot.data$chromosome)){
    plot.data$snp[plot.data$chromosome==i] <- granice$x[i] +
      plot.data$snp[plot.data$chromosome==i]
  }

  a<- plot.data$snp
  b <- plot.data$val
  viewport_name <- unname(grid.ls(print = FALSE)[[1]][[4]])
  downViewport(viewport_name)
  # showViewport()
  # popViewport()
  pushViewport(viewport(xscale=c(0, max(granice_max)+1), yscale=c(0, 1.1*max(plot.data$val))))
  tmp = grid.locator()
  tmp.n <- as.numeric(tmp)
  diff.a <- (a-tmp.n[1])^2
  diff.b <- (b-tmp.n[2])^2
  paste("Selected SNP is in clump",
        plot.data$clump[which.min(diff.a/max(diff.a) + diff.b/max(diff.b))])
}

#' Identify clump number in geneSlopeResult class plot
#'
#' @param x geneSlopeResult class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#'
#' @export
identify_clump.geneSlopeResult <- function(x, ...) {
  plot.data <- NULL
  for(i in 1L:length(x$selectedClumps)){
    plot.data <- rbind(plot.data,
                       cbind(as.numeric(x$X_info[x$selectedSnpsClumpingNumbers[x$selectedClumps[[i]]],1]),
                             as.numeric(x$X_info[x$selectedSnpsClumpingNumbers[x$selectedClumps[[i]]],3]),
                             i, (x$effects[i]^2/var(x$y))))
  }
  rownames(plot.data) <- NULL
  plot.data <- data.frame(plot.data)
  colnames(plot.data) <- c("chromosome", "snp", "clump", "val")
  plot.data <- cbind(plot.data,
                     representatives = unlist(x$selectedClumps) %in% unlist(x$selectedSNPs))
  granice <- aggregate(x$X_info[,3], list(x$X_info[,1]), max)
  granice_max <- cumsum(granice$x)
  granice$x <- c(0,head(cumsum(granice$x),-1))
  for(i in unique(plot.data$chromosome)){
    plot.data$snp[plot.data$chromosome==i] <- granice$x[i] +
      plot.data$snp[plot.data$chromosome==i]
  }
  plot.data$val[plot.data$representatives] <- (x$effects^2/var(x$y))

  a<- plot.data$snp
  b <- plot.data$val
  viewport_name <- unname(grid.ls(print = FALSE)[[1]][[4]])
  downViewport(viewport_name)
  # showViewport()
  # popViewport()
  pushViewport(viewport(xscale=c(0, max(granice_max)+1), yscale=c(0, 1.1*max(plot.data$val))))
  tmp = grid.locator()
  tmp.n <- as.numeric(tmp)
  paste("Selected SNP is in clump",
        plot.data$clump[which.min((a-tmp.n[1])^2 + (b-tmp.n[2])^2)])
}
