#' identify_clump
#'
#' Enable interactive selection of snps in plot.
#' Return clump number.
#'
#' @title identify_clump
#' @param x appropriate class object
#' @param ... other arguments
#' @return No return value, called for side effects
#'
#' @rdname identify_clump
#' @export identify_clump
identify_clump <- function(x, ...){
  UseMethod("identify_clump")
}


#' Identify clump number in \code{\link{clumpingResult}} class plot
#'
#' @param x \code{\link{clumpingResult}} class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @return No return value, called for side effects
#' @export
#'
identify_clump.clumpingResult <- function(x, ...) {
  plot.data <- create_clumping_plot_data(x)

  if(length(unique(x$X_info[,3])) == 1){
    chromosome_limits <- aggregate(x$X_info[,4], list(x$X_info[,1]), max)
  } else {
    chromosome_limits <- aggregate(x$X_info[,3], list(x$X_info[,1]), max)
  }

  chromosome_limits_max <- cumsum(chromosome_limits$x)
  chromosome_limits$x <- c(0,head(cumsum(chromosome_limits$x),-1))

  a <- plot.data$snp
  b <- plot.data$val

  grid::downViewport("layout")
  tmp = as.numeric(grid::grid.locator(unit = "npc"))
  tmp.n <- as.numeric(tmp)*c(max(chromosome_limits_max)*1.33, 1.33*max(plot.data$val))
  diff.a <- (a-tmp.n[1])^2
  diff.b <- (b-tmp.n[2])^2
  grid::upViewport(n = 0)
  paste("Selected SNP is in clump",
        plot.data$clump[which.min(diff.a/max(diff.a) + diff.b/max(diff.b))])
}

#' Identify clump number in \code{\link{selectionResult}} class plot
#'
#' @param x \code{\link{selectionResult}} class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @return No return value, called for side effects
#' @export
identify_clump.selectionResult <- function(x, ...) {
  plot.data <- create_slope_plot_data(x)

  if(length(unique(x$X_info[,3])) == 1){
    chromosome_limits <- aggregate(x$X_info[,4], list(x$X_info[,1]), max)
  } else {
    chromosome_limits <- aggregate(x$X_info[,3], list(x$X_info[,1]), max)
  }

  chromosome_limits_max <- cumsum(chromosome_limits$x)
  chromosome_limits$x <- c(0,head(cumsum(chromosome_limits$x),-1))

  a <- plot.data$snp
  b <- plot.data$val

  grid::downViewport("layout")
  tmp <- as.numeric(grid::grid.locator(unit = "npc"))
  tmp.n <- as.numeric(tmp)*c(max(chromosome_limits_max)*1.33, 1.33*max(plot.data$val))
  diff.a <- (a-tmp.n[1])^2
  diff.b <- (b-tmp.n[2])^2
  grid::upViewport(n = 0)
  paste("Selected SNP is in clump",
        plot.data$clump[which.min((a - tmp.n[1])^2 + (b - tmp.n[2])^2)])
}
