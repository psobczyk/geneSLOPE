#' phenotypeData class
#'
#' Phenotype data
#'
#' @details Always a named list of ten elements
#' \enumerate{
#' \item \code{y} numeric vector, phenotype
#' \item \code{yInfo} matrix with additional information about observations
#' provied in .fam file
#' }
#'
#' @seealso \code{\link{readPhenotype}}
#' @name phenotypeData
NULL

#' Print phenotypeData class object
#'
#' @param x phenotypeData class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method print phenotypeData
print.phenotypeData <- function(x, ...){
  cat("Object of class phenotypeData\n")
  cat("$y: vector of size", length(x$y), "\n")
  cat("$yInfo: matrix with information\n")
  cat("\t", nrow(x$yInfo), " observations\n")
  cat("\t", ncol(x$yInfo), " variables\n")
}

#' Summary phenotypeData class object
#'
#' @param object phenotypeData class object
#' @param ... Further arguments to be passed to or from other methods. They are ignored in this function.
#' @export
#'
#' @method summary phenotypeData
summary.phenotypeData <- function(object, ...){
  cat("Object of class phenotypeData\n")
  cat("$y: vector of size", length(x$y), "\n")
  cat("$yInfo: matrix with information\n")
  cat("\t", nrow(x$yInfo), " observations\n")
  cat("\t", ncol(x$yInfo), " variables\n")
}
