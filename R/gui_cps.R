#' GUI for GWAS with SLOPE
#'
#' A graphical user interface for performing Genome-wide
#' Association Study with SLOPE
#'
#' @return null
#' @export
gui_cps <- function() {
  appDir <- system.file("shiny-examples", "myapp", package = "cps")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `cps`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
