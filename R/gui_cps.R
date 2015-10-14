#' @export
gui_cps <- function() {
  appDir <- system.file("shiny-examples", "myapp", package = "cps")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `cps`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
