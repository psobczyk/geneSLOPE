#' Transforming .ped to .raw using PLINK executable
#'
#' Performs PLINK in command line to transform .ped
#' into R compatible .raw file.
#'
#' @export
#' @param plinkExec path to PLINK executable
#' @param file name of .ped file
transformPedToRaw <- function(plinkExec = NULL, file){
  if(is.null(plinkExec))
    stop("Path to PLINK needs to be specified")
  if(!file.exists(file))
    stop(paste0("Error. ", file, " - no such file."))

  splittedName <- strsplit(file, "\\.")[[1]]
  if(tail(splittedName,1)!="ped")
    stop(paste0("Error. ", file, " - .ped file is required."))

  file = paste(head(splittedName, -1), collapse = ".")
  command <- paste0("./", plinkExec, " --file ", file,
                    " --out ", file, " --recodeA --noweb")
  print(command)
  message("Transforming .ped file into R compatible .raw")
  system(command)
  message("File transformed")
}
