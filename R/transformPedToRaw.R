#' Transformin .ped to .raw using plink command line
#'
#' @param file, character, name of .raw file
transformPedToRaw <- function(file){
  if(!file.exists(file))
    error(paste0("Error. ", file, " - no such file."))

  splittedName <- strsplit(file, "\\.")[[1]]
  if(tail(splittedName,1)!="ped")
    error(paste0("Error. ", file, " - .ped file is required."))

  file = splittedName[1]
  command <- paste("./plink --file", file, "--out", file,
                   "--recodeA --noweb")
  message("Transforming .ped file into R compatible .raw")
  system(command)
  message("File transformed")
}
