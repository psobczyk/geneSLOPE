#' GWAS with clumped SLOPE
#'
#' Files with phenotype and snps are selected interactively.
#'
#' @param phenotypeFile file containing phentotype as its first column
#' @param snpFiles vector of filenames with SNPs e.g. different chromosomes
#' @param pValMax threshold p-value in marginal snp test
#' @param header logical, whether file with phenotype contains header
#' @param sep seperator for file with phenotype
#'
#' @export
#' @details Files with SNPs should be of .raw file format
main <- function(phenotypeFile = NULL, snpFiles = NULL, pValMax = 0.1, header = T, sep = ","){
  if(is.null(phenotypeFile))
    phenotypeFile <- tcltk::tk_choose.files(caption = "Choose file with phenotype")
  phe <- read.table(phenotypeFile, header = header, sep = sep, stringsAsFactors = FALSE)
  y <- phe[,1]
  Filters=matrix(c(".raw file", ".raw"), nrow=1)
  if(is.null(snpFiles))
    snpFiles <- tcltk::tk_choose.files(caption = "Choose SNP files",
                                multi = TRUE, filters = Filters,
                                index = nrow(Filters))
  data_all_files <- NULL
  data_all_files_info <- NULL
  for(file in snpFiles){
    data_single_file <- readPLINK(file)
    # replace missing values with column mean
    data_single_file$snps <- apply(data_single_file$snps, 2, replace_na_with_mean)
    message("Missing values were replaced by column mean")
    #remove all SNPs with no variablity
    nonZeroSd <- apply(data_single_file$snps, 2, sd)!=0
    data_single_file$snps <- data_single_file$snps[,nonZeroSd]
    data_single_file$snpInfo <- data_single_file$snpInfo[nonZeroSd,]
    message(paste(sum(!nonZeroSd), "variables with zero variance were removed"))
    #filter columns with large p-value
    message(paste("Filtering SNPs based on marginal tests.",
                  "Depending on size of data, this may take few minutes"))
    # super optimized p-value computation
    suma = sum((y-mean(y))^2)
    n = length(y) - 2
    pVals <- apply(data_single_file$snps, 2,
                   function(x) pValComp(x,y,n,suma))
    data_single_file$snps <- data_single_file$snps[,pVals<pValMax]
    data_single_file$snpInfo <- data_single_file$snpInfo[pVals<pValMax,]
    message(paste(sum(pVals>=pValMax), "variables with marginal p-value larger than",
                  pValMax, "are excluded from proceeding analysis"))
    data_all_files <- cbind(data_all_files, data_single_file$snps)
    data_all_files_info <- rbind(data_all_files_info, data_single_file$snpInfo)
    message(paste("File:", file, "was sucesfully read"))
  }
  message("All SNP data sucesfully read")

  #clump SNPs to remove highly correlated and to reduce dimenion
  message(paste("Clumping procedure started on", ncol(data_all_files), "snps"))
  clumpedSNPs <- clumpProcedure(y, data_all_files, rho = 0.3)
  message(paste(length(clumpedSNPs$SNPnumber), "clumps extracted"))

  slopeResult <- SLOPE::SLOPE(clumpedSNPs$SNPs, y)
  selectedSNPs <- unlist(clumpedSNPs$SNPnumber)[slopeResult$selected]
  if(is.null(data_all_files_info)){
    result <- colnames(data_all_files)[selectedSNPs]
  } else{
    result <- data_all_files_info[selectedSNPs]
  }

  filename <- paste0("results_clumpedSlope", format(Sys.time(),  "%y_%m_%y_%H_%M_%S"), ".Rdata")
  save(result, file=filename)
  result
}

