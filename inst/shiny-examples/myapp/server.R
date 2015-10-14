
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(cps)

options(shiny.maxRequestSize=50*1024^2)

shinyServer(function(input, output) {

  phenotype <- reactive({
    name <- input$fileY
    validate(
      need(name != "", "Please upload phenotype data")
    )
    phenotype <- read.table(name$datapath, header = input$header,
                                     sep = ",", stringsAsFactors = FALSE)
    phenotype
  })

  output$phenotypeOk <- reactive({
    ncol(phenotype())
  })

  data <- eventReactive(input$go, {
    DataSets <- input$file
    validate(
      need(DataSets != "", "Please upload snp data")
    )
    data_all_files <- NULL
    for(p in 1:nrow(input$file)){
      data_single_file <- cps:::readPLINK(DataSets[[p, 'datapath']])
      data_single_file$snps <- apply(data_single_file$snps, 2, cps:::replace_na_with_mean)
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
                     function(x) cps:::pValComp(x,y,n,suma))
      data_single_file$snps <- data_single_file$snps[,pVals<input$pValCutoff]
      data_all_files <- cbind(data_all_files, data_single_file$snps)
    }
    data_all_files
  })

  output$summary <- renderText({
    paste("Phenotype data loaded.", nrow(phenotype()), "observations.")
  })

  slopeResult <- reactive({
    y <- phenotype()[,1]
    X <- data()
    validate(
      need(dim(X) == length(y), "Please upload snp data")
    )
    clumpedSLOPE(y, X)
  })

  output$slopeTable <- renderDataTable({
    slopeResult()
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("results_clumpedSlope", format(Sys.time(),  "%y_%m_%y_%H_%M_%S"), ".csv")
    },
    content = function(file) {
      write.csv(slopeResult(), file, row.names = FALSE)
    }
  )
})
