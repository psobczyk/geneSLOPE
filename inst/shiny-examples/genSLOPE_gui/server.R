
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(cps)

options(shiny.maxRequestSize=300*1024^2)

shinyServer(function(input, output) {

  phenotype <- reactive({
    name <- input$fileY
    validate(
      need(name != "", "Please upload phenotype data")
    )
    readPhenotype(name$datapath, sep=input$sep, header = input$header)
  })

  output$phenotypeOk <- reactive({
    length(phenotype()$y)
  })

  screening <- eventReactive(input$go, {
    DataSet <- input$file
    Map <- input$map.file
    validate(
      need(DataSet != "", "Please upload snp data")
    )
    readSNPs(DataSet$datapath, Map$datapath, y = phenotype()$y,
                pValMax = input$pValCutoff, chunk_size = 1e2, verbose = FALSE)
  })

  output$summary <- renderText({
    paste("Phenotype data loaded.", length(phenotype()$y), "observations.")
  })

  clumping <-  eventReactive(input$go, {
    clumpProcedure(screening(), input$rho, verbose = FALSE)
  })

  slopeResult <- eventReactive(input$go, {
    genSLOPE(clumping(), input$fdr, verbose = FALSE)
  })

  output$clumpSummary <- renderPrint({
    summary(clumping())
  })

  output$slopePlot <- renderPlot({
    plot(slopeResult())
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("results_genSlope", format(Sys.time(),  "%y_%m_%y_%H_%M_%S"), ".csv")
    },
    content = function(file) {
      write.csv(slopeResult()$X, file, row.names = FALSE)
    }
  )
})
