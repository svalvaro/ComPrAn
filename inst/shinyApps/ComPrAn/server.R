###################
# server.R
#
###################
server <- function(input, output, session) {

  ########################
  # import tab server side ####
  ########################
  # Import and do mandatory cleaning of the data only after the process button has been pressed:
  peptides_full <- eventReactive(input$processRaw, {
    readFile <- input$inputfile
    if (is.null(readFile)) {
      dataFile <- system.file("extdata", "data.txt", package = "ComPrAn")
      data.table::fread(dataFile)
    } else {
      data.table::fread(readFile$datapath)
    }
  })

  observeEvent(input$inputfile,{
      click("processRaw")
  })
  
  peptides <- reactive({cleanData(peptides_full())})
  
  # Import normalized data only after the process button has been pressed:
  compiledNorm_import <- eventReactive(input$processNorm, {
    readNorm <- input$inputfileNorm
    if (is.null(readNorm)) {
        dataFile <- system.file("extdata", "dataNormProts.txt", package = "ComPrAn")
        protInportForAnalysis(data.table::fread(dataFile))
    } else {
        protInportForAnalysis(data.table::fread(readNorm$datapath))
    }
  })

  vImportMessage <- reactiveValues(data = "No file uploaded and no example file chosed.
                           Please upload a file or click on one of the process buttons in the above tabs.")
  
  observeEvent(input$processRaw, {
      vImportMessage$data <- "Using raw data file. Proceed to part 1."
  })  
  observeEvent(input$processNorm, {
      vImportMessage$data <- "Using normalized vales file. Proceed to part 2."
  })  
  
  output$useCase <- renderText({
      vImportMessage$data
  })
  
  
  # output$useCase <- renderText({
  #   if (input$processRaw != 0) {
  #     "Using raw data file. Proceed to part 1."
  #   } else if (input$processNorm != 0) {
  #     "Using normalized vales file. Proceed to part 2."
  #   } else {
  #     "No file uploaded and no example file chosed. Please upload a file or click on one of the process buttons in the above tabs."
  #   }
  # })

  # output$NormInputTest_0 <- renderText({
  # 
  #   if(input$processNorm != 0) {
  #     # if(input$processRaw != 0) {
  #       # names(normalized_full())
  #     is.object(compiledNorm_import())
  #     # names(peptides())
  #     # print("hello, there.")
  #   } else {
  #     print("no norm file")
  #   }
  # 
  # })


  # # To test a widgets value
  # output$test_widget_value <- renderText({names(peptides())[2]})

  ########################
  # summary tab server side ####
  ########################

  # Plot split between labeled and unlabeled
  output$totalSplit <- renderPlot({

    totalSplit <- data.frame(value = c(nrow(peptides_full()) - nrow(peptides()), nrow(peptides())),
                             label = c("val1", "val2"))
    # totalSplit <- data.frame(value = c(330, 670),
    #                          label = c("val1", "val2"))
    percentUsed <- round(totalSplit$value[2]/sum(totalSplit$value) * 100, 2)

    totalSplit %>%
      ggplot(aes(x = 1, y = value, fill = label)) +
      geom_col(position = "stack") +
      coord_flip(expand = 0) +
      annotate("text", y = totalSplit$value[2], label = paste(percentUsed, "% of total\nused in analysis"), x = 1, hjust = 1.1) +
      scale_x_continuous("") +
      scale_y_continuous("Peptides") +
      scale_fill_manual(values = c("grey90", "skyblue")) +
      theme(axis.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "none",
            rect = element_rect(fill = plotBkg, color = plotBkg))

  })

  # Display number of peptides to the user:
  output$test_input <- renderText({
    if (is.null(input$inputfile)) {
      paste("The example dataset containing", nrow(peptides_full()), " peptides and ", max(peptides()$Fraction) ," fractions will be used. After mandatory cleaning", nrow(peptides()), " (",round(nrow(peptides())/nrow(peptides_full())*100,2) ,"% ) remain. You may now proceed to filtering your data.")
    } else {
      paste0("Your dataset containing ", nrow(peptides_full()), " peptides and ", max(peptides()$Fraction) ," fractions has been imported. After mandatory cleaning", nrow(peptides()), " (",round(nrow(peptides())/nrow(peptides_full())*100,2) ,"% ) remain. You may now proceed to filtering your data.")
    }
  })

  # Plot split between labeled and unlabeled
  output$labUnlabSplit <- renderPlot({
    labUnlabSplit <- data.frame(y = peptides()$isLabel)
    percentLab <- round(table(labUnlabSplit$y)/nrow(labUnlabSplit)*100,2)[2]
    percentUnlab <- round(table(labUnlabSplit$y)/nrow(labUnlabSplit)*100,2)[1]

    labUnlabSplit %>%
      ggplot(aes(x = 1, fill = y)) +
      geom_bar(position = "fill") +
      coord_flip(expand = 0) +
      annotate("text", y = 0, label = paste0(input$labelledName, ", ", percentLab, "%"), x = 1, hjust = -0.1) +
      annotate("text", y = 1, label = paste0(input$unlabelledName, ", ", percentUnlab, "%"), x = 1, hjust = 1.1) +
      scale_x_continuous("") +
      scale_y_continuous("", labels = scales::percent) +
      scale_fill_brewer(palette = "Set1") +
      theme(axis.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = "none",
            rect = element_rect(fill = plotBkg, color = plotBkg))

  })

  ########################
  # filter tab server side ####
  ########################

  output$UI_rank <- renderUI({
      req(peptides())
      sliderInput("rank", label = "Keep peptides ranked below or equal to:",
                  min = 1, max = max(peptides()$Rank), value = 1, step = 1)
  })
  
  # Pre-filtering
  # A in the venn diagrams
  setdiff_LU_unfiltered <- reactive({setdiff(peptides()$`Protein Group Accessions`[peptides()$isLabel],
                                             peptides()$`Protein Group Accessions`[!peptides()$isLabel])
  })

  # B in the venn diagrams
  LU_intersect_unfiltered <- reactive({intersect(peptides()$`Protein Group Accessions`[peptides()$isLabel],
                                                 peptides()$`Protein Group Accessions`[!peptides()$isLabel])
  })

  # C in the venn diagrams
  setdiff_UL_unfiltered <- reactive({setdiff(peptides()$`Protein Group Accessions`[!peptides()$isLabel],
                                             peptides()$`Protein Group Accessions`[peptides()$isLabel])
  })

  # D in the venn diagrams
  LU_union_unfiltered <- reactive({union(peptides()$`Protein Group Accessions`[peptides()$isLabel],
                                         peptides()$`Protein Group Accessions`[!peptides()$isLabel])
  })

  # Display text for venn diagrams, pre:
  output$preFilterVennText <- renderText({

    # Intersect
    B <- length(LU_intersect_unfiltered())

    # Union
    D <- length(LU_union_unfiltered())

    paste("The dataset contains", D, "proteins.", B, "occur in both the", input$labelledName, "and", input$unlabelledName, "samples.")
  })

  # Display venn diagrams, pre:
  output$preFilterVenn <- renderPlot({

    # Labelled but not unlabelled
    A <- length(setdiff_LU_unfiltered())

    # Intersection
    B <- length(LU_intersect_unfiltered())

    # Unlabelled but not labelled
    C <- length(setdiff_UL_unfiltered())

    grid.rect(gp=gpar(fill = plotBkg, col = NA))
    draw.pairwise.venn(A + B,
                       B + C,
                       B,
                       fill = rbCol,
                       lty = 0,
                       category = c(input$labelledName, input$unlabelledName),
                       fontfamily = "sans",
                       cat.fontfamily = "sans",
                       cat.col = rbCol)

  })

  # Filtering...
  peptides_filtered <- eventReactive(input$filter, {

    withProgress(message = 'Filtering and reformatting data...', value = 0, {

      # Step 1:
      n <- 4
      #Sys.sleep(1)
      incProgress(1/n, detail = "Filtering...")
      peptides_filtered <- toFilter(peptides(), rank = as.numeric(input$rank), cl = input$checkGroup)

      # Step 2:
      #Sys.sleep(1)
      incProgress(1/n, detail = "Formating Modifications and labels...")
      peptides_filtered <- splitModLab(peptides_filtered)

      # Step 3:
     # Sys.sleep(1)
      incProgress(1/n, detail = "Simplifying...")
      if(input$simplify) {
        peptides_filtered <- simplifyProteins(peptides_filtered)
      }

      # Step 4:
      #Sys.sleep(1)
      incProgress(1/n, detail = "Ready")
      #Sys.sleep(1)

    })

    return(peptides_filtered)
  })
  
  
  # Post-filtering
  # A in the venn diagrams
  setdiff_LU_filtered <- reactive({setdiff(peptides_filtered()$`Protein Group Accessions`[peptides_filtered()$isLabel],
                                           peptides_filtered()$`Protein Group Accessions`[!peptides_filtered()$isLabel])
  })

  # B in the venn diagrams
  LU_intersect_filtered <- reactive({intersect(peptides_filtered()$`Protein Group Accessions`[peptides_filtered()$isLabel],
                                               peptides_filtered()$`Protein Group Accessions`[!peptides_filtered()$isLabel])
  })

  # C in the venn diagrams
  setdiff_UL_filtered <- reactive({setdiff(peptides_filtered()$`Protein Group Accessions`[!peptides_filtered()$isLabel],
                                           peptides_filtered()$`Protein Group Accessions`[peptides_filtered()$isLabel])
  })

  # D in the venn diagrams
  LU_union_filtered <- reactive({union(peptides_filtered()$`Protein Group Accessions`[peptides_filtered()$isLabel],
                                       peptides_filtered()$`Protein Group Accessions`[!peptides_filtered()$isLabel])})

  # Display text for venn diagrams, post:
  output$postFilterVennText <- renderText({
    B <- length(LU_intersect_filtered())

    D <- length(LU_union_filtered())

    paste("The filtered dataset contains", D, "proteins.", B, "occur in both the", input$labelledName, "and", input$unlabelledName, "samples.")
  })

  # Display venn diagrams, post:
  output$postFilterVenn <- renderPlot({
    # Labelled but not unlabelled
    A <- length(setdiff_LU_filtered())

    # Intersect
    B <- length(LU_intersect_filtered())

    # Unlabelled but not labelled
    C <- length(setdiff_UL_filtered())

    grid.rect(gp=gpar(fill = plotBkg, col = NA))
    draw.pairwise.venn(A + B,
                       B + C,
                       B,
                       fill = rbCol,
                       lty = 0,
                       category = c(input$labelledName, input$unlabelledName),
                       fontfamily = "sans",
                       cat.fontfamily = "sans",
                       cat.col = rbCol)

  })
  
  
  # observeEvent(input$filter,{
  #     output$postFilterVenn = NULL
  #     output$postFilterVennText = NULL
  #     })
  # 
  # observeEvent(input$processRaw,{
  #     output$postFilterVenn = NULL
  #     output$postFilterVennText = NULL
  # })
  
  
  # output$pepsFilteredButton <- renderUI({
  #     if(is.null(peptides_filtered())) {
  #         return(0)
  #     } else {
  #         actionButton("pickPepsNow", "Select peptides")
  #     }
  # })
  observeEvent(peptides_filtered(),{
      output$pepsFilteredButton <- renderUI({
          actionButton("pickPepsNow", "Select peptides")
      })
  })
  
  # observeEvent(input$filter,{
  #     output$pepsFilteredButton <- renderUI({
  #         req(peptides_filtered())
  #         actionButton("pickPepsNow", "Select peptides")
  #     })
  # })
  # 
  
  # Pick representative peptide...
  peptide_index <- eventReactive(input$pickPepsNow, {
      withProgress(message = 'Picking representative peptides...',
                   value = 0, {
                       
                      n <- 3
                       # Sys.sleep(3)
                     incProgress(1/n, detail = "Choosing Representative peptides")
                     peptide_index <- pickPeptide(peptides_filtered())
                    #Sys.sleep(1)
                })
      return(peptide_index)
  })
  
  
  vFilterTabPostVenn <- reactiveValues(data = "")

  observeEvent(input$filter,{
      vFilterTabPostVenn$data <- "show"
  })
  
  observeEvent(input$processRaw,{
      vFilterTabPostVenn$data <- ""
  })


  output$filterTabPostVenn <- renderText({
      vFilterTabPostVenn$data
  })

  outputOptions(output, "filterTabPostVenn", suspendWhenHidden = FALSE)
  
  
  vOpenRepPep <- reactiveValues(data = "")

  observeEvent(input$processNorm,{
      vOpenRepPep$data <- ""
  }
               )
  
  observeEvent(input$processRaw,{
      vOpenRepPep$data <- ""
  }
  )
  
  observeEvent(input$filter,{
      vOpenRepPep$data <- ""
  }
  )

  observeEvent(input$pickPepsNow, {
      req(peptide_index())
      vOpenRepPep$data <- "Representative peptides selected, proceed to next section!"
  })

  output$openRepPep <- renderText({
      vOpenRepPep$data
  })
  
  outputOptions(output, "openRepPep", suspendWhenHidden = FALSE)
  
  
  # output$peptidesSelected <- renderText({
  #     if (!is.null(peptide_index())) {
  #         "Representative peptides selected, proceed to next section!"
  #     }else {
  #         ""
  #     }
  # })


  ########################
  # analysis tab server side ####
  ########################

  proteinLists <- eventReactive(peptide_index(), {
      onlyInOneLabelState_ENV(peptide_index())
  })
  
  # List of filtered values for selecting proteins
  output$dt <- renderUI({
      selectInput("proteinsUnion", label = "Choose protein species",
                  choices = peptides_filtered()$`Protein Group Accessions`,
                  multiple = FALSE)
  })
  
  v <- reactiveValues(data = NULL)

  observeEvent(input$chooseUnlabeled, {
    v$data <- proteinLists()$onlyUnlabelled
  })

  observeEvent(input$chooseBoth, {
    v$data <- proteinLists()$both
  })

  observeEvent(input$chooseLabeled, {
    v$data <- proteinLists()$onlyLabelled
  })

  output$trace_table <- renderDataTable({
    if (is.null(iris)) return()
    DT::datatable(iris, options = list(paging = FALSE))
  })

  # # Display table for selections
  # output$x1 = DT::renderDataTable(cars, server = FALSE)
  # 
  # # highlight selected rows in the scatterplot
  # output$x2 = renderPlot({
  #   s = input$x1_rows_selected
  #   par(mar = c(4, 4, 1, .1))
  #   plot(cars)
  #   if (length(s)) points(cars[s, , drop = FALSE], pch = 19, cex = 2)
  # })

  # print info text:
  # To test an object's value
  # output$test_object_value <- renderText({length(peptides_plot())})
  output$test_object_value <- renderText({length(v$data)})
  output$test_object_value_2 <- renderText({length(peptides_full_parameters())})

  # Display table for selections
  # output$allPeptides_choose = DT::renderDataTable(data.frame(Protein = peptides_plot()), server = FALSE)
  # output$allPeptides_choose = DT::renderDataTable(data.frame(Proteins = peptides_plot()),
  #                                                 options = list(paging = FALSE),
  #                                                 server = FALSE,
  #                                                 selection = "single")

  output$allPeptides_choose = DT::renderDataTable({
    if (is.null(v$data)) {
      return()
    } else {
      DT::datatable(data.frame(Proteins = v$data),
                    options = list(paging = FALSE),
                    selection = "single")
    }
  }, server = FALSE
  )

  # allPeptidesPlot
  output$allPeptidesPlot = renderPlot({

    req(peptide_index())

    if (input$tabset1 == "All Proteins") {
      if (is.null(input$proteinsUnion)) {
        return()
      } else {
        protein <- input$proteinsUnion
      }
    } else if (is.null(input$allPeptides_choose_rows_selected)) {
      return()
    } else {
      protein <- v$data[input$allPeptides_choose_rows_selected]
    }
    req(protein)
    allPeptidesPlot(peptide_index(), protein, max(peptides()$Fraction),
                    meanLine = input$allPeptidesPlot_mean,
                    repPepLine = input$allPeptidesPlot_reppep,
                    separateLabStates = input$allPeptidesPlot_sepstates,
                    grid = input$allPeptidesPlot_removegrid,
                    titleLabel = input$allPeptidesPlot_title,
                    titleAlign = input$allPeptidesPlot_align,
                    xlabel = input$allPeptidesPlot_xaxis,
                    ylabel = input$allPeptidesPlot_yaxis,
                    labelled = input$labelledName,
                    unlabelled = input$unlabelledName)

  })

  ## Allow user to downloadtables of proteins in both/only one data set
  output$downloadDataBoth <- downloadHandler(
      filename = function() {
          paste0("proteinsIn_Both_samples.txt")
      },
      content = function(file) {
          # Write to a file specified by the 'file' argument
          write(proteinLists()$both, file, sep = "\n")
      }
  )
  
  output$downloadDataLabeled <- downloadHandler(
      filename = function() {
          paste0("proteinsOnlyIn_",input$labelledName, "_sample.txt")
      },
      content = function(file) {
          # Write to a file specified by the 'file' argument
          write(proteinLists()$onlyLabelled, file, sep = "\n")
      }
  )
  
  output$downloadDataUnlabeled <- downloadHandler(
      filename = function() {
          paste0("proteinsOnlyIn_",input$unlabelledName, "_sample.txt")
      },
      content = function(file) {
          # Write to a file specified by the 'file' argument
          write(proteinLists()$onlyUnlabelled, file, sep = "\n")
      }
  )
  
  # Display download link when data are available
  
  output$dl_both <- renderUI({
      req(proteinLists())
      downloadLink("downloadDataBoth", "List of proteins present in both samples")
  })
  output$dl_onlyLab <- renderUI({
      req(proteinLists())
      downloadLink("downloadDataLabeled", "List of proteins present only in labelled sample")
  })
  output$dl_onlyUnlab <- renderUI({
      req(proteinLists())
      downloadLink("downloadDataUnlabeled", "List of proteins present only in unlabelled sample")
  })
  # 
  ###############################
  # Normalization tab server side ####
  ###############################

  compiledExport <- eventReactive(input$normData,{
      getNormTable(peptide_index(),purpose = "export")
    })

  compiledNorm <- eventReactive(input$normData,{
      getNormTable(peptide_index(),purpose = "analysis")
    })

  
  observeEvent(compiledExport(),{
      output$dl_Norm <- renderUI({
              downloadLink("downloadNormData", "Download Normalized Values")
      })
  })
  # output$dl_Norm <- renderUI({
  #   if(is.null(compiledExport())) {
  #     return(0)
  #   } else {
  #     downloadLink("downloadNormData", "Download Normalized Values")
  #   }
  # })

  output$downloadNormData <- downloadHandler(
    filename = function() {
      "NormalizedProteinData.txt"
    },
    content = function(file) {
      export(compiledExport(), file)
    }
  )
  
  vNormProtDownload <- reactiveValues(data = "")
  
  observeEvent(input$normData,{
      vNormProtDownload$data <- "show"
  })
  
  observeEvent(input$processRaw,{
      vNormProtDownload$data <- ""
  })
  
  observeEvent(input$pickPepsNow,{
      vNormProtDownload$data <- ""
  })
  

  output$normProtDownload <- renderText({
      vNormProtDownload$data
  })
  
  outputOptions(output, "normProtDownload", suspendWhenHidden = FALSE)

  # output$NormTest <- renderText({
  #   if(is.null(input$normData)) {
  #     return()
  #   } else {
  #     names(compiledNorm())
  #   }
  # 
  # })

  vOpenPart2 <- reactiveValues(data = "")
  
  observeEvent(input$processNorm,{
      vOpenPart2$data <- "Part2 ready!"
  }
  )
  
  observeEvent(input$normData, {
      vOpenPart2$data <- "Part2 ready!"
  })
  
  observeEvent(input$processRaw, {
      vOpenPart2$data <- ""
  })
  
  observeEvent(input$filter,{
      vOpenPart2$data <- ""
  }
  )
  
  observeEvent(input$pickPepsNow,{
      vOpenPart2$data <- ""
  }
  )
  
  output$openPart2 <- renderText({
      vOpenPart2$data
  })
  
  outputOptions(output, "openPart2", suspendWhenHidden = FALSE)
  
  ###############################
  # proteinNormViz tab server side ####
  ###############################

  # Check for a dataset. if one has been imported, or if the example file is being used
  # then take it. If not, use the result from part 1
  compiledNorm_plot <- reactive({
    if(input$processNorm != 0) {
      compiledNorm_import()
    } else {
      compiledNorm()
    }

  })
  
  # vNormProts <- reactiveValues(data = NULL)
  # 
  # observeEvent(compiledNorm_plot(), {
  #     vNormProts$data <- TRUE
  # })
  # observeEvent(input$processNorm, {
  #     vNormProts$data <- TRUE
  # }) 
  # 
  # output$normDataPresent <- reactive({
  #     vNormProts$data
  # })
  
  # 
  # observeEvent(compiledNorm_plot(), {
  #     vNormProts$data <- "dataIn"
  # })
  # observeEvent(input$processNorm, {
  #     vNormProts$data <- "dataIn"
  # }) 
  # 
  # output$normDataPresent <- renderText({
  #     vNormProts$data
  # })
  # outputOptions(output, "normDataPresent", suspendWhenHidden = FALSE)
  # Just some output testing here
  output$NormInputTest <- renderText({
    length(unique(compiledNorm_plot()[compiledNorm_plot()$scenario == "B",]$`Protein Group Accessions`))
    })

  # List of filtered values for selecting proteins
  output$dt_2 <- renderUI({
    selectInput("proteinNormChoose", label = "Choose protein species",
                choices = unique(compiledNorm_plot()[compiledNorm_plot()$scenario == "B",]$`Protein Group Accessions`),
                multiple = FALSE)
  })

  # Produce the plot as a reactive object

  normPlotSingle <- reactive({
    req(input$proteinNormChoose)

    if (is.null(input$proteinNormChoose)) {
      return()
    } else {
      protein <- input$proteinNormChoose
    }

    # proteinPlot(compiledNorm_plot()[compiledNorm_plot()$scenario == "B",], protein, max(peptides()$Fraction),
       proteinPlot(compiledNorm_plot(), protein, as.numeric(max(compiledNorm_plot()$Fraction)),
 
                grid = input$allProteinPlot_removegrid,
                titleLabel = input$allProteinPlot_title,
                titleAlign = input$allProteinPlot_align,
                xlabel = input$allProteinPlot_xaxis,
                ylabel = input$allProteinPlot_yaxis,
                legendLabel = input$allProteinPlot_legend,
                labelled = input$labelledName,
                unlabelled = input$unlabelledName)
  })

  # Send the reactive plot to the UI
  output$proteinPlot = renderPlot({
      normPlotSingle()
  })

  # Display download link only if a protein is selected
  output$dl_Norm_Plot <- renderUI({
    if(is.null(input$proteinNormChoose)) {
      return(0)
    } else {
      downloadLink("downloadNormPlotSingle", "Download plot")
    }
  })

  # Download handler for the reactive normalized single protein plot
  output$downloadNormPlotSingle <- downloadHandler(
    filename = function() {
      paste0(input$proteinNormChoose_normValues, ".pdf")
    },
    content = function(file) {
      ggsave(file, normPlotSingle())
    }
  )


  # Input for multiple normalized plots
  normPlotsMultipleInput <- reactive({
    strsplit(x = input$normPlotsMultipleInput, split = "\\s")[[1]]
  })

  # Save multiplots to a reactive object (list)
  g_norm <- reactive({
    normPlotsMultipleInput() %>%
    # c("Q16540", "P52815", "P09001", "Q13405", "Q9H2W6", "Q9NYK5", "Q96DV4") %>%
      map(~ proteinPlot(compiledNorm_plot(), ., as.numeric(max(compiledNorm_plot()$Fraction)),

                        grid = input$allProteinPlot_removegrid,
                        titleLabel = input$allProteinPlot_title,
                        titleAlign = input$allProteinPlot_align,
                        xlabel = input$allProteinPlot_xaxis,
                        ylabel = input$allProteinPlot_yaxis,
                        legendLabel = input$allProteinPlot_legend,
                        labelled = input$labelledName,
                        unlabelled = input$unlabelledName))
  })

      output$allgraphsNorm = downloadHandler(
        filename = 'multiPlotsNorm.pdf',
        content = function(file) {
          ggsave(file, marrangeGrob(grobs = g_norm(), nrow=1, ncol=1))
          })

  ###############################
  # heatMaps tab server side ####
  ###############################
  # Just some test output for me
  output$HeatTest <- renderText({
    input$heatMapFile$datapath })

  # df <- eventReactive(input$heatMapFile,
  #       {
  #           tryCatch(
  #               {
  #                   df <- read_tsv(input$heatMapFile$datapath)
  #               },
  #               error = function(e) {
  #                   # return a safeError if a parsing error occurs
  #                   stop(safeError(e))
  #               }
  #           )
  #       })
  # 
  # df <- eventReactive(input$exampleGroup,
  #                     {
  #                         df <- read_tsv(system.file("extdata", "exampleGroup.txt", package = "ComPrAn"))
  #                     })
  # 
  # # Make reactive plot object
  # heatMapPlotObject <- reactive({
  #     req(df())
  #     
  #     if (is.null(input$heatMapGroupNameColumn)) {
  #         groupHeatMap(compiledNorm_plot()[compiledNorm_plot()$scenario == "B",], 
  #                      df(), input$heatMapGroupName,
  #                      titleAlign = "center",
  #                      grid = F, colNumber = input$showSamplesHeatMap,
  #                      labelled = input$labelledName,
  #                      unlabelled = input$unlabelledName)
  #     } else {
  #         groupHeatMap(compiledNorm_plot()[compiledNorm_plot()$scenario == "B",], 
  #                      df(), input$heatMapGroupName,
  #                      titleAlign = "center",
  #                      newNamesCol = input$heatMapGroupNameColumn,
  #                      grid = F, colNumber = input$showSamplesHeatMap,
  #                      labelled = input$labelledName,
  #                      unlabelled = input$unlabelledName)
  #     }
  #     
  # })
  
  ### reactive value to get a list of all column names in group data frame
  vPresentColumns <- reactiveValues(data = NULL)
  ## read in group data frame
  vGroupDF <- reactiveValues(data = NULL)

  observeEvent(input$exampleGroup,{
      vGroupDF$data <- read_tsv(system.file("extdata", "exampleGroup.txt", package = "ComPrAn"))
      vPresentColumns$data <- names(vGroupDF$data)
  })
  
  observeEvent(input$heatMapFile,{
      tryCatch(
          {
              vGroupDF$data <- read_tsv(input$heatMapFile$datapath)
          },
          error = function(e) {
              # return a safeError if a parsing error occurs
              stop(safeError(e))
          }
      )
      vPresentColumns$data <- names(vGroupDF$data)
  })
  

  
  output$testCols <- renderText(
      vPresentColumns$data
  )
      

  #     
  # observeEvent(input$exampleGroup,{
  #     df <- read_tsv(system.file("extdata", "exampleGroup.txt", package = "ComPrAn"))
  #     vPresentColumns$data <- names(df)
  # })
  # 
  # observeEvent(input$heatMapFile,{
  #     tryCatch(
  #         {
  #             df <- read_tsv(input$heatMapFile$datapath)
  #         },
  #         error = function(e) {
  #             # return a safeError if a parsing error occurs
  #             stop(safeError(e))
  #         }
  #     )
  #     vPresentColumns$data <- names(df)
  # })
  
  # output$testColNames <- renderText({
  #     vPresentColumns$data
  #   })
  ### reactive value to return as a new column name
  
  vHeatmapGroupColumn <- reactiveValues(data = NULL)

  # observeEvent(input$heatMapGroupNameColumn,{
  #     if(input$heatMapGroupNameColumn %in%  vPresentColumns$data){
  #         vHeatmapGroupColumn$data <- input$heatMapGroupNameColumn
  #         output$colNotFound <- renderText({
  #             NULL})
  #         
  #     } 
  #     #if(!input$heatMapGroupNameColumn %in%  vPresentColumns$data)
  #      else   {
  #         output$colNotFound <- renderText({
  #             paste("Note: Column", input$heatMapGroupNameColumn,
  #                 "is not present in group data!",sep = " ")
  #     })}
  # 
  # })
  # # observeEvent(input$renameProteinsHeatMap ,{
  #     output$colNotFound <- renderText({
  #         NULL
  #     })
  # })
  
  observeEvent(input$renameProteinsHeatMap ,{
      vHeatmapGroupColumn$data <- NULL
  })
  
  # Text input with new column name
  # output$dt_3 <- renderUI({
  #     if (input$renameProteinsHeatMap) {
  #         textInput("heatMapGroupNameColumn", label = "Column with new names", value = vHeatmapGroupColumn$data)
  #     }
  # })
  
  #List of columns to select a column to rename the data with
  output$HeatmapGroupColList <- renderUI({
      if (input$renameProteinsHeatMap) {
      selectInput("heatGroupNameCol", label = "Choose column with new names",
                  choices = vPresentColumns$data,
                  multiple = FALSE)
      }
  })
  
  observe({
      if (input$renameProteinsHeatMap) {
          vHeatmapGroupColumn$data <- input$heatGroupNameCol
      }
  })
  
  
  
  vHeatmapOrderColumn <- reactiveValues(data = NULL)
  
  observeEvent(input$reorderProteinsHeatMap ,{
      vHeatmapOrderColumn$data <- NULL
  })
  

  #List of columns to select a column to reorder the data with
  output$HeatmapGroupColList_2 <- renderUI({
      req(vPresentColumns$data)
      if (input$reorderProteinsHeatMap) {
          selectInput("heatGroupOrderCol", label = "Choose column with order",
                      choices = vPresentColumns$data,
                      multiple = FALSE)
      }
  })
  
  observeEvent(input$heatGroupOrderCol,{
          req(vGroupDF$data)
         if(typeof(vGroupDF$data[[input$heatGroupOrderCol]]) == "double"){
          output$testColType <- renderText({typeof(vGroupDF$data[[input$heatGroupOrderCol]])})
          vHeatmapOrderColumn$data <- input$heatGroupOrderCol
         }
          

  })
  
  
  output$testDF <- renderText({
      #req(input$reorderProteinsHeatMap)
      vHeatmapOrderColumn$data
  })
  
  # Make reactive plot object from example group
  heatMapPlotObject_example <- reactive({
      # req(input$exampleGroup)
      # df <- read_tsv(system.file("extdata", "exampleGroup.txt", package = "ComPrAn"))
      req(vGroupDF$data)
      groupHeatMap(compiledNorm_plot()[compiledNorm_plot()$scenario == "B",],
                   vGroupDF$data, input$heatMapGroupName,
                   titleAlign = "center",
                   newNamesCol = vHeatmapGroupColumn$data,
                   grid = F, colNumber = input$showSamplesHeatMap,
                   labelled = input$labelledName,
                   unlabelled = input$unlabelledName)
      
  })
  
  # Make reactive plot object
  heatMapPlotObject <- reactive({
      # req(input$heatMapFile)
      # tryCatch(
      #     {
      #         df <- read_tsv(input$heatMapFile$datapath)
      #     },
      #     error = function(e) {
      #         # return a safeError if a parsing error occurs
      #         stop(safeError(e))
      #     }
      # )
      # 
      
      req(vGroupDF$data)
      groupHeatMap(compiledNorm_plot()[compiledNorm_plot()$scenario == "B",],
                   vGroupDF$data, input$heatMapGroupName,
                   titleAlign = "center",
                   newNamesCol = vHeatmapGroupColumn$data,
                   grid = F, colNumber = input$showSamplesHeatMap,
                   labelled = input$labelledName,
                   unlabelled = input$unlabelledName,
                   orderColumn = vHeatmapOrderColumn$data)
      
      
  })

  # Make reactive plot object from example group - currently in use
  # heatMapPlotObject_example <- reactive({
  #     req(input$exampleGroup)
  #     df <- read_tsv(system.file("extdata", "exampleGroup.txt", package = "ComPrAn"))
  #     
  #     if (is.null(input$heatMapGroupNameColumn)) {
  #         groupHeatMap(compiledNorm_plot()[compiledNorm_plot()$scenario == "B",],
  #                      df, input$heatMapGroupName,
  #                      titleAlign = "center",
  #                      grid = F, colNumber = input$showSamplesHeatMap,
  #                      labelled = input$labelledName,
  #                      unlabelled = input$unlabelledName)
  #     } else {
  #         groupHeatMap(compiledNorm_plot()[compiledNorm_plot()$scenario == "B",],
  #                      df, input$heatMapGroupName,
  #                      titleAlign = "center",
  #                      newNamesCol = input$heatMapGroupNameColumn,
  #                      grid = F, colNumber = input$showSamplesHeatMap,
  #                      labelled = input$labelledName,
  #                      unlabelled = input$unlabelledName)
  #     }
  #     
  # })

  
  observeEvent(input$exampleGroup,{
               output$heatMapPlot = renderPlot({
                   heatMapPlotObject_example()
               })
  })
  
  observeEvent(input$heatMapFile,{
               output$heatMapPlot = renderPlot({
                   heatMapPlotObject()
               })
  })
  
  observeEvent(input$processNorm,{
      output$heatMapPlot = NULL
  })
  # # Send heatmap to GUI
  # output$heatMapPlot = renderPlot({
  #     heatMapPlotObject()
  # })
  # 
  # # Send heatmap to GUI
  # output$heatMapPlot_example = renderPlot({
  #     heatMapPlotObject_example()
  # })
  
  # Save heatmap
  # Display download link only if a protein is selected
  output$dl_Heat_Plot <- renderUI({
    if(is.null(heatMapPlotObject())) {
      return(0)
    } else {
      downloadLink("downloadHeatPlot", "Download Heatmap")
    }
  })

  # Download handler for the reactive normalized single protein plot
  output$downloadHeatPlot <- downloadHandler(
    filename = function() {
      "heatmap.pdf"
    },
    content = function(file) {
      ggsave(file, heatMapPlotObject())
    }
  )

  ###############################
  # coMigration tab server side ####
  ###############################

  # Comigration plot 1:
  # Comigration plot 2, reactive:
  coMig_2_plot <- reactive({
    req(input$groupData_coMig2_g1)
    req(input$groupData_coMig2_g2)

    twoGroupsWithinLabelCoMigration(compiledNorm_plot(), as.numeric(max(compiledNorm_plot()$Fraction)),
                                    # group1Data = group1DataVector,
                                    # group1Name = group1Name,
                                    # group2Data = group2DataVector,
                                    # group2Name = group2Name,

                                    group1Data = strsplit(x = input$groupData_coMig2_g1, split = "\\s")[[1]],
                                    group1Name = input$groupName_coMig2_g1,

                                    group2Data = strsplit(x = input$groupData_coMig2_g2, split = "\\s")[[1]],
                                    group2Name = input$groupName_coMig2_g2,

                                    grid = input$grid_coMig2,
                                    meanLine = input$meanLine_coMig2,
                                    medianLine = input$medianLine_coMig2,

                                    jitterPoints = input$jitterPoints_coMig2,
                                    pointSize = input$pointSize_coMig2,
                                    alphaValue = input$alphaValue_coMig2,
                                    titleAlign = input$titleAlign_coMig2,

                                    ylabel = input$ylabel_coMig2,
                                    xlabel = input$xlabel_coMig2,
                                    legendLabel = input$legendLabel_coMig2,

                                    labelled = input$labelledName,
                                    unlabelled = input$unlabelledName
    )

  })

  # Send plot to GUI
  
  vComig <- reactiveValues(data = NULL)
  
  observeEvent(input$groupData_coMig1, {
      vComig$data <- strsplit(x = input$groupData_coMig1, split = "\\s")[[1]]
  })
  
  observeEvent(input$mtLSU, {
      vComig$data <- mtLSUProts
  })  
  observeEvent(input$mtSSU, {
      vComig$data <- mtSSUProts
  })  
  
  coMig1Plot_proteins <- reactive({
      if(is.null(vComig$data)){
          return(NULL)
      } else {
          return(vComig$data)
      } 
  })

  coMig_1_plot <- reactive({
     
    req(coMig1Plot_proteins())

    oneGroupTwoLabelsCoMigration(compiledNorm_plot(), as.numeric(max(compiledNorm_plot()$Fraction)),
                                                         # groupData = groupDataVector,
                                                         # groupName = groupName,

                                                         groupData = coMig1Plot_proteins(),
                                                         groupName = input$groupName_coMig1,

                                                         grid = input$grid_coMig1,
                                                         meanLine = input$meanLine_coMig1,
                                                         medianLine = input$medianLine_coMig1,

                                                         jitterPoints = input$jitterPoints_coMig1,
                                                         pointSize = input$pointSize_coMig1,
                                                         alphaValue = input$alphaValue_coMig1,
                                                         titleAlign = input$titleAlign_coMig1,

                                                         ylabel = input$ylabel_coMig1,
                                                         xlabel = input$xlabel_coMig1,
                                                         legendLabel = input$legendLabel_coMig1,

                                                         labelled = input$labelledName,
                                                         unlabelled = input$unlabelledName)
    })

  output$coMig_1 <- renderPlot({
    coMig_1_plot()
  })

  # Save comig2
  # Display download link only if a protein is selected
  output$dl_Comig1_Plot <- renderUI({
    if(is.null(coMig_1_plot())) {
      return(0)
    } else {
      downloadLink("downloadComig1Plot", "Download comigration (type 1) plot")
    }
  })

  # Download handler for the reactive normalized single protein plot
  output$downloadComig1Plot <- downloadHandler(
    filename = function() {
      "comigration1.pdf"
    },
    content = function(file) {
      ggsave(file, coMig_1_plot())
    }
  )

  # Comigration plot 2, reactive:
  coMig_2_plot <- reactive({
    req(input$groupData_coMig2_g1)
    req(input$groupData_coMig2_g2)

    twoGroupsWithinLabelCoMigration(compiledNorm_plot(), as.numeric(max(compiledNorm_plot()$Fraction)),
                                    # group1Data = group1DataVector,
                                    # group1Name = group1Name,
                                    # group2Data = group2DataVector,
                                    # group2Name = group2Name,

                                    group1Data = strsplit(x = input$groupData_coMig2_g1, split = "\\s")[[1]],
                                    group1Name = input$groupName_coMig2_g1,

                                    group2Data = strsplit(x = input$groupData_coMig2_g2, split = "\\s")[[1]],
                                    group2Name = input$groupName_coMig2_g2,

                                    grid = input$grid_coMig2,
                                    meanLine = input$meanLine_coMig2,
                                    medianLine = input$medianLine_coMig2,

                                    jitterPoints = input$jitterPoints_coMig2,
                                    pointSize = input$pointSize_coMig2,
                                    alphaValue = input$alphaValue_coMig2,
                                    titleAlign = input$titleAlign_coMig2,

                                    ylabel = input$ylabel_coMig2,
                                    xlabel = input$xlabel_coMig2,
                                    legendLabel = input$legendLabel_coMig2,

                                    labelled = input$labelledName,
                                    unlabelled = input$unlabelledName
    )

  })

  # Send plot to GUI
  output$coMig_2 <- renderPlot({
    coMig_2_plot()
  })

  # Save comig2
  # Display download link only if a protein is selected
  output$dl_Comig2_Plot <- renderUI({
    if(is.null(coMig_2_plot())) {
      return(0)
    } else {
      downloadLink("downloadComig2Plot", "Download comigration (type 2) plot")
    }
  })

  # Download handler for the reactive normalized single protein plot
  output$downloadComig2Plot <- downloadHandler(
    filename = function() {
      "comigration2.pdf"
    },
    content = function(file) {
      ggsave(file, coMig_2_plot())
    }
  )

  ###############################
  # cluster tab server side ####
  ###############################
  
  ## Create components necessary for clustering
  clusteringDF <- eventReactive(input$distCentered,{
      clusterComp(compiledNorm_plot(),
                              scenar = "A", 
                              PearsCor = input$distCentered)})
  
  
  # Generate slider for cutoff
  output$UI_distCutoff <- renderUI({
      req(clusteringDF())
      
      sliderInput("distCutoff", "Cutoff for distance matrix",
                  min = round(min(clusteringDF()$labDistM,
                                  clusteringDF()$unlabDistM),2),
                  max = round(max(clusteringDF()$labDistM,
                                  clusteringDF()$unlabDistM),2),
                  value = round(max(clusteringDF()$labDistM,
                                    clusteringDF()$unlabDistM)/2,2),
                  step = 0.05)
  })
  
  labelledTable_clust <- reactive({
      req(input$distCutoff)
      assignClusters(.listDf = clusteringDF(),
                                                  sample = "labeled", 
                                        method = input$distMethod, 
                                        cutoff = input$distCutoff)})

  unlabelledTable_clust <- reactive({
      req(input$distCutoff)
      assignClusters(.listDf = clusteringDF(),
                                                    sample = "unlabeled", 
                                                    method = input$distMethod, 
                                                    cutoff = input$distCutoff)})

  #make bar plots summarizing numbers of proteins per cluster
  labeledBar <- reactive({
    req(labelledTable_clust())

    makeBarPlotClusterSummary(labelledTable_clust(), name = input$labelledName)
    })

  unlabeledBar <- reactive({
    req(unlabelledTable_clust())

    makeBarPlotClusterSummary(unlabelledTable_clust(), name = input$unlabelledName)
    })

  # Send plots to GUI
  output$labeledBar_plot <- renderPlot({
    labeledBar()
  })
  output$unlabeledBar_plot <- renderPlot({
    unlabeledBar()
  })

  # GUI for labeledBar_plot
  # Display download link only if a protein is selected
  output$dl_labeledBar_Plot <- renderUI({
    if(is.null(labeledBar())) {
      return(0)
    } else {
      downloadLink("downloadlabeledBar_Plot", "Download plot")
    }
  })

  # GUI for unlabeledBar_plot
  # Display download link only if a protein is selected
  output$dl_unlabeledBar_Plot <- renderUI({
    if(is.null(unlabeledBar())) {
      return(0)
    } else {
      downloadLink("downloadunlabeledBar_Plot", "Download plot")
    }
  })

  # Save labeledBar_plot
  # Download handler for the reactive normalized single protein plot
  output$downloadlabeledBar_Plot <- downloadHandler(
    filename = function() {
      "labeledClusters.pdf"
    },
    content = function(file) {
      ggsave(file, labeledBar())
    }
  )

  # Save unlabeledBar_plot
  # Download handler for the reactive normalized single protein plot
  output$downloadunlabeledBar_Plot <- downloadHandler(
    filename = function() {
      "labeledClusters.pdf"
    },
    content = function(file) {
      ggsave(file, unlabeledBar())
    }
  )



  # Datatable generation and download
  #create table for export
  tableForClusterExport <- reactive({exportClusterAssignments(labelledTable_clust(),unlabelledTable_clust())})

  # GUI for cluster table
  # Display download link only if a protein is selected
  output$dl_clustertable <- renderUI({
    req(tableForClusterExport())
    if(is.null(tableForClusterExport())) {
      return(0)
    } else {
      downloadLink("downloadClusters", "Download table of cluster IDs")
    }
  })

  # Save unlabeledBar_plot
  # Download handler for the reactive normalized single protein plot
  output$downloadClusters <- downloadHandler(
    filename = function() {
      "Clusters.txt"
    },
    content = function(file) {
      export(tableForClusterExport(), file)
    }
  )


}
