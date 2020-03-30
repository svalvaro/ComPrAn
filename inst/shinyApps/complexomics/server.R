###################
# server.R
#
###################
server <- function(input, output, session) {

  ########################
  # import tab server side
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

  peptides <- reactive({cleanData(peptides_full())})

  # Import normalized data only after the process button has been pressed:
  normalized_full <- eventReactive(input$processNorm, {
    readNorm <- input$inputfileNorm
    if (is.null(readNorm)) {
      readr::read_tsv("./data/NormValues.txt")
    } else {
      readr::read_tsv(readNorm$datapath)

    }
  })

  compiledNorm_import <- reactive({normalized_full() %>%
      gather(Fraction, `Precursor Area`, -c(`Protein Group Accessions`, `Protein Descriptions`, scenario, label)) %>%
      rename("isLabel" = "label") %>%
      select(`Protein Group Accessions`, `Protein Descriptions`, Fraction, isLabel, `Precursor Area`, scenario) %>%
      mutate(`Precursor Area` = na_if(`Precursor Area`, 0)) %>%
      mutate(isLabel = as.character(isLabel),
             Fraction = as.integer(Fraction))

      })

  # Report on use case
  # output$useCase <- renderText({
  #   print(input$processRaw)
  # })

  output$useCase <- renderText({
    if (input$processRaw != 0) {
      "Using raw data file. Proceed to part 1."
    } else if (input$processNorm != 0) {
      "Using normalized vales file. Proceed to part 2."
    } else {
      "No file uploaded and no example file chosed. Please upload a file or click on one of the process buttons in the above tabs."
    }
  })

  output$NormInputTest_0 <- renderText({

    if(input$processNorm != 0) {
      # if(input$processRaw != 0) {
        # names(normalized_full())
      is.object(compiledNorm_import())
      # names(peptides())
      # print("hello, there.")
    } else {
      print("no norm file")
    }

  })


  # To test a widgets value
  output$test_widget_value <- renderText({names(peptides())[2]})

  ########################
  # summary tab server side
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
  # filter tab server side
  ########################

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
      Sys.sleep(1)
      incProgress(1/n, detail = "Filtering...")
      peptides_filtered <- toFilter(peptides(), rank = as.numeric(input$rank), cl = input$checkGroup)

      # Step 2:
      Sys.sleep(1)
      incProgress(1/n, detail = "Formating Modifications and labels...")
      peptides_filtered <- splitModLab(peptides_filtered)

      # Step 3:
      Sys.sleep(1)
      incProgress(1/n, detail = "Simplifying...")
      if(input$simplify) {
        peptides_filtered <- simplifyProteins(peptides_filtered)
      }

      # Step 4:
      Sys.sleep(1)
      incProgress(1/n, detail = "Ready")
      Sys.sleep(1)

    })

    # makeEnv(peptides)
    return(peptides_filtered)
  })

  # Make the Environment...
  observeEvent(peptides_filtered(),{
    withProgress(message = 'Making environment...', value = 0, {

      n <- 6
      # Sys.sleep(3)
      incProgress(1/n, detail = "Making environment")
      makeEnv(peptides_filtered())

      #run new version of pickPeptide_new on all
      incProgress(1/n, detail = "Choosing Representative peptides")
      for (i in names(peptide_index)) {

        # assign(i, pickPeptide_new(peptide_index[[i]]), envir = peptide_index)
        assign(i, pickPeptide(peptide_index[[i]]), envir = peptide_index)

      }

      #maybe write out in how many fractions was a gived protein detected
      # count(unique(peptide_index[[listOnlyOneLabState[["onlyLabelled"]][30]]]['Fraction']))

    })
  })

  listOnlyOneLabState <- eventReactive(peptides_filtered(), {
    onlyInOneLabelState_ENV(peptide_index)
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

  # List of filtered values for selecting proteins
  output$dt <- renderUI({
    selectInput("proteinsUnion", label = "Choose protein species",
                choices = peptides_filtered()$`Protein Group Accessions`[1:10],
                multiple = FALSE)
  })

  ########################
  # analysis tab server side
  ########################

  proteinLists <- eventReactive(input$filter, {
    onlyInOneLabelState_ENV(peptide_index)
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

  # Display table for selections
  output$x1 = DT::renderDataTable(cars, server = FALSE)

  # highlight selected rows in the scatterplot
  output$x2 = renderPlot({
    s = input$x1_rows_selected
    par(mar = c(4, 4, 1, .1))
    plot(cars)
    if (length(s)) points(cars[s, , drop = FALSE], pch = 19, cex = 2)
  })

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

    allPeptidesPlot(peptide_index, protein, max(peptides()$Fraction),
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

  ###############################
  # Normalization tab server side
  ###############################

  protNormLab <- eventReactive(input$normData, {

    withProgress(message = 'Normalizing labeled...', value = 0, {
      incProgress(0.3, detail = "working...")

      names(peptide_index) %>%
        map_df(~ extractRepPeps(peptide_index[[.]],scenario = 'A', label = T)) %>%
        normalizeTable()

    })
  })

  protNormUnlab <- eventReactive(input$normData, {
    withProgress(message = 'Normalizing unlabeled...', value = 0, {
      incProgress(0.3, detail = "working...")
      names(peptide_index) %>%
        map_df(~ extractRepPeps(peptide_index[[.]],scenario = 'A', label = F)) %>%
        normalizeTable()
    })
  })

  protNormComb <- eventReactive(input$normData, {
    withProgress(message = 'Normalizing combined...', value = 0, {
      incProgress(0.3, detail = "working...")
      names(peptide_index) %>%
        map_df(~ extractRepPeps(peptide_index[[.]],scenario = 'B')) %>%
        normalizeTable()
    })
  })

  compiledExport <- reactive({
    normTableForExport(protNormLab(), protNormUnlab(), protNormComb())
    })

  compiledNorm <- reactive({
      normTableWideToLong(protNormLab(), protNormUnlab(), protNormComb())
    })

  output$dl_Norm <- renderUI({
    if(is.null(compiledExport())) {
      return(0)
    } else {
      downloadLink("downloadNormData", "Download Normalized Values")
    }
  })

  output$downloadNormData <- downloadHandler(
    filename = function() {
      "data.txt"
    },
    content = function(file) {
      export(compiledExport(), file)
    }
  )

  output$NormTest <- renderText({
    if(is.null(input$normData)) {
      return()
    } else {
      names(compiledNorm())
    }

  })

  ###############################
  # proteinNormViz tab server side
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
  # heatMaps tab server side
  ###############################

  # Just some test output for me
  output$HeatTest <- renderText({

    input$heatMapFile$datapath

  })

  # List of filtered values for selecting proteins
  output$dt_3 <- renderUI({
    if (input$renameProteinsHeatMap) {
      textInput("heatMapGroupNameColumn", label = "Column with new names", value = "Protein Group Accessions")
    }
  })

  # Make reactive plot object
  heatMapPlotObject <- reactive({
    req(input$heatMapFile)

    tryCatch(
      {
        df <- read_tsv(input$heatMapFile$datapath)
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        stop(safeError(e))
      }
    )

    if (is.null(input$heatMapGroupNameColumn)) {
      groupHeatMap(compiledNorm_plot()[compiledNorm_plot()$scenario == "B",], df, input$heatMapGroupName,
                   titleAlign = "center",
                   grid = F, colNumber = input$showSamplesHeatMap,
                   labelled = input$labelledName,
                   unlabelled = input$unlabelledName)
    } else {
      groupHeatMap(compiledNorm_plot()[compiledNorm_plot()$scenario == "B",], df, input$heatMapGroupName,
                   titleAlign = "center",
                   newNamesCol = input$heatMapGroupNameColumn,
                   grid = F, colNumber = input$showSamplesHeatMap,
                   labelled = input$labelledName,
                   unlabelled = input$unlabelledName)
    }

  })

  # Send heatmap to GUI
  output$heatMapPlot = renderPlot({
    heatMapPlotObject()
  })

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
  # coMigration tab server side
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
  coMig_1_plot <- reactive({
    req(input$groupData_coMig1)

    oneGroupTwoLabelsCoMigration(compiledNorm_plot(), as.numeric(max(compiledNorm_plot()$Fraction)),
                                                         # groupData = groupDataVector,
                                                         # groupName = groupName,

                                                         groupData = strsplit(x = input$groupData_coMig1, split = "\\s")[[1]],
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
  # cluster tab server side
  ###############################

  forClustering <- reactive({
    compiledNorm_plot() %>%
    # as_tibble() %>%
      filter(scenario == "A") %>%
      select(-scenario) %>%
      mutate(`Precursor Area` = replace_na(`Precursor Area`, 0)) %>%
      spread(Fraction, `Precursor Area`)
  })

  labelledTable <- reactive({forClustering()[forClustering()$isLabel==TRUE,]})
  unlabelledTable <- reactive({forClustering()[forClustering()$isLabel==FALSE,]})

  #create distance matrix
  labDist <- reactive({
    withProgress(message = "Calculating distance matrix...", value = 0, {
      incProgress(0.5, detail = "Working...")
      makeDist(t(select(labelledTable(),-c(1:3))), centered = input$distCentered)
    })
    })

  unlabDist <- reactive({
    withProgress(message = "Calculating distance matrix...", value = 0, {
      incProgress(0.5, detail = "Working...")
      makeDist(t(select(unlabelledTable(),-c(1:3))), centered = input$distCentered)
    })
    })

  # Generate slider for cutoff
  output$UI_distCutoff <- renderUI({
    req(labDist())
    req(unlabDist())

    sliderInput("distCutoff", "Cutoff for distance matrix",
                min = round(min(labDist(),unlabDist()),2),
                max = round(max(labDist(),unlabDist()),2),
                value = 0.8, step = 0.05)
  })

  labelledTable_clust <- reactive({assignClusters(labelledTable(), labDist(),
                                  method = input$distMethod, cutoff = input$distCutoff)})

  unlabelledTable_clust <- reactive({assignClusters(unlabelledTable(),unlabDist(),
                                    method = input$distMethod, cutoff = input$distCutoff)})

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
