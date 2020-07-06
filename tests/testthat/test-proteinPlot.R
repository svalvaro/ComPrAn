test_that("output is a ggplot object", {
    
    ## Use example normalised proteins file
    inputFile <- system.file("extdata", "dataNormProts.txt", package = "ComPrAn")
    forAnalysis <- protImportForAnalysis(inputFile)
    protein <- "P52815"
    max_frac <- 23
    output <- proteinPlot(forAnalysis[forAnalysis$scenario == "B",], protein, max_frac)
    
    expect_s3_class(output, "ggplot")
})
