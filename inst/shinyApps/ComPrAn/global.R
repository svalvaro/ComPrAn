options(shiny.maxRequestSize = 150*1024^2)

plotBkg <- "#ecf0f5"
rbCol <- brewer.pal(3, "Set1")[1:2]

# Read in data
dataFile <- system.file("extdata", "data.txt", package = "ComPrAn")
peptides <- fread(dataFile)