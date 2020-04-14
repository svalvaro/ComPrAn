options(shiny.maxRequestSize = 150*1024^2)

plotBkg <- "#ecf0f5"
rbCol <- brewer.pal(3, "Set1")[1:2]

col_vector_peptides <- c('#ffc125', '#a020f0')
names(col_vector_peptides) <- c('TRUE', 'FALSE')

col_vector_proteins <- c('#ff9d2e', '#07b58a')
names(col_vector_proteins) <- c('TRUE', 'FALSE')

# Read in data
dataFile <- system.file("extdata", "data.txt", package = "ComPrAn")
peptides <- fread(dataFile)