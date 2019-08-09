# Test script

#   Install Package:           'Cmd + Shift + B'
load_all()

# Read in data
peptides <- data.table::fread("../complexomics_app/data/data.txt")

peptides <- cleanData(peptides, fCol = "Search ID")
peptides <- toFilter(peptides, rank = 1)
peptides <- splitModLab(peptides)
peptides <- simplifyProteins(peptides)

####### PART 2 ########
makeEnv(peptides)

####### PART 3 ########
#run pickPeptide on all
for (i in names(peptide_index)) {
  assign(i, pickPeptide(peptide_index[[i]]), envir = peptide_index)
}

listOnlyOneLabState <- onlyInOneLabelState_ENV(peptide_index)

names(peptide_index) %>%
  map_df(~ extractRepPeps(peptide_index[[.]], scenario = 'A', label = T))  %>%
  normalizeTable() -> protNormLab

names(peptide_index) %>%
  map_df(~ extractRepPeps(peptide_index[[.]], scenario = 'A', label = F))  %>%
  normalizeTable() -> protNormUnlab

names(peptide_index) %>%
  map_df(~ extractRepPeps(peptide_index[[.]], scenario = 'B')) %>%
  normalizeTable() -> protNormComb

forExport <- normTableForExport(protNormLab, protNormUnlab, protNormComb)
forAnalysis <- normTableWideToLong(protNormLab, protNormUnlab, protNormComb)

col_vector_peptides <- c('#ffc125', '#a020f0')
names(col_vector_peptides) <- c('TRUE', 'FALSE')

col_vector_proteins <- c('#ff9d2e', '#07b58a')
names(col_vector_proteins) <- c('TRUE', 'FALSE')


proteinPlot(forAnalysis[forAnalysis$scenario == "B",], "Q9Y2R5", 23)


groupData <- read_tsv("../complexomics_app/data/mitoRiboLSU_uniprotID_filtered.txt")
groupName <- 'mtLSU' # input panel for user to type in group name or by default name of the group file
colNumber <- 2

groupHeatMap(dataFrame = forAnalysis[forAnalysis$scenario == "B",], groupData, groupName,
             titleAlign = "center",
             #newNamesCol = "Gene name",  # if newNamesCol is not specified, names stay default = 'Protein Group Accessions'
             grid = F, colNumber = 2)

groupHeatMap(forAnalysis[forAnalysis$scenario == "B",], groupData, groupName,
             titleAlign = "center",
             newNamesCol = "Gene name",
             grid = F, colNumber = 2)

# Co-migration plots

groupData <- read_tsv("../complexomics_app/data/mitoRiboLSU_uniprotID_filtered.txt")
groupName <- 'mtLSU' # input panel for user to type in group name or by default name of the group file

groupDataVector <- groupData$`Protein Group Accessions` # this is the neccessary input format for group

max_frac <- 23 # number of fractions

oneGroupTwoLabelsCoMigration(forAnalysis,
                             max_frac = max_frac,
                             groupDataVector,groupName,
                             meanLine = T) # example plot


group1DataVector <- groupDataVector # from above
group1Name <- groupName # from above

group2Data <- read_tsv("../complexomics_app/data/mitoRiboSSU_uniprotID_filtered.txt")
group2DataVector <- group2Data$`Protein Group Accessions`
group2Name <- 'mtSSU'


twoGroupsWithinLabelCoMigration(dataFrame = forAnalysis, max_frac = max_frac,
                                group1Data = group1DataVector, group1Name = group1Name,
                                group2Data = group2DataVector, group2Name = group2Name,
                                grid = T, titleAlign = 'center',showTitle = T,
                                jitterPoints = 0.3, alphaValue = 0.3, meanLine = F, medianLine = T) # example plot

twoGroupsWithinLabelCoMigration(forAnalysis,
                                max_frac = max_frac) # return error message when no groupData provided


# Clustering

forAnalysis %>%
  as_tibble() %>%
  filter(scenario == "A") %>%
  select(-scenario) %>%
  mutate(`Precursor Area` = replace_na(`Precursor Area`, 0)) %>%
  spread(Fraction, `Precursor Area`) -> forClustering

forClustering[is.na(forClustering)] <- 0




forAnalysis[forAnalysis$scenario == "A",] %>%
  select(-scenario) %>%
  spread(Fraction, `Precursor Area`) -> forClustering

forClustering[is.na(forClustering)] <- 0

labelledTable <- forClustering[forClustering$isLabel==TRUE,]
unlabelledTable <- forClustering[forClustering$isLabel==FALSE,]

# Create distance matrix
labDist <- makeDist(t(select(labelledTable,-c(1:3))), centered = T)
unlabDist <- makeDist(t(select(unlabelledTable,-c(1:3))), centered = T)

# # Assign clusters to DFs
# min(labDist,unlabDist)
# max(labDist,unlabDist)


#cutoff - can this be slidebar from min(labDist,unlabDist) to max(labDist,unlabDist)?
labelledTable_clust <- assignClusters(labelledTable, labDist,method = 'average', cutoff = 0.8)
unlabelledTable_clust <- assignClusters(unlabelledTable,unlabDist ,method = 'average', cutoff = 0.8)

#make bar plots summarizing numbers of proteins per cluster
makeBarPlotClusterSummary(labelledTable_clust)
makeBarPlotClusterSummary(unlabelledTable_clust)

#create table for export

tableForClusterExport <- exportClusterAssignments(labelledTable_clust,unlabelledTable_clust)
