#' Title
#'
#' @param df data frame, contains columns:
#'              `Protein Group Accessions` character
#'              `Protein Descriptions` character
#'               isLabel character ('TRUE'/'FALSE')
#'               columns 1 to n, numeric, n is the total number of 
#'               fractions/slices, each of this columns
#'               contains `Precursor Area` values in a given fraction(columns) 
#'               for a protein(rows)
#'               cluster integer
#' @param name character, specifies the name of the sample
#'
#' @importFrom stats reorder
#'
#' @return plot
#' @export
#' 
#' @examples
#' 
#' ##Use example normalised proteins file
#' inputFile <- system.file("extdata", "dataNormProts.txt", package = "ComPrAn")
#' #read file in and change structure of table to required format
#' forAnalysis <- protImportForAnalysis(inputFile)
#' # create components necessary for clustering
#' clusteringDF <- clusterComp(forAnalysis,scenar = "A", PearsCor = "centered")
#' #assign clusters
#' labTab_clust <- assignClusters(.listDf = clusteringDF,sample = "labeled",
#' method = 'complete', cutoff = 0.5)
#' unlabTab_clust <- assignClusters(.listDf = clusteringDF,sample = "unlabeled",
#'                                method = 'complete', cutoff = 0.5)
#' #Make bar plots for labeled and unlabeled samples
#' makeBarPlotClusterSummary(labTab_clust, name = 'labeled')
#' makeBarPlotClusterSummary(unlabTab_clust, name = 'unlabeled')
makeBarPlotClusterSummary <- function(df, name = 'sample 1') {
    df %>%
        group_by(cluster) %>%
        mutate(n = n()) %>%
        ggplot(aes(reorder(cluster,-n))) +
        geom_bar()+
        theme_minimal()+
        labs(x= 'Cluster number', 
                y ='Number of proteins in cluster', 
                title=name)
}