#' Covert clustered tables into format for export
#'
#' @param labClustTable output: data frame containing columns:
#'              `Protein Group Accessions` character
#'              `Protein Descriptions` character
#'              `Cluster number - unlabeled` integer
#'              `Cluster number - labeled` integer
#' @param unlabClustTable labClustTable, unlabClustTable: data frames, 
#' contain columns:
#'              `Protein Group Accessions` character
#'              `Protein Descriptions` character
#'               isLabel character ('TRUE'/'FALSE') - here in one data frame all
#'                are TRUE in second one all are FALSE
#'               columns 1 to n, numeric, n is the total number of 
#'               fractions/slices, each of this columns
#'              contains `Precursor Area` values in a given fraction(columns)
#'               for a protein(rows) cluster integer
#'
#' @importFrom tidyr spread
#'
#' @return dataframe
#' @export
#'
#' @examples 
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
#' #make table of cluster assginment
#' tableClusterExport <- exportClusterAssignments(labTab_clust,unlabTab_clust)
#' 
exportClusterAssignments <- function(labClustTable, unlabClustTable){
    clustTable <- rbind(labClustTable, unlabClustTable)
    clustTable  %>%
        select(`Protein Group Accessions`, 
                `Protein Descriptions`,
                isLabel, cluster) %>%
        spread(isLabel,cluster) %>%
        rename(`Cluster number - labeled` = 'TRUE',
                `Cluster number - unlabeled` = 'FALSE') ->clustTable
    return(clustTable)
}
