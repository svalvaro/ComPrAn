#' Create a data frames with cluster assignment
#'
#' This function creates a data frame with column specifying clusters assigned
#' ot each protein using the table and distance matrix produced by clusterComp()
#' function.
#'
#' @param .listDf list of data frames produced by clusterComp() function
#' @param sample which of the two samples you want to apply the function 
#' to (labeled/unlabeled).
#' @param method character, One of 'average', 'single' or 'complete' (default),
#'  specifies the linkage method to be used inside R hclust() function
#' @param cutoff numeric, specifies the h value in R cutree() function, 
#' height at which to 'cut the tree', everything with distance below this 
#' value is assigned into same cluster everything with larger distance is 
#' in a different cluster extreme possible values are 0 to 2 (might not be 
#' reached for all data sets)
#'
#' @importFrom stats hclust cutree
#'
#' @return dataframe
#' @export
#' 
#' @seealso \code{\link{clusterComp}}
#' @examples
#' 
#' ##Use example normalised proteins file
#' inputFile <- system.file("extData", "dataNormProts.txt", package = "ComPrAn")
#' #read file in and change structure of table to required format
#' forAnalysis <- protImportForAnalysis(inputFile)
#' # create components necessary for clustering
#' clusteringDF <- clusterComp(forAnalysis,scenar = "A", PearsCor = "centered")
#' #assign clusters
#' labTab_clust <- assignClusters(.listDf = clusteringDF,sample = "labeled",
#' method = 'complete', cutoff = 0.5)
#' unlabTab_clust <- assignClusters(.listDf = clusteringDF,sample = "unlabeled",
#'                                method = 'complete', cutoff = 0.5)
assignClusters <- function(.listDf, sample , method = 'complete', cutoff = 0.5){
    if(sample == "labelled"|sample == "labeled"){
        DF <- .listDf[["labTable"]]
        distObject <- .listDf[["labDistM"]]
    } else if (sample == "unlabelled"|sample == "unlabeled"){
        DF <- .listDf[["unlabTable"]]
        distObject <- .listDf[["unlabDistM"]]
    } else{
        stop("Valid values for \"sample\" are \"labeled\" or \"unlabeled\".")
    }
    
    if(!method %in% c("complete","average","single"))
        stop("Valid values for \"method\" are \"complete\", \"average\"
                or \"single\".")
    #creates a new column in a DF containing cluster numbers
    hcData <- stats::hclust(distObject, method = method)
    clusters <- stats::cutree(hcData, h = cutoff)
    return(mutate(DF,cluster = clusters))
}