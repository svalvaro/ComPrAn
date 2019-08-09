#' Make a distmatrix
#'
#' This function modifies the input DF data frame by performing clustering on distObject
#' and adding a 'cluster' column to the DF
#'
#' @param DF DF: data frame, contains columns: `Protein Group Accessions` character `Protein Descriptions` character
#'                isLabel character ('TRUE'/'FALSE')
#'                               columns 1 to n, numeric, n is the total number of fractions/slices, each of this columns
#'                                                          contains `Precursor Area` values in a given fraction(columns) for a protein(rows)
#' @param distObject input data
#' @param method character One of 'average', 'single' or 'complete' (default), specifies the linkage method to be used inside R hclust() function
#' @param cutoff numeric, specifies the h value in R cutree() function, height at which to 'cut the tree', everything with distance below this value is assigned into same cluster everything with larger distance is in a different cluster extreme possible values are 0 to 2 (might not be reached for all data sets)
#'
#' @return dataframe
#' @export
assignClusters <- function(DF, distObject, method = 'complete', cutoff = 0.5){

  #creates a new column in a DF containing cluster numbers
  hcData <- hclust(distObject, method = method)
  clusters <- cutree(hcData, h = cutoff)
  return(mutate(DF,cluster = clusters))
}

