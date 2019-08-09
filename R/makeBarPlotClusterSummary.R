#' Title
#'
#' @param df data frame, contains columns:
#'              `Protein Group Accessions` character
#'              `Protein Descriptions` character
#'               isLabel character ('TRUE'/'FALSE')
#'               columns 1 to n, numeric, n is the total number of fractions/slices, each of this columns
#'                           contains `Precursor Area` values in a given fraction(columns) for a protein(rows)
#'               cluster integer
#' @param name character, specifies the name of the sample
#'
#' @return plot
#' @export
makeBarPlotClusterSummary <- function(df, name = 'sample 1') {
  df %>%
    group_by(cluster) %>%
    mutate(n = n()) %>%
    ggplot(aes(reorder(cluster,-n))) +
    geom_bar()+
    theme_minimal()+
    labs( x= 'Cluster number', y ='Number of proteins in cluster', title=name)
}
