#' Create components necessary for clustering
#'
#' Reformat the table for the one neccessary for assignClusters function.
#' Calculate the distance matirx using selected variant of correlation.
#' 
#' @param .df data frame, table of normalised protein values
#' @param scenar character, scenario intended for clustering, either "A" or "B"
#' @param PearsCor character, pearsons correlation variant (centered/uncentered)
#' 
#' @importFrom tidyr replace_na
#' @return list of data frames
#' @export
#' 
#' @examples
#' 
#' ##Use example normalised proteins file
#' inputFile <- system.file("extdata", "dataNormProts.txt", package = "ComPrAn")
#' #read file in and change structure of table to required format
#' forAnalysis <- protInportForAnalysis(data.table::fread(inputFile))
#' # create components necessary for clustering
#' clusteringDF <- clusterComp(forAnalysis,scenar = "A", PearsCor = "centered")
#' 
clusterComp <- function(.df, scenar = "A", PearsCor = "centered"){
    .df %>% 
        filter(scenario == scenar ) %>% 
        select(-scenario) %>%
        mutate(`Precursor Area` = replace_na(`Precursor Area`, 0)) %>% 
        spread(Fraction, `Precursor Area`) -> .df
    
    labTab <- .df[.df$isLabel==TRUE,]
    unlabTab <- .df[.df$isLabel==FALSE,]
    
    if (PearsCor == "centered"|PearsCor == "centred"){
        labDist <- makeDist(t(select(labTab,-c(1:3))), 
                                     centered = TRUE)
        unlabDist <- makeDist(t(select(unlabTab,-c(1:3))), 
                                       centered = TRUE)
    } else if(PearsCor == "uncentered"|PearsCor == "uncentred"){
        labDist <- makeDist(t(select(labTab,-c(1:3))), 
                                     centered = FALSE)
        unlabDist <-makeDist(t(select(unlabTab,-c(1:3))), 
                                      centered = FALSE)
    } else(
        stop("Valid values for \"PearsCor\" are \"centered\" or 
             \"uncentered\".")
    )
    return(list(labTable = labTab,
                unlabTable = unlabTab,
                labDistM = labDist,
                unlabDistM =unlabDist))
}