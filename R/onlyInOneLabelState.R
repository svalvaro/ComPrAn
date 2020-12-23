#' Report Proteins Present In Only One Label State
#'
#' This function returns NAMES of proteins present in only labelled/only 
#' unlabelld or both label states
#'
#' @param .data An environment containing dataframes
#'
#' @return a list with 3 items, each item is a vector containing names
#' belonging to one of 3 groups
#' @export
#' 
#' @examples 
#' 
#' ##Use example peptide data set, read in and clean data
#' inputFile <- system.file("extData", "data.txt", package = "ComPrAn")
#' peptides <- peptideImport(inputFile)
#' peptides <- cleanData(peptides, fCol = "Search ID")
#' ## separate chemical modifications and labelling into separate columns
#' peptides <- splitModLab(peptides) 
#' ## remove unneccessary columns, simplify rows
#' peptides <- simplifyProteins(peptides) 
#' ## Pick representative peptide for each protein for both scenarios
#' peptide_index <- pickPeptide(peptides)
#' ## extract list of names of proteins present in one/both samples
#' oneStateList <- onlyInOneLabelState(peptide_index)
onlyInOneLabelState <- function(.data){
    results <- vapply(.data, function(x)
        if(all(!x$isLabel)){
            'onlyUnlabelled'
        } else if(all(x$isLabel)){
            'onlyLabelled'
            #add to vecA
        } else {
            'both'
        },"a")
    
    results <- list(onlyLabelled = names(results[results == "onlyLabelled"]),
                    onlyUnlabelled =names(results[results == "onlyUnlabelled"]),
                    both = names(results[results == "both"]))
    return(results)
}