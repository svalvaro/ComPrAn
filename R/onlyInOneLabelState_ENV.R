#' Report Proteins Present In Only One Label State
#'
#' This function returns NAMES of proteins present in only labelled/only unlabelld or both label states
#'
#' @param .data An environment containing dataframes
#'
#' @return a list with 3 items, each item is a vector containing names belonging to one of 3 groups
#' @export
#' 
#' @examples
#' ##Use example peptide data set, read in and clean data
#' inputFile <- system.file("extdata", "data.txt", package = "ComPrAn")
#' peptides <- cleanData(data.table::fread(inputFile), fCol = "Search ID")
#' ## separate chemical modifications and labelling into separate columns
#' peptides <- splitModLab(peptides) 
#' ## remove unneccessary columns, simplify rows
#' peptides <- simplifyProteins(peptides) 
#' ##make environment
#' peptide_index <- makeEnv(peptides)
#' ## Pick representative peptide for each protein for both scenarios
#' for(i in names(peptide_index)){assign(i,pickPeptide(peptide_index[[i]]),envir=peptide_index)}
#' ## Create a list of proteins present in both/only in one label state
#' listOnlyOneLabState <- onlyInOneLabelState_ENV(peptide_index)
onlyInOneLabelState_ENV <- function(.data){
  #This function returns NAMES of proteins present in only labelled/only unlabelld or both label states
  #As an input this function takes an environment
  #Output is a list with 3 items, each item is a vector containing names belonging to one of 3 groups

  results <- sapply(.data, function(x)
    if(all(!x$isLabel)){
      'onlyUnlabelled'
    } else if(all(x$isLabel)){
      'onlyLabelled'
      #add to vecA
    } else {
      'both'
    })

  results <- list(onlyLabelled = names(results[results == "onlyLabelled"]),
                  onlyUnlabelled = names(results[results == "onlyUnlabelled"]),
                  both = names(results[results == "both"]))

  return(results)

}
