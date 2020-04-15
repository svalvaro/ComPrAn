#' Report Proteins Present In Only One Label State
#'
#' This function returns NAMES of proteins present in only labelled/only unlabelld or both label states
#'
#' @param .data An environment containing dataframes
#'
#' @return a list with 3 items, each item is a vector containing names belonging to one of 3 groups
#' @export
#' 
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
