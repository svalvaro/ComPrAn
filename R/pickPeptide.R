#' Select Top Peptide For Various Scenarios
#'
#' This function selects a single unique peptide to represent each Protein Group Accession
#' There are 2 ways of selecting peptides, both are perform as they are needed for different tasks later on
#' 1) Scenario A: select peptide occuring in most fractions, do this individually for labelled/unlabelled
#' (max value for any peptide is equal to number of fractions) in case of ties, pick
#' peptide whith highest Precursor Area in any fraction
#' 2) Scenario B: select peptide occuring in most fractions counting both label states together
#' (max value for any peptide is equal to twice the number of fractions)
#' in case of ties, pick peptide with highest Precursor Area in any fraction.
#' Representative peptide in Scenario B is picked only for proteins that have
#' at least one shared peptide between label states.
#'
#' @param .data a dataframe
#'
#' @return a dataframe
#' @export
#' 
#' @examples 
#' 
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
#'
pickPeptide <- function(.data) {

  ##############scenarioA
  .data %>%
    dplyr::group_by(UniqueCombinedID_A, isLabel)%>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(isLabel) %>%
    dplyr::mutate(repPepA = ifelse(n == max(n), TRUE, FALSE)) %>%
    dplyr::mutate(maxArea = max(`Precursor Area`[repPepA])) %>%
    dplyr::mutate(maxAreaPeptide = UniqueCombinedID_A[(`Precursor Area` == maxArea & repPepA)][1]) %>%
    dplyr::mutate(repPepA = ifelse(repPepA & UniqueCombinedID_A == maxAreaPeptide, TRUE, FALSE)) %>%
    dplyr::select(-n, -maxArea, -maxAreaPeptide) -> .data


  .data %>%
    dplyr::group_by(UniqueCombinedID_B) %>%
    dplyr::mutate(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(UniqueCombinedID_A) %>%
    dplyr::mutate(n_2 = sum(n[isLabel][1], n[!isLabel][1], na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(maxN = ifelse(identical(n_2[n_2 != n], integer(0)),9999,max(n_2[n_2 != n]))) %>%
    dplyr::mutate(repPepB = ifelse(n_2 == maxN, TRUE, FALSE))-> .data

  if(all(.data$isLabel)|all(!.data$isLabel)){
    .data$repPepB = FALSE
  }

  if(!all(!.data$repPepB)){
    .data %>%
      dplyr::mutate(maxArea = UniqueCombinedID_A[`Precursor Area` == max(`Precursor Area`[repPepB])& repPepB][1]) %>%
      dplyr::mutate(repPepB = ifelse(repPepB & UniqueCombinedID_A == maxArea, TRUE, FALSE)) %>%
      dplyr::select(-n, -n_2, -maxN, -maxArea) -> .data
  }else{
    .data %>%
      dplyr::select(-n, -n_2, -maxN) -> .data
  }
  return(.data)
}
