#' Select Top Peptide For Various Scenarios
#'
#' This function selects a single unique peptide to represent each 
#' `Protein Group Accession`
#' There are 2 ways of selecting peptides, both are perform as they are needed
#' for different tasks later on.
#' \enumerate{
#' \item Scenario A: select peptide occuring in most fractions, do this 
#' individually for labelled/unlabelled (max value for any peptide is equal to 
#' number of fractions) in case of ties, pick peptide whith highest 
#' `Precursor Area` in any fraction.
#' \item Scenario B: select peptide occuring in most fractions counting both
#' label states together (max value for any peptide is equal to twice the 
#' number of fractions) in case of ties, pick peptide with highest
#' `Precursor Area` in any fraction. Representative peptide in Scenario B is 
#' picked only for proteins that have shared peptide between label states.
#' }
#' 
#' @param .data a dataframe
#'
#' @return list of data frames
#' @export
#' 
#' @examples 
#' 
#' ##Use example peptide data set, read in and clean data
#' inputFile <- system.file("extdata", "data.txt", package = "ComPrAn")
#' peptides <- peptideImport(inputFile)
#' peptides <- cleanData(peptides, fCol = "Search ID")
#' ## separate chemical modifications and labelling into separate columns
#' peptides <- splitModLab(peptides) 
#' ## remove unneccessary columns, simplify rows
#' peptides <- simplifyProteins(peptides) 
#' ## Pick representative peptide for each protein for both scenarios
#' peptide_index <- pickPeptide(peptides)
pickPeptide <- function(.data){ 
    .data %>%        ############## scenario A
        dplyr::group_by(`Protein Group Accessions`,
                        UniqueCombinedID_A,isLabel)%>%
        dplyr::mutate(n = n()) %>% dplyr::ungroup() %>%
        dplyr::group_by(`Protein Group Accessions`,isLabel) %>%
        dplyr::mutate(repPepA = ifelse(n == max(n), TRUE, FALSE)) %>%
        dplyr::mutate(maxArea = max(`Precursor Area`[repPepA])) %>%
        dplyr::mutate(maxAreaPeptide = 
            UniqueCombinedID_A[(`Precursor Area` == maxArea & repPepA)][1]) %>%
        dplyr::mutate(repPepA = ifelse(repPepA & 
                                        UniqueCombinedID_A == maxAreaPeptide,
                                        TRUE,FALSE))%>%
        dplyr::select(-n, -maxArea, -maxAreaPeptide) -> .data
    .data %>%        ############### scenarion B
        dplyr::group_by(`Protein Group Accessions`,UniqueCombinedID_B) %>%
        dplyr::mutate(n = n()) %>% dplyr::ungroup() %>%
        dplyr::group_by(`Protein Group Accessions`,UniqueCombinedID_A) %>%
        dplyr::mutate(n_2 = sum(n[isLabel][1], n[!isLabel][1],na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(`Protein Group Accessions`) %>% 
        dplyr::mutate(maxN = ifelse(identical(n_2[n_2 != n], integer(0)),9999,
                                    max(n_2[n_2 != n]))) %>%
        #add TRUE to all peptides that are present in maximum number of 
        #fractions this can happen for multiple peptides in any given protein
        #here we don't check whether peptide is present in BOTH label states
        dplyr::mutate(repPepB = ifelse(n_2 == maxN, TRUE, FALSE)) %>%
        dplyr::ungroup() %>% 
        # check whether each peptide is present in both label states, if there 
        # are peptides present only for one label state repPepB MUST be FALSE 
        # for this `Protein Group Accessions` and `UniqueCombinedID_A`
        dplyr::group_by(`Protein Group Accessions`,`UniqueCombinedID_A`) %>% 
        dplyr::mutate(repPepB = ifelse(length(unique(isLabel)) == 1,FALSE,
                                        repPepB)) %>% dplyr::ungroup() %>% 
        #when there is TRUE in any of repPepB rows:if multpile peptides selected
        #pick the one with highest precursor area, if there are more pick first 
        dplyr::group_by(`Protein Group Accessions`) %>% 
        dplyr::mutate(maxArea = ifelse(any(repPepB), 
                                UniqueCombinedID_A[`Precursor Area` == max(
                                `Precursor Area`[repPepB])& repPepB][1],NA)) %>%
        dplyr::mutate(repPepB = ifelse(repPepB & UniqueCombinedID_A == maxArea, 
                                TRUE,FALSE)) %>%  dplyr::ungroup() %>% 
        dplyr::select(-n, -n_2, -maxN, -maxArea) -> .data
    ###change data to named list(https://github.com/tidyverse/dplyr/issues/4223)
    .data %>% dplyr::group_by(`Protein Group Accessions`) %>% 
        dplyr::group_keys(`Protein Group Accessions`) %>% pull(1) -> group_names
    .data %>% dplyr::group_split(`Protein Group Accessions`) %>% 
        rlang::set_names(group_names)->.data
    return(.data)
}