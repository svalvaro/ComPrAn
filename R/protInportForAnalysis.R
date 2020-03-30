#' Modify import protein data
#' This function converts imported protien table into format compatible with downstream analysis
#' @param .data - data frame, contains columns:
#' @param "Protein Group Accessions" -  character/factor
#' @param "Protein Descriptions" - character
#' @param "scenario" - character/factor
#' @param "label" - logical
#' @param columns "1" to "n" - numeric

protInportForAnalysis <- function(.data){
  .data %>%  gather(Fraction, `Precursor Area`, -c(`Protein Group Accessions`, `Protein Descriptions`, scenario, label)) %>%
    rename("isLabel" = "label") %>%
    select(`Protein Group Accessions`, `Protein Descriptions`, Fraction, isLabel, `Precursor Area`, scenario) %>% 
    mutate(`Precursor Area` = na_if(`Precursor Area`, 0)) %>%
    mutate(isLabel = as.character(isLabel),
           Fraction = as.integer(Fraction)) -> .data
  return(.data)
}