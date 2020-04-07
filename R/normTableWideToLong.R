#' Convert Normalized Dataframe To Long format
#'
#' This is a convenient function for plotting
#'
#' @param labTab a dataframe
#' @param unlabTab a dataframe
#' @param comboTab a dataframe
#'
#' @return a dataframe
#' @export
normTableWideToLong <- function(labTab, unlabTab, comboTab) {

  # convert normalized data frames from wide to long format - will be used for plotting
  labTab %>%
    gather("condition", `Precursor Area`, -c(`Protein Group Accessions`, `Protein Descriptions`)) %>%
    separate("condition", c('Fraction', 'isLabel'), sep = '_') %>%
    mutate(scenario = "A") -> labTab

  unlabTab %>%
    gather("condition", `Precursor Area`, -c(`Protein Group Accessions`, `Protein Descriptions`)) %>%
    separate("condition", c('Fraction', 'isLabel'), sep = '_') %>%
    mutate(scenario = "A") -> unlabTab

  comboTab %>%
    gather("condition", `Precursor Area`, -c(`Protein Group Accessions`, `Protein Descriptions`)) %>%
    separate("condition", c('Fraction', 'isLabel'), sep = '_') %>%
    mutate(scenario = "B") -> comboTab

  # .data$Fraction <- as.numeric(as.character(.data$Fraction))

  compiled  <- labTab %>%
    bind_rows(unlabTab) %>%
    bind_rows(comboTab)

  compiled$Fraction <- as.integer(as.character(compiled$Fraction))

  compiled$`Precursor Area`[compiled$`Precursor Area` == 0] <- NA

  return(compiled)
}
